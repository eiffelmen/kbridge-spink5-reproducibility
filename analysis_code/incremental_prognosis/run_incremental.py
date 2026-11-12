#!/usr/bin/env python3
"""Locked paired incremental Cox analysis; see analysis_plan.md."""
from pathlib import Path
import hashlib,json,sys,platform,warnings
import numpy as np
import pandas as pd
from scipy.stats import chi2
import scipy, lifelines, statsmodels
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from statsmodels.duration.hazard_regression import PHReg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
SEEDS = list(range(42, 52))
NBOOT = 1000

manifest = json.loads((ROOT / "source_manifest.json").read_text())
sources = {}
input_checks = []
for entry in manifest["sources"]:
    file = ROOT / entry["copy"]
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    if digest != entry["sha256"]:
        raise ValueError(f"Source SHA-256 mismatch for {entry['key']}: {file}")
    sources[entry["key"]] = file
    input_checks.append({"key": entry["key"], "sha256": digest, "status": "verified"})

if "--check-inputs" in sys.argv:
    print(json.dumps({"status": "PASS", "checks": input_checks}, indent=2))
    raise SystemExit(0)

# The original TCGA model input is row ordered but lacks sample IDs; recover them by
# exactly reconstructing its existing inclusion/design, then asserting all values match.
tcga=pd.read_csv(sources['tcga_model'],sep='\t')
a=pd.read_csv(sources['tcga_patient_join'],sep='\t')
req=['time','event','SPINK5_per_SD','age10','clinical_stage_num','male','ever_smoked','site_group']
a=a[a.event.notna() & a.time.notna() & a.time.gt(0)].dropna(subset=req).copy()
z=a[['time','event','SPINK5_per_SD','age10','male','ever_smoked']].copy()
for s in (2,3,4): z[f'stage_{s}_vs_1']=(a.clinical_stage_num==s).astype(float)
for s in ['larynx','oropharynx','hypopharynx']: z[f'site_{s}_vs_oral_cavity']=(a.site_group==s).astype(float)
assert list(z.columns)==list(tcga.columns)
assert np.allclose(z.to_numpy(float),tcga.to_numpy(float),rtol=1e-10,atol=1e-10)
tcga.insert(0,'patient_id',a.sampleID.to_numpy())
tcga.rename(columns={'SPINK5_per_SD':'SPINK5'},inplace=True)
ext=pd.read_csv(sources['external_model']).rename(columns={'geo':'patient_id','os_native':'time','death':'event','SPINK5_z':'SPINK5'})
assert len(tcga)==491 and tcga.event.sum()==210 and tcga.patient_id.is_unique
assert len(ext)==233 and ext.event.sum()==73 and ext.patient_id.is_unique
for label,d in [('TCGA',tcga),('GSE65858',ext)]: d.to_csv(INPUT/f'{label}_verified_analysis_input.tsv',sep='\t',index=False)
configs=[('TCGA',tcga,False),('GSE65858',ext,False),('GSE65858_sex_stratified',ext,True)]
allfull=[];allfullcoef=[];allfolds=[];allpred=[];allreps=[];allcvcoef=[];allscaling=[];allmembers=[];allboot=[];allsum=[];allchecks=[]


def compare_pairs(d,risk,stratified):
    n=len(d); i,j=np.triu_indices(n,1)
    t=d.time.to_numpy();e=d.event.to_numpy()
    good=((t[i]<t[j])&(e[i]==1)) | ((t[j]<t[i])&(e[j]==1)) | ((t[i]==t[j])&(e[i]!=e[j]))
    if stratified:
        sex=d.male.to_numpy();good &= sex[i]==sex[j]
    i,j=i[good],j[good]
    earlier_i=(t[i]<t[j]) | ((t[i]==t[j])&(e[i]==1))
    early=np.where(earlier_i,i,j);late=np.where(earlier_i,j,i)
    score=np.where(risk[early]>risk[late],1.,np.where(risk[early]==risk[late],.5,0.))
    return i,j,score


def fit(df,features,stratified):
    cols=['time','event']+features+(['male'] if stratified else [])
    m=CoxPHFitter(penalizer=0)
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        m.fit(df[cols].astype(float),'time','event',strata=['male'] if stratified else None,
              fit_options={'precision':1e-9,'r_precision':1e-10,'max_steps':200})
    return m,'; '.join(str(w.message) for w in ws)

for label,data,stratified in configs:
    print('START',label,flush=True)
    data=data.reset_index(drop=True).copy()
    data['case_index']=np.arange(len(data))
    allfeatures=[c for c in data if c not in ['patient_id','time','event','case_index']]
    if stratified: allfeatures.remove('male')
    baseline=[c for c in allfeatures if c!='SPINK5']
    fullmodels={}
    for typ,features in [('clinical',baseline),('clinical_SPINK5',allfeatures)]:
        m,warn=fit(data,features,stratified); fullmodels[typ]=m
        sm=PHReg(data.time,data[features],status=data.event,ties='efron',strata=data.male if stratified else None).fit(disp=False)
        diffcoef=float(np.max(np.abs(m.params_.reindex(features).to_numpy()-sm.params)))
        difflog=float(abs(m.log_likelihood_-sm.llf))
        assert diffcoef<1e-4 and difflog<1e-6,(label,typ,diffcoef,difflog)
        allchecks.append({'cohort':label,'model':typ,'max_abs_coefficient_difference':diffcoef,'abs_loglik_difference':difflog})
        allfull.append({'cohort':label,'model':typ,'n':len(data),'deaths':int(data.event.sum()),'n_coefficients':len(features),'loglik':m.log_likelihood_,'loglik_statsmodels':sm.llf,'warnings':warn,'features':';'.join(features)})
        cf=m.summary.reset_index();cf.insert(0,'cohort',label);cf.insert(1,'model',typ);allfullcoef.append(cf)
    lr=2*(fullmodels['clinical_SPINK5'].log_likelihood_-fullmodels['clinical'].log_likelihood_)
    assert lr>=-1e-8
    # Each repetition stores paired concordance contributions from only same-fold pairs.
    caches=[]
    for repeat,seed in enumerate(SEEDS):
        rng=np.random.default_rng(seed); fold=np.full(len(data),-1,dtype=int)
        for event in [0,1]:
            ids=np.flatnonzero(data.event.to_numpy()==event);rng.shuffle(ids)
            for f,part in enumerate(np.array_split(ids,5)): fold[part]=f
        assert np.all(fold>=0)
        rep_i=[];rep_j=[];rep_c=[];rep_e=[]
        for f in range(5):
            train=data.loc[fold!=f].copy();test=data.loc[fold==f].copy()
            for c in ['age10','SPINK5']:
                mu=float(train[c].mean());sd=float(train[c].std(ddof=0));assert sd>0
                train[c]=(train[c]-mu)/sd;test[c]=(test[c]-mu)/sd
                allscaling.append({'cohort':label,'repeat':repeat,'seed':seed,'fold':f,'feature':c,'training_mean':mu,'training_population_sd':sd})
            dropped=[c for c in allfeatures if train[c].nunique()<=1]
            pred={};foldwarnings={}
            for typ,features0 in [('clinical',baseline),('clinical_SPINK5',allfeatures)]:
                features=[c for c in features0 if c not in dropped]
                m,warn=fit(train,features,stratified);foldwarnings[typ]=warn
                pred[typ]=m.predict_log_partial_hazard(test).to_numpy()
                for term,value in m.params_.items():allcvcoef.append({'cohort':label,'repeat':repeat,'seed':seed,'fold':f,'model':typ,'feature':term,'coefficient':value})
            i,j,cs=compare_pairs(test,pred['clinical'],stratified)
            ie,je,es=compare_pairs(test,pred['clinical_SPINK5'],stratified)
            assert len(cs)>0
            assert np.array_equal(i,ie) and np.array_equal(j,je)
            # Independent concordance implementation validation, matched target.
            for typ,sc in [('clinical',cs),('clinical_SPINK5',es)]:
                if stratified:
                    numer=denom=0.
                    for sex in [0,1]:
                        mask=test.male.eq(sex).to_numpy();s=test.loc[mask]
                        ii,jj,ss=compare_pairs(s,pred[typ][mask],False)
                        if len(ss)==0: continue
                        lc=concordance_index(s.time,-pred[typ][mask],s.event)
                        assert np.isclose(lc,ss.mean(),atol=1e-12)
                        numer+=lc*len(ss);denom+=len(ss)
                    checked=numer/denom
                else: checked=concordance_index(test.time,-pred[typ],test.event)
                assert np.isclose(checked,sc.mean(),atol=1e-12)
            gids=test.case_index.to_numpy()
            rep_i.append(gids[i]);rep_j.append(gids[j]);rep_c.append(cs);rep_e.append(es)
            allfolds.append({'cohort':label,'repeat':repeat,'seed':seed,'fold':f,'train_n':len(train),'test_n':len(test),'test_deaths':int(test.event.sum()),'comparable_pairs':len(cs),'clinical_correct_plus_half_ties':float(cs.sum()),'extended_correct_plus_half_ties':float(es.sum()),'C_clinical':float(cs.mean()),'C_extended':float(es.mean()),'delta_C':float(es.mean()-cs.mean()),'dropped_training_constant_features':';'.join(dropped),'warnings_clinical':foldwarnings['clinical'],'warnings_extended':foldwarnings['clinical_SPINK5']})
            for k,(idx,row) in enumerate(test.iterrows()):
                allpred.append({'cohort':label,'repeat':repeat,'seed':seed,'fold':f,'case_index':int(row.case_index),'patient_id':row.patient_id,'time':row.time,'event':int(row.event),'male':int(row.male),'risk_clinical':float(pred['clinical'][k]),'risk_extended':float(pred['clinical_SPINK5'][k])})
        for k,row in data.iterrows():allmembers.append({'cohort':label,'repeat':repeat,'seed':seed,'patient_id':row.patient_id,'test_fold':int(fold[k])})
        cache=(np.concatenate(rep_i),np.concatenate(rep_j),np.concatenate(rep_c),np.concatenate(rep_e));caches.append(cache)
        C=cache[2].mean();E=cache[3].mean()
        allreps.append({'cohort':label,'repeat':repeat,'seed':seed,'comparable_pairs':len(cache[2]),'C_clinical':float(C),'C_extended':float(E),'delta_C':float(E-C)})
        print(label,'seed',seed,'C',round(C,6),round(E,6),'delta',round(E-C,6),flush=True)
    # Bootstrap fixed cross-validated paired predictions, within-fold distinct-patient pairs.
    rng=np.random.default_rng(20260926);weights=rng.multinomial(len(data),np.repeat(1/len(data),len(data)),size=NBOOT)
    bc=[];be=[]
    for i,j,c,e in caches:
        w=weights[:,i].astype(float)*weights[:,j];den=w.sum(axis=1)
        assert np.all(den>0)
        bc.append(np.sum(w*c[None,:],axis=1)/den)
        be.append(np.sum(w*e[None,:],axis=1)/den)
        assert np.all(np.isfinite(bc[-1])) and np.all(np.isfinite(be[-1]))
    Bc=np.mean(bc,axis=0);Be=np.mean(be,axis=0);Bd=Be-Bc
    for b in range(NBOOT):allboot.append({'cohort':label,'bootstrap':b,'C_clinical':float(Bc[b]),'C_extended':float(Be[b]),'delta_C':float(Bd[b])})
    reps=pd.DataFrame([r for r in allreps if r['cohort']==label]);d=reps.delta_C
    allsum.append({'cohort':label,'n':len(data),'deaths':int(data.event.sum()),'LR_chi2_df1':float(lr),'LR_p':float(chi2.sf(lr,1)),'C_clinical_mean':float(reps.C_clinical.mean()),'C_extended_mean':float(reps.C_extended.mean()),'delta_C_mean':float(d.mean()),'delta_C_split_sd':float(d.std(ddof=1)),'delta_C_split_min':float(d.min()),'delta_C_split_max':float(d.max()),'conditional_bootstrap_delta_low':float(np.quantile(Bd,.025)),'conditional_bootstrap_delta_high':float(np.quantile(Bd,.975)),'bootstrap_repeats':NBOOT,'bootstrap_refits_models':False,'concordance_pairs':'within-sex and within-test-fold' if stratified else 'within-test-fold','stratified_by_sex':stratified})

summary=pd.DataFrame(allsum);summary.to_csv(OUT/'incremental_summary.tsv',sep='\t',index=False)
for name,records in [('full_model_statistics',allfull),('implementation_validation',allchecks),('cv_fold_metrics',allfolds),('cv_patient_predictions',allpred),('cv_repetition_metrics',allreps),('cv_coefficients',allcvcoef),('cv_scaling',allscaling),('cv_membership',allmembers),('conditional_bootstrap_fixed_predictions',allboot)]:
    pd.DataFrame(records).to_csv(OUT/f'{name}.tsv',sep='\t',index=False)
pd.concat(allfullcoef).to_csv(OUT/'full_model_coefficients.tsv',sep='\t',index=False)
versions={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'lifelines':lifelines.__version__,'statsmodels':statsmodels.__version__,'matplotlib':matplotlib.__version__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'plan_sha256':manifest['plan_sha256']}
(ROOT/'runtime_versions.json').write_text(json.dumps(versions,indent=2))
(OUT/'incremental_summary.json').write_text(json.dumps(allsum,indent=2))

# Two panels: paired CV discrimination; paired delta with conditional intervals.
plt.rcParams.update({'font.family':'sans-serif','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
s=summary.iloc[:2];fig,(ax,bx)=plt.subplots(1,2,figsize=(8.6,3.3),gridspec_kw={'width_ratios':[1,1]})
for k,row in s.iterrows():
    yy=1-k
    ax.plot([row.C_clinical_mean,row.C_extended_mean],[yy,yy],color='#9FB1BA',lw=2)
    ax.scatter(row.C_clinical_mean,yy,marker='s',s=32,color='#657A87',label='Clinical' if k==0 else None)
    ax.scatter(row.C_extended_mean,yy,marker='o',s=32,color='#136B80',label='Clinical + SPINK5' if k==0 else None)
    bx.errorbar(row.delta_C_mean,yy,xerr=[[row.delta_C_mean-row.conditional_bootstrap_delta_low],[row.conditional_bootstrap_delta_high-row.delta_C_mean]],fmt='o',color='#136B80',capsize=3)
    bx.text(.98,yy+.20,f"ΔC {row.delta_C_mean:+.3f}; LR P={row.LR_p:.3g}",transform=bx.get_yaxis_transform(),ha='right',fontsize=8)
ax.set_yticks([1,0],['TCGA-HNSC\n491 patients / 210 deaths','GSE65858\n233 patients / 73 deaths'])
ax.set_xlim(.5,.8);ax.set_ylim(-.6,1.55);ax.set_xlabel('Internally validated Harrell C')
ax.legend(loc='upper center',bbox_to_anchor=(.52,1.17),frameon=False,ncol=2,fontsize=8)
bx.axvline(0,color='#888888',lw=.8,ls='--');bx.set_yticks([1,0],[]);bx.set_ylim(-.6,1.55);bx.set_xlabel('Paired increment in Harrell C')
for a in [ax,bx]:a.spines[['top','right']].set_visible(False)
fig.text(.5,.02,'Means from 10 × 5-fold paired CV. Intervals resample patients with fixed CV predictions; models are not refitted.',ha='center',fontsize=7)
fig.subplots_adjust(left=.21,right=.98,top=.80,bottom=.24,wspace=.38)
fig.savefig(OUT/'incremental_discrimination.pdf');fig.savefig(OUT/'incremental_discrimination.png',dpi=240);plt.close(fig)
print(summary.to_string(index=False),flush=True)
print('ALL DONE',flush=True)
