#!/usr/bin/env python3
"""Reanalyse public TCGA-HNSC clinical and SPINK5 expression data."""
from pathlib import Path
import json, hashlib, sys, platform, warnings
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from statsmodels.duration.hazard_regression import PHReg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lifelines, scipy, statsmodels

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
clinical=pd.read_csv(DATA/'HNSC_clinicalMatrix.tsv',sep='\t',na_values=['[Discrepancy]','[Not Available]','[Not Applicable]','[Unknown]'])
sample_ids=json.loads((DATA/'expression_sample_ids.json').read_text())
values=json.loads((DATA/'SPINK5_all_samples.json').read_text())[0]
expr=pd.DataFrame({'sampleID':sample_ids,'SPINK5_log2':pd.to_numeric(values,errors='coerce')})
expr['sample_code']=expr.sampleID.str[-2:]
expr['patient']=expr.sampleID.str[:12]
expr.to_csv(OUT/'expression_spink5_all_samples.tsv',sep='\t',index=False)
assert expr.sampleID.is_unique
assert clinical.sampleID.is_unique
x=expr[expr.sample_code.eq('01')].merge(clinical,on='sampleID',how='left',validate='one_to_one')
assert x.patient.is_unique
x['event']=x.vital_status.map({'DECEASED':1,'LIVING':0})
x['time']=np.where(x.event.eq(1),pd.to_numeric(x.days_to_death,errors='coerce'),pd.to_numeric(x.days_to_last_followup,errors='coerce'))
x['age10']=pd.to_numeric(x.age_at_initial_pathologic_diagnosis,errors='coerce')/10
stage_map={'Stage I':1,'Stage II':2,'Stage III':3,'Stage IVA':4,'Stage IVB':4,'Stage IVC':4}
x['clinical_stage_num']=x.clinical_stage.map(stage_map)
x['pathologic_stage_num']=x.pathologic_stage.map(stage_map)
x['stage_combined_num']=x.pathologic_stage_num.fillna(x.clinical_stage_num)
x['male']=x.gender.map({'MALE':1,'FEMALE':0})
x['ever_smoked']=x.tobacco_smoking_history.map({1:0,2:1,3:1,4:1,5:1})
x['HPV_p16_positive']=x.hpv_status_by_p16_testing.map({'Negative':0,'Positive':1})
x['HPV_ISH_positive']=x.hpv_status_by_ish_testing.map({'Negative':0,'Positive':1})
site_map={s:'oral_cavity' for s in ['Oral Tongue','Oral Cavity','Floor of mouth','Buccal Mucosa','Alveolar Ridge','Hard Palate','Lip']}
site_map.update({s:'oropharynx' for s in ['Tonsil','Base of tongue','Oropharynx']})
site_map.update({'Larynx':'larynx','Hypopharynx':'hypopharynx'})
x['site_group']=x.anatomic_neoplasm_subdivision.map(site_map)
# A single population SD is kept across all restricted models.
spink_mean=x.SPINK5_log2.mean(); spink_sd=x.SPINK5_log2.std(ddof=0)
x['SPINK5_per_SD']=(x.SPINK5_log2-spink_mean)/spink_sd
base=x[x.event.notna() & x.time.notna() & (x.time>0)].copy()
flow={'expression_samples':len(expr),'primary_tumour_code01':int((expr.sample_code=='01').sum()),'normal_code11':int((expr.sample_code=='11').sum()),'metastatic_code06':int((expr.sample_code=='06').sum()),'primary_with_expression_clinical':len(x),'missing_or_nonpositive_os':int(len(x)-len(base)),'os_eligible_primary':len(base),'os_eligible_deaths':int(base.event.sum()),'SPINK5_mean_log2':float(spink_mean),'SPINK5_sd_log2_population':float(spink_sd)}
misscols=['SPINK5_log2','event','time','age10','clinical_stage_num','pathologic_stage_num','stage_combined_num','male','ever_smoked','HPV_p16_positive','HPV_ISH_positive','site_group']
missing=pd.DataFrame([{'variable':c,'available':int(x[c].notna().sum()),'missing':int(x[c].isna().sum()),'missing_pct':float(x[c].isna().mean()*100)} for c in misscols])
missing.to_csv(OUT/'missingness_primary_tumours.tsv',sep='\t',index=False)
x.to_csv(OUT/'analysis_primary_tumours.tsv',sep='\t',index=False)
records=[]; allcoefs=[]; phrecords=[]; validation=[]

def design(df, stage=None, extras=(), ordinal=False):
    continuous=['SPINK5_per_SD'] + (['age10'] if stage else []) + [c for c in extras if c!='site_group']
    required=['time','event']+continuous+([stage] if stage else [])+(['site_group'] if 'site_group' in extras else [])
    d=df.dropna(subset=required).copy()
    z=d[['time','event']+continuous].astype(float)
    if stage:
        if ordinal: z['stage_ordinal']=d[stage].astype(float)
        else:
            for s in (2,3,4): z[f'stage_{s}_vs_1']=(d[stage]==s).astype(float)
    if 'site_group' in extras:
        for site in ['larynx','oropharynx','hypopharynx']:
            z[f'site_{site}_vs_oral_cavity']=(d.site_group==site).astype(float)
    return d,z

def fit(name,df,stage=None,extras=(),ordinal=False):
    d,z=design(df,stage,extras,ordinal)
    model=CoxPHFitter(penalizer=0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model.fit(z,duration_col='time',event_col='event',show_progress=False)
    s=model.summary.loc['SPINK5_per_SD']
    record={'model':name,'n':len(z),'deaths':int(z.event.sum()),'coefficients':len(model.params_),'HR_per_SD':float(s['exp(coef)']),'CI95_low':float(s['exp(coef) lower 95%']),'CI95_high':float(s['exp(coef) upper 95%']),'p':float(s['p']),'log_likelihood':float(model.log_likelihood_),'covariates':';'.join(model.params_.index),'warnings':'; '.join(str(w.message) for w in caught)}
    records.append(record)
    summary=model.summary.reset_index(); summary.insert(0,'model',name); allcoefs.append(summary)
    ph=proportional_hazard_test(model,z,time_transform='rank').summary.reset_index(); ph.insert(0,'model',name); phrecords.append(ph)
    z.to_csv(OUT/f'model_data_{name}.tsv',sep='\t',index=False)
    # Independent numerical implementation in statsmodels, same Efron ties.
    sm=PHReg(z.time,z.drop(columns=['time','event']),status=z.event,ties='efron').fit(disp=False)
    validation.append({'model':name,'max_abs_coef_difference_lifelines_statsmodels':float(np.max(np.abs(sm.params-model.params_.to_numpy()))),'max_abs_se_difference':float(np.max(np.abs(sm.bse-model.standard_errors_.to_numpy())))})
    return d,z,model

fit('M0_unadjusted',base)
fit('M1_age_clinical_stage',base,'clinical_stage_num')
full_d,full_z,full_model=fit('M2_age_clinical_stage_sex_smoking_site',base,'clinical_stage_num',('male','ever_smoked','site_group'))
fit('M1_same_complete_cases_as_M2',full_d,'clinical_stage_num')
fit('S1_age_pathologic_stage',base,'pathologic_stage_num')
fit('S2_age_pathologic_then_clinical_stage',base,'stage_combined_num')
# HPV sensitivity has sparse outcomes; use a parsimonious ordinal-stage model and retain its matched comparator.
hpv=base[base.HPV_p16_positive.notna()].copy()
fit('S3_p16_subset_age_stage_ordinal',hpv,'clinical_stage_num',ordinal=True)
fit('S4_p16_subset_age_stage_ordinal_HPV',hpv,'clinical_stage_num',('HPV_p16_positive',),ordinal=True)
# Exclude oropharyngeal sites to assess sensitivity to HPV-associated case composition.
fit('S5_nonoropharynx_age_clinical_stage_sex_smoking',base[base.site_group.notna() & base.site_group.ne('oropharynx')],'clinical_stage_num',('male','ever_smoked'),ordinal=False)

res=pd.DataFrame(records); res.to_csv(OUT/'cox_model_summary.tsv',sep='\t',index=False)
pd.concat(allcoefs).to_csv(OUT/'cox_all_coefficients.tsv',sep='\t',index=False)
pd.concat(phrecords).to_csv(OUT/'proportional_hazards_rank_tests.tsv',sep='\t',index=False)
pd.DataFrame(validation).to_csv(OUT/'independent_numerical_validation.tsv',sep='\t',index=False)

# Paired tumour-normal comparison among patients with both specimens.
tum=expr[expr.sample_code=='01']; normal=expr[expr.sample_code=='11']
pairs=tum[['patient','SPINK5_log2']].merge(normal[['patient','SPINK5_log2']],on='patient',suffixes=('_tumour','_normal'),validate='one_to_one')
pairs['tumour_minus_normal']=pairs.SPINK5_log2_tumour-pairs.SPINK5_log2_normal
pairs.to_csv(OUT/'paired_tumour_normal.tsv',sep='\t',index=False)
wx=stats.wilcoxon(pairs.SPINK5_log2_tumour,pairs.SPINK5_log2_normal,alternative='two-sided')
comparison={'paired_n':len(pairs),'paired_median_log2_difference':float(pairs.tumour_minus_normal.median()),'wilcoxon_statistic':float(wx.statistic),'wilcoxon_two_sided_p':float(wx.pvalue),'unpaired_tumour_n':len(tum),'unpaired_normal_n':len(normal),'tumour_median':float(tum.SPINK5_log2.median()),'normal_median':float(normal.SPINK5_log2.median())}
flow['comparison']=comparison
flow['versions']={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'lifelines':lifelines.__version__,'statsmodels':statsmodels.__version__}
(OUT/'analysis_summary.json').write_text(json.dumps(flow,indent=2))

# Modest, vector, reproducible forest plot. Every planned analysis is retained.
plt.rcParams.update({'font.family':'sans-serif','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
labels=['Unadjusted','Age + clinical stage','+ sex, smoking, anatomical site','Age + stage, same complete cases','Age + pathological stage','Age + stage with clinical fallback','p16 subset: age + ordinal stage','p16 subset: + HPV p16','Non-oropharynx: age, stage, sex, smoking']
fig,ax=plt.subplots(figsize=(9.0,5.2)); y=np.arange(len(res))[::-1]
ax.errorbar(res.HR_per_SD,y,xerr=np.vstack([res.HR_per_SD-res.CI95_low,res.CI95_high-res.HR_per_SD]),fmt='o',color='#1A4D6E',ecolor='#456A80',capsize=3,markersize=4)
ax.axvline(1,color='#999999',ls='--',lw=.8)
ax.set_yticks(y,labels); ax.set_xlabel('Overall-survival hazard ratio per 1 SD higher SPINK5 expression (95% CI)')
ax.set_xscale('log'); ax.set_xlim(.35,1.65); ax.set_xticks([.4,.6,.8,1,1.5],[.4,.6,.8,1,1.5])
for i,r in res.iterrows(): ax.text(1.04,y[i],f"{r.HR_per_SD:.2f} ({r.CI95_low:.2f}–{r.CI95_high:.2f})\n{r.n} patients / {r.deaths} deaths",transform=ax.get_yaxis_transform(),va='center',fontsize=8)
ax.spines[['top','right']].set_visible(False); fig.subplots_adjust(left=.38,right=.75,bottom=.16,top=.94)
fig.savefig(OUT/'spink5_cox_sensitivity.pdf'); fig.savefig(OUT/'spink5_cox_sensitivity.png',dpi=240); plt.close(fig)
print(res[['model','n','deaths','HR_per_SD','CI95_low','CI95_high','p']].to_string(index=False))
print(json.dumps(flow,indent=2))
print('Independent fit checks:',validation)
