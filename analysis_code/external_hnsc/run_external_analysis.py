"""Locked single-gene external association check; no outcome-based selection."""
from pathlib import Path
import csv,gzip,io,json,hashlib,sys,warnings
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from statsmodels.duration.hazard_regression import PHReg
BASE=Path(__file__).resolve().parent

def parse():
    meta={}; rows=[];processing='';in_table=False
    with gzip.open(BASE/'GSE65858_series_matrix.txt.gz','rt') as f:
        for l in f:
            if l.startswith('!series_matrix_table_begin'):in_table=True;continue
            if l.startswith('!series_matrix_table_end'):break
            if in_table:rows.append(l);continue
            row=next(csv.reader([l],delimiter='\t'))
            if len(row)<2:continue
            if row[0]=='!Sample_geo_accession':meta['geo']=row[1:]
            if row[0]=='!Sample_characteristics_ch1':
                k=row[1].split(': ',1)[0];meta[k]=[v.split(': ',1)[1] for v in row[1:]]
            if row[0]=='!Sample_data_processing':processing=row[1]
    ph=pd.DataFrame(meta).set_index('geo').replace('NA',np.nan)
    assert ph.index.is_unique and ph.ID.is_unique
    ph.to_csv(BASE/'sample_metadata.csv')
    expr=pd.read_csv(io.StringIO(''.join(rows)),sep='\t',index_col=0)
    assert expr.index.is_unique and expr.columns.is_unique
    assert expr.shape==(31330,270)
    assert set(expr.columns)==set(ph.index)
    (BASE/'processing.txt').write_text(processing)
    # The exact GEO-deposited Illumina manufacturer annotation. Header detected explicitly.
    annfile=BASE/'GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz'
    with gzip.open(annfile,'rt',errors='replace') as f:
        lines=f.readlines()
    header=next(i for i,l in enumerate(lines) if 'Probe_Id' in l and 'Symbol' in l)
    end=next((i for i in range(header+1,len(lines)) if lines[i].startswith('[')),len(lines))
    ann=pd.read_csv(io.StringIO(''.join(lines[header:end])),sep='\t',dtype=str)
    probes=ann.loc[ann['Symbol'].fillna('').str.strip().eq('SPINK5')].copy()
    assert len(probes)>0
    probes['present_in_series_matrix']=probes['Probe_Id'].isin(expr.index)
    probes.to_csv(BASE/'SPINK5_probe_mapping.csv',index=False)
    probe_ids=probes.loc[probes.present_in_series_matrix,'Probe_Id'].unique().tolist()
    assert probe_ids
    ph['SPINK5']=expr.loc[probe_ids].median(axis=0).reindex(ph.index)
    expr.loc[probe_ids].to_csv(BASE/'SPINK5_probe_expression.csv')
    ph['os_native']=pd.to_numeric(ph['os'],errors='coerce')
    ph['death']=ph['os_event'].map({'TRUE':1,'FALSE':0})
    ph['age10']=pd.to_numeric(ph.age,errors='coerce')/10
    ph['male']=ph.gender.map({'M':1,'F':0})
    ph['stage']=ph.uicc_stage.map({'I':'I','II':'II','III':'III','IVA':'IV','IVB':'IV','IVC':'IV'})
    ph['hpv_active']=ph.hpv16_dna_rna.map({'DNA+RNA+':1,'DNA-':0,'DNA+RNA-':0})
    return ph,probe_ids,probes

ph,probe_ids,mapping=parse()
info={'n_total':len(ph),'n_events_total':int(ph.death.sum()),'probes':probe_ids,'n_probes_annotation':len(mapping),
      'time_definition':'GEO os and os_event; native recorded follow-up units; unit not explicit in matrix',
      'processing':'GEO deposited log2 + RSN + ComBat processed values; no further log transform',
      'checks':{'unique_geo':bool(ph.index.is_unique),'unique_patient_ID':bool(ph.ID.is_unique)}}
for k in ['tumor_type','gender','stage','hpv16_dna_rna','hpv_dna']:
    info[k+'_counts']={str(a):int(b) for a,b in ph[k].value_counts(dropna=False).items()}
models=[('primary_unadjusted',True,[]),('primary_clinical',True,['age10','male','stage']),
        ('primary_clinical_hpv',True,['age10','male','stage','hpv_active']),
        ('all_tumors_clinical_hpv',False,['age10','male','stage','hpv_active'])]
results=[];validation=[];phs=[];missing=[];all_observations=[]
for name,primary,covs in models:
    cohort=ph.loc[ph.tumor_type.eq('Primary')].copy() if primary else ph.copy()
    mean=cohort.SPINK5.mean();sd=cohort.SPINK5.std(ddof=1)
    cohort['SPINK5_z']=(cohort.SPINK5-mean)/sd
    needed=['os_native','death','SPINK5_z']+covs
    valid=cohort[needed].notna().all(axis=1)&cohort.os_native.gt(0)&cohort.death.isin([0,1])
    for col in needed:
        missing.append({'model':name,'variable':col,'missing':int(cohort[col].isna().sum()),'cohort_n':len(cohort)})
    d=cohort.loc[valid,needed].copy()
    if 'stage' in d:
        stage=pd.Categorical(d.pop('stage'),categories=['I','II','III','IV'],ordered=True)
        stage=pd.get_dummies(stage,prefix='stage',drop_first=True,dtype=float);stage.index=d.index
        d=pd.concat([d,stage],axis=1)
    d=d.astype(float)
    fit=CoxPHFitter().fit(d,duration_col='os_native',event_col='death')
    fit.summary.to_csv(BASE/f'{name}_all_coefficients.csv')
    row=fit.summary.loc['SPINK5_z']
    results.append({'model':name,'n':len(d),'deaths':int(d.death.sum()),'excluded':int(len(cohort)-len(d)),
                    'hr':row['exp(coef)'],'ci_low':row['exp(coef) lower 95%'],'ci_high':row['exp(coef) upper 95%'],
                    'p':row['p'],'beta':row['coef'],'se':row['se(coef)'],'c_index_apparent':fit.concordance_index_,
                    'SPINK5_mean':mean,'SPINK5_sd':sd,'covariates':';'.join(covs)})
    ph_test=proportional_hazard_test(fit,d,time_transform='rank').summary.reset_index().rename(columns={'index':'covariate'})
    ph_test['model']=name;phs.append(ph_test)
    x=d.drop(columns=['os_native','death'])
    check=PHReg(d.os_native,x,status=d.death,ties='efron').fit()
    j=x.columns.get_loc('SPINK5_z')
    validation.append({'model':name,'beta_lifelines':row['coef'],'beta_statsmodels':check.params[j],
                       'beta_abs_difference':abs(row['coef']-check.params[j]),'se_statsmodels':check.bse[j],
                       'p_statsmodels':check.pvalues[j]})
    assert abs(row['coef']-check.params[j])<1e-4
    incl=cohort.copy();incl['included']=valid;incl['model']=name;all_observations.append(incl)
    d.to_csv(BASE/f'{name}_model_input.csv')
res=pd.DataFrame(results);res.to_csv(BASE/'external_survival_results.csv',index=False)
pd.DataFrame(validation).to_csv(BASE/'independent_PHReg_verification.csv',index=False)
pd.concat(phs).to_csv(BASE/'proportional_hazards_tests.csv',index=False)
pd.DataFrame(missing).to_csv(BASE/'missingness_by_model.csv',index=False)
pd.concat(all_observations).to_csv(BASE/'analysis_observations_and_exclusions.csv')
info['expression_summary_primary']=ph.loc[ph.tumor_type.eq('Primary'),'SPINK5'].describe().to_dict()
info['models']=results
(BASE/'analysis_summary.json').write_text(json.dumps(info,indent=2))
files=['GSE65858_series_matrix.txt.gz','GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz','analysis_plan.md','run_external_analysis.py']
(BASE/'SHA256SUMS.txt').write_text('\n'.join(hashlib.sha256((BASE/f).read_bytes()).hexdigest()+'  '+f for f in files)+'\n')
print(res.to_string(index=False))
print('probe IDs:',probe_ids)
