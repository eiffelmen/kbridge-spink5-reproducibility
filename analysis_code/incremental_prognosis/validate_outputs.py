"""Offline independent checks of exported membership, scaling, predictions and metrics."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from lifelines.utils import concordance_index
R=Path(__file__).resolve().parent; O=R/'results'
p=pd.read_csv(O/'cv_patient_predictions.tsv',sep='\t')
f=pd.read_csv(O/'cv_fold_metrics.tsv',sep='\t')
s=pd.read_csv(O/'cv_scaling.tsv',sep='\t')
c=pd.read_csv(O/'cv_coefficients.tsv',sep='\t')
m=pd.read_csv(O/'cv_membership.tsv',sep='\t')
rep=pd.read_csv(O/'cv_repetition_metrics.tsv',sep='\t')
summary=pd.read_csv(O/'incremental_summary.tsv',sep='\t')
assert len(p)==9570 and len(f)==150 and len(rep)==30
assert p.groupby(['cohort','repeat','patient_id']).size().eq(1).all()
assert m.groupby(['cohort','repeat','patient_id']).size().eq(1).all()
assert not f.filter(like='warnings_').notna().any().any()
assert not f.dropped_training_constant_features.notna().any()
max_scaling=0.;max_pred_const_resid=0.;max_metric_diff=0.
for label in ['TCGA','GSE65858','GSE65858_sex_stratified']:
    d=pd.read_csv(R/'inputs'/('TCGA_verified_analysis_input.tsv' if label=='TCGA' else 'GSE65858_verified_analysis_input.tsv'),sep='\t').set_index('patient_id')
    assert len(d)==(491 if label=='TCGA' else 233)
    for repeat_id in range(10):
        members=m[(m.cohort==label)&(m.repeat==repeat_id)].set_index('patient_id')
        folds=f[(f.cohort==label)&(f.repeat==repeat_id)]
        for fold_id in range(5):
            train=d.loc[members.index[members.test_fold!=fold_id]]
            test=d.loc[members.index[members.test_fold==fold_id]].copy()
            pred=p[(p.cohort==label)&(p.repeat==repeat_id)&(p.fold==fold_id)].set_index('patient_id').loc[test.index]
            for term in ['age10','SPINK5']:
                stat=s[(s.cohort==label)&(s.repeat==repeat_id)&(s.fold==fold_id)&(s.feature==term)].iloc[0]
                max_scaling=max(max_scaling,abs(train[term].mean()-stat.training_mean),abs(train[term].std(ddof=0)-stat.training_population_sd))
                test[term]=(test[term]-stat.training_mean)/stat.training_population_sd
            for typ,col in [('clinical','risk_clinical'),('clinical_SPINK5','risk_extended')]:
                coef=c[(c.cohort==label)&(c.repeat==repeat_id)&(c.fold==fold_id)&(c.model==typ)].set_index('feature').coefficient
                # lifelines subtracts a constant training-centering offset. Verify prediction differences.
                reconstructed=np.sum(test[coef.index].to_numpy()*coef.to_numpy()[None,:],axis=1)
                resid=reconstructed-pred[col].to_numpy()
                max_pred_const_resid=max(max_pred_const_resid,float(np.ptp(resid)))
                if label!='GSE65858_sex_stratified':
                    metric=concordance_index(test.time,-pred[col],test.event)
                    saved=folds.loc[folds.fold==fold_id,'C_clinical' if typ=='clinical' else 'C_extended'].iloc[0]
                    max_metric_diff=max(max_metric_diff,abs(metric-saved))
        for typ in ['clinical','extended']:
            weighted=(folds[f'C_{typ}']*folds.comparable_pairs).sum()/folds.comparable_pairs.sum()
            saved=rep.loc[(rep.cohort==label)&(rep.repeat==repeat_id),f'C_{typ}'].iloc[0]
            assert abs(weighted-saved)<1e-12
    sr=summary.loc[summary.cohort==label].iloc[0]
    rr=rep[rep.cohort==label]
    assert abs(sr.delta_C_mean-rr.delta_C.mean())<1e-12
    assert abs(sr.C_clinical_mean-rr.C_clinical.mean())<1e-12
    assert abs(sr.C_extended_mean-rr.C_extended.mean())<1e-12
assert max_scaling<1e-12 and max_pred_const_resid<1e-9 and max_metric_diff<1e-12
manifest=json.loads((R/'source_manifest.json').read_text())
for x in manifest['sources']:
    assert hashlib.sha256((R/x['copy']).read_bytes()).hexdigest()==x['sha256']
checks={'status':'PASS','patient_prediction_rows':len(p),'fold_rows':len(f),'repetition_rows':len(rep),'unique_once_per_repetition':True,'source_copy_hashes_match':True,'fit_warnings':0,'training_constant_drops':0,'max_train_scaling_difference':max_scaling,'max_reconstructed_prediction_difference_after_constant_centering':max_pred_const_resid,'max_concordance_difference_vs_lifelines':max_metric_diff}
(O/'validation_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
