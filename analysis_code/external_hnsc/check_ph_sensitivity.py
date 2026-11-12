"""Diagnostic-triggered sensitivity: primary model rank PH test flagged sex (P=.0271).
Retain primary model; stratify baseline hazard by sex without estimating its coefficient.
"""
from pathlib import Path
import pandas as pd,json
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from statsmodels.duration.hazard_regression import PHReg
B=Path(__file__).resolve().parent
d=pd.read_csv(B/'primary_clinical_hpv_model_input.csv',index_col=0)
c=CoxPHFitter().fit(d,'os_native','death',strata=['male'])
c.summary.to_csv(B/'sex_stratified_PH_sensitivity_coefficients.csv')
p=proportional_hazard_test(c,d,time_transform='rank').summary
p.to_csv(B/'sex_stratified_PH_sensitivity_tests.csv')
r=c.summary.loc['SPINK5_z'];x=d.drop(columns=['os_native','death','male'])
s=PHReg(d.os_native,x,status=d.death,strata=d.male,ties='efron').fit()
j=x.columns.get_loc('SPINK5_z')
v={'n':len(d),'deaths':int(d.death.sum()),'hr':r['exp(coef)'],'ci_low':r['exp(coef) lower 95%'],'ci_high':r['exp(coef) upper 95%'],'p':r['p'],'beta_statsmodels':s.params[j],'beta_lifelines':r['coef'],'beta_abs_difference':abs(s.params[j]-r['coef']),'SPINK5_rank_PH_p':p.loc['SPINK5_z','p']}
assert v['beta_abs_difference']<1e-4
(B/'sex_stratified_PH_sensitivity.json').write_text(json.dumps(v,indent=2))
print(json.dumps(v,indent=2))
