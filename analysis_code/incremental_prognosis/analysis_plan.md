# Locked incremental prognosis analysis

Written 2026-09-26 before fitting any model in this analysis.

## Question and scope

For identical patients within each of two existing cohorts, does adding one continuous SPINK5 term to prespecified clinical Cox covariates add conditional survival information and internally validated risk discrimination?

This is a new secondary analysis of previously inspected public cohorts. It is not prospective preregistration. The following specification is fixed before examining the incremental results. It refits models independently within each cohort. It is not external validation of a frozen TCGA prediction model, drug-response prediction, or treatment-effect estimation.

## Input populations and fixed comparisons

- TCGA-HNSC: the exact saved `analysis/public_hnsc/results/model_data_M2_age_clinical_stage_sex_smoking_site.tsv`, 491 primary-tumour patients, 210 deaths. Clinical baseline: age per 10 years, sex, ever smoking, categorical clinical stage I–IV, anatomical site (oral cavity reference; larynx, oropharynx, hypopharynx indicators). Extended: the same baseline plus one continuous SPINK5 term.
- GSE65858: the exact saved `analysis/external_hnsc/primary_clinical_hpv_model_input.csv`, 233 primary-tumour patients, 73 deaths. Baseline: age per 10 years, sex, categorical UICC stage I–IV, transcriptionally active HPV16. Extended: the same baseline plus one continuous SPINK5 term.
- No outcome-based probe choice, expression cutoff, feature selection, additional patient exclusion, or selection of the better cohort/model/split. No time-unit conversion, pooled-cohort survival fit, or fixed-year outcome metric.
- Primary output retains both cohorts regardless of effect size or significance. Existing missingness/exclusion rules remain unchanged.

## Full-cohort conditional-information comparison

Fit unpenalized Cox proportional-hazards models with Efron ties on exactly the same patients. Report `2*(loglik_extended-loglik_clinical)` with one degree of freedom, P value, and both complete coefficient tables. Independently verify likelihoods and coefficients with lifelines and statsmodels PHReg. This tests one additional conditional association term and is not, by itself, a prediction-improvement test.

## Paired internal validation

Use 10 repetitions of five-fold cross-validation, event-stratified, with locked seeds 42–51. Each patient's record appears in exactly one held-out fold per repetition. Baseline and extended models share every fold. Standardize continuous age and SPINK5 on training-fold means and population SDs only; apply the frozen fold parameters to held-out patients. Binary and categorical indicators are unchanged. Models use no penalties or tuning. A training-constant feature is dropped for that fold, recorded, and identically omitted from prediction; this rule is defined before fitting.

Compute Harrell's C on held-out patients within each fold. Risk scores from independently fitted folds are not pooled for a cross-fold concordance calculation. Combine fold estimates by the number of comparable survival pairs. Report each cohort's mean C_clinical, mean C_extended, paired mean delta_C, SD and full range across the 10 repetitions, together with each fold and repetition result. Repeated partitions are not independent patients: the SD/range is split sensitivity, not a confidence interval. Count censoring/time ties consistently with the saved implementation and verify concordance against lifelines.

Retain per-patient per-fold baseline and extended predictions, train/test IDs, fold scaling, coefficients, fit warnings and every result. Fit failure is reported; do not switch algorithms or select successful seeds to improve estimates.

## Paired prediction uncertainty (secondary)

As a compact secondary assessment, bootstrap patients 1000 times with seed 20260926, using the **already saved cross-validated predictions** and preserving a resampled patient's multiplicity and paired predictions for both models across all repetitions. Compute weighted concordance from comparable distinct-patient pairs within each held-out fold; never compare predictions across folds or count copies of the same individual as comparable independent outcomes. Recompute the average paired delta_C across repetitions for each resample. Report the percentile interval as a **conditional patient-resampling interval for fixed cross-validated predictions**, not a full model-redevelopment interval. Models are not refitted in this bootstrap. This interval excludes uncertainty from retraining and reselecting the evaluation design and cannot establish clinical utility.

## Sex proportional-hazards sensitivity in GSE65858

An existing diagnostic flagged sex in the GSE65858 primary model. Preserve the prespecified main analysis above. Separately compare the same two models with a sex-stratified baseline hazard, retaining age, stage and active HPV16. For this diagnostic sensitivity, use the same patients and CV splits, compare only within-sex pairs (partial hazards from separate baseline strata do not rank absolute cross-sex risk), and label that discrimination target separately. Verify the full-cohort stratified likelihood ratio using statsmodels PHReg strata. This analysis is retained regardless of its result and does not replace the main estimate.

## Presentation

One table and a compact plot will show the primary two-cohort C values/differences and full-data likelihood-ratio results. Any conditional bootstrap interval will be explicitly qualified. Neither an observed association nor a C-index change is described as clinical net benefit, therapy selection, or AI drug discovery. Dependencies, source paths, source hashes, input copies, exact code and outputs are saved under this directory.
