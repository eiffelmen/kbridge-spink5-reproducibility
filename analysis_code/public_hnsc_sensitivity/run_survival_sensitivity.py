"""Exploratory sensitivity analyses of the two published HNSCC cohorts."""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from patsy import dmatrix
from scipy.stats import chi2, spearmanr


def fit(data, strata=None):
    return CoxPHFitter().fit(data, 'time', 'event', strata=strata)


def summarise(results, cohort, model, data, fit_object):
    row = fit_object.summary.loc['SPINK5']
    results.append({
        'cohort': cohort, 'model': model, 'n': len(data), 'deaths': int(data.event.sum()),
        'HR_per_SD': row['exp(coef)'], 'CI_low': row['exp(coef) lower 95%'],
        'CI_high': row['exp(coef) upper 95%'], 'P': row['p'], 'Wald_z': row['z'],
    })


def analyse(args):
    source = args.inputs
    results, nonlinearity, diagnostics, inclusion, correlations = [], [], [], [], []
    tc_meta = pd.read_csv(source / 'tcga_clinical.tsv', sep='\t', index_col='sampleID')
    tc = pd.read_csv(source / 'tcga_main.tsv', sep='\t', index_col='patient_id')
    gs_meta = pd.read_csv(source / 'gse65858_metadata.csv', index_col='geo')
    gs = pd.read_csv(source / 'gse65858_main.csv', index_col='geo').rename(
        columns={'os_native': 'time', 'death': 'event', 'SPINK5_z': 'SPINK5'}
    )

    for name, data in [('TCGA-HNSC', tc), ('GSE65858', gs)]:
        baseline = fit(data)
        summarise(results, name, 'Primary model', data, baseline)
        ph = proportional_hazard_test(baseline, data, time_transform='rank').summary
        for var, row in ph.iterrows():
            diagnostics.append({
                'cohort': name, 'variable': var, 'chi_square': row['test_statistic'],
                'P': row['p'], 'degrees_of_freedom': 1,
            })

        knots = data.SPINK5.quantile([.05, .35, .65, .95]).to_numpy()
        basis = dmatrix('cr(x, knots=interior, lower_bound=lo, upper_bound=hi, constraints="center") - 1', {
            'x': data.SPINK5, 'interior': knots[1:3], 'lo': knots[0], 'hi': knots[3],
        }, return_type='dataframe')
        basis.columns = [f'spline_{i}' for i in range(basis.shape[1])]
        spline = fit(data.drop(columns='SPINK5').join(basis))
        statistic = 2 * (spline.log_likelihood_ - baseline.log_likelihood_)
        degrees = basis.shape[1] - 1
        nonlinearity.append({
            'cohort': name, 'n': len(data), 'deaths': int(data.event.sum()),
            'LR_chi_square': statistic, 'degrees_of_freedom': degrees,
            'P': chi2.sf(statistic, degrees), **{f'knot_{j+1}': v for j, v in enumerate(knots)},
        })

        metadata = tc_meta if name == 'TCGA-HNSC' else gs_meta
        eligible = metadata.loc[metadata.time.gt(0) & metadata.event.isin([0, 1])] if name == 'TCGA-HNSC' else metadata.loc[metadata.tumor_type.eq('Primary')].assign(
            event=lambda frame: frame.os_event.astype(str).str.upper().map({'TRUE': 1, 'FALSE': 0}))
        subgroups = {'Included': eligible.index.isin(data.index), 'Excluded': ~eligible.index.isin(data.index)}
        if name == 'TCGA-HNSC':
            subgroups['p16 known'] = eligible.HPV_p16_positive.notna()
            subgroups['p16 unknown'] = eligible.HPV_p16_positive.isna()
        for group, selected in subgroups.items():
            frame = eligible.loc[selected]
            inclusion.append({'cohort': name, 'group': group, 'n': len(frame), 'deaths': int(frame.event.sum())})

    grade = tc_meta.neoplasm_histologic_grade.map({'G1':'G1', 'G2':'G2', 'G3':'G3-4', 'G4':'G3-4'})
    grade_cases = tc.loc[grade.reindex(tc.index).notna()].copy()
    summarise(results, 'TCGA-HNSC', 'Primary model on known-grade cases', grade_cases, fit(grade_cases))
    indicators = pd.get_dummies(grade.reindex(grade_cases.index), prefix='grade', dtype=float).drop(columns='grade_G1')
    grade_model = grade_cases.join(indicators)
    summarise(results, 'TCGA-HNSC', 'Add histologic grade', grade_model, fit(grade_model))

    recorded_treatment = gs_meta.treatment.reindex(gs.index)
    treatment_cases = gs.loc[recorded_treatment.isin(['mono', 'multi'])].copy()
    summarise(results, 'GSE65858', 'Primary model on mono-/multimodal cases', treatment_cases, fit(treatment_cases))
    treatment_model = treatment_cases.assign(treatment=recorded_treatment.reindex(treatment_cases.index))
    summarise(results, 'GSE65858', 'Stratify by recorded treatment category', treatment_model, fit(treatment_model, ['treatment']))
    sex_model = gs.copy()
    summarise(results, 'GSE65858', 'Stratify by sex', sex_model, fit(sex_model, ['male']))

    gene_definitions = pd.read_csv(source / 'mcp_counter_genes.tsv', sep='\t')
    signatures = {group: set(gene_definitions.loc[gene_definitions['Cell population'].eq(group), 'HUGO symbols'])
                  for group in ['B lineage', 'Fibroblasts']}
    markers = set.union(*signatures.values(), {'SPINK5'})
    if args.xena_expression:
        with gzip.open(args.xena_expression, 'rt') as stream:
            header = next(stream)
            rows = [line for line in stream if line.split('\t', 1)[0] in markers]
        (source / 'tcga_marker_expression.tsv').write_text(header + ''.join(rows), encoding='utf-8')
    expr = pd.read_csv(source / 'tcga_marker_expression.tsv', sep='\t', index_col='sample').T
    matched = expr.reindex(tc_meta.index)
    # The Xena API export retains four significant digits, unlike the full matrix.
    np.testing.assert_allclose(matched.SPINK5, tc_meta.SPINK5_log2, rtol=5e-4, atol=1e-6)
    score_frame = pd.DataFrame(index=matched.index)
    score_frame['SPINK5_log2'] = matched.SPINK5
    for population, genes in signatures.items():
        available = sorted(genes.intersection(matched.columns))
        values = matched[available].mean(axis=1, skipna=False)
        score_frame[population] = values
        usable = values.notna() & matched.SPINK5.notna()
        rho, p = spearmanr(matched.loc[usable, 'SPINK5'], values[usable], alternative='two-sided')
        correlations.append({
            'cohort': 'TCGA-HNSC', 'population': population, 'n': int(usable.sum()),
            'markers_used': ','.join(available), 'markers_expected': len(genes), 'rho': rho, 'P': p,
        })
    scores_for_model = score_frame[list(signatures)].reindex(tc.index).rename(
        columns={'B lineage': 'B_lineage_score', 'Fibroblasts': 'Fibroblast_score'})
    score_model = tc.join(scores_for_model).dropna()
    if len(score_model) != len(tc):
        matched_cases = tc.loc[score_model.index]
        summarise(results, 'TCGA-HNSC', 'Primary model on composition-score cases', matched_cases, fit(matched_cases))
    summarise(results, 'TCGA-HNSC', 'Add B-lineage and fibroblast scores', score_model, fit(score_model))

    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(args.out / 'cox_sensitivity.tsv', sep='\t', index=False)
    pd.DataFrame(nonlinearity).to_csv(args.out / 'nonlinearity.tsv', sep='\t', index=False)
    pd.DataFrame(diagnostics).to_csv(args.out / 'ph_tests.tsv', sep='\t', index=False)
    pd.DataFrame(inclusion).to_csv(args.out / 'case_inclusion.tsv', sep='\t', index=False)
    pd.DataFrame(correlations).to_csv(args.out / 'composition_correlations.tsv', sep='\t', index=False)
    score_frame.to_csv(args.out / 'tcga_composition_scores.tsv', sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=Path(__file__).parent / 'inputs')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'results')
    parser.add_argument('--xena-expression', type=Path, help='Optional full Xena HiSeqV2 TSV gzip to regenerate the marker extract.')
    analyse(parser.parse_args())
