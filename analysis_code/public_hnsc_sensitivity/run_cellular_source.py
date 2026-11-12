"""Summarize deposited SPINK5 signal by cell compartment in GSE103322."""
import argparse
import csv
import gzip
from pathlib import Path

import pandas as pd

COMPARTMENTS = {
    'B cell': 'B cells', 'T cell': 'Other immune', 'Mast': 'Other immune',
    'Macrophage': 'Other immune', 'Dendritic': 'Other immune',
    'Fibroblast': 'Stromal', '-Fibroblast': 'Stromal',
    'Endothelial': 'Stromal', 'myocyte': 'Stromal',
}


def extract(matrix):
    wanted = {'Lymph node', 'classified  as cancer cell',
              'classified as non-cancer cells', 'non-cancer cell type', 'SPINK5'}
    with gzip.open(matrix, 'rt', newline='') as stream:
        reader = csv.reader(stream, delimiter='\t')
        names = next(reader)[1:]
        data = {}
        for label, *row in reader:
            label = label.strip("'")
            if label in wanted:
                assert len(names) == len(row) and label not in data
                data[label] = row
                if data.keys() == wanted:
                    break
    assert data.keys() == wanted and len(set(names)) == len(names)
    return pd.DataFrame(data, index=pd.Index(names, name='cell'))


def calculate(data, target):
    assert len(data) == 5902
    data['patient'] = data.index.to_series().str.extract(r'^HN(?:SCC)?_?(\d+)', expand=False)
    pooled = data.patient.isna()
    assert data.index[pooled].str.startswith('HNSCC_combo').all()
    data = data.loc[~pooled].copy()
    data['compartment'] = data['non-cancer cell type'].map(COMPARTMENTS)
    data.loc[data['classified  as cancer cell'].eq(1), 'compartment'] = 'Malignant epithelial'
    primary = data.loc[data['Lymph node'].eq(0) & data.compartment.notna()].copy()
    primary['detectable'] = primary.SPINK5.gt(0)
    per_patient = primary.groupby(['compartment', 'patient']).detectable.agg(['size', 'sum']).reset_index()
    per_patient['detected_fraction'] = per_patient['sum'] / per_patient['size']
    total = per_patient.groupby('compartment').agg(
        patients=('patient', 'nunique'), cells=('size', 'sum'),
        detected=('sum', 'sum'), median_patient_fraction=('detected_fraction', 'median'),
        min_patient_fraction=('detected_fraction', 'min'), max_patient_fraction=('detected_fraction', 'max'))
    total['positive_patients'] = per_patient.loc[per_patient['sum'].gt(0)].groupby('compartment').size()
    total['positive_patients'] = total['positive_patients'].fillna(0).astype(int)
    total['cell_fraction'] = total.detected / total.cells
    target.mkdir(parents=True, exist_ok=True)
    primary.to_csv(target / 'puram_primary_cell_expression.tsv', sep='\t')
    per_patient.to_csv(target / 'puram_patient_detection.tsv', sep='\t', index=False)
    total.to_csv(target / 'puram_summary.tsv', sep='\t')
    print(total.to_string())
    print('Distinct primary-tumour patient IDs:', primary.patient.nunique())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path)
    parser.add_argument('--input', type=Path, default=Path(__file__).parent / 'inputs' / 'puram_gene_extract.tsv')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'results')
    args = parser.parse_args()
    if args.matrix:
        extract(args.matrix).to_csv(args.input, sep='\t')
    calculate(pd.read_csv(args.input, sep='\t', index_col='cell'), args.out)
