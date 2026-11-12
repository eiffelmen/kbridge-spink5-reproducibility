"""External FACS-labelled human tonsil validation of the frozen B-to-ASC score."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

import kbridge_model as source


def load_tonsil(path, file_id, genes):
    counts = pd.read_csv(path / f"{file_id}_filtered_endogenes_UMI.csv.gz", index_col=0)
    meta = pd.read_excel(path / f"{file_id}_filtered_metadata.xlsx",
                         sheet_name=1, index_col="UniqueCellID")
    assert counts.index.is_unique and meta.index.is_unique
    assert counts.index.equals(meta.index), "Inconsistent expression and phenotype cell identifiers"
    selected = meta.loc[:, ["Donor", "Plate", "PhenoOriginal"]].copy()
    gene_names = counts.columns.str.upper()
    present = sorted(set(genes) & set(gene_names))
    keep = counts.columns[gene_names.isin(present)]
    denom = counts.sum(axis=1)
    assert denom.gt(0).all()
    panel = counts.loc[:, keep]
    panel.columns = panel.columns.str.upper()
    panel = panel.T.groupby(level=0).sum().T
    log_cpm = np.log1p(panel.div(denom, axis=0) * 1e6)
    return selected, log_cpm


def apply_models(log_cpm, model, baseline):
    genes = model.genes
    std = log_cpm.std(axis=0, ddof=1).replace(0, 1)
    target_scaled = (log_cpm - log_cpm.mean(axis=0)).div(std, axis=1).reindex(
        columns=genes, fill_value=0
    )
    source_center = pd.Series(model.scaler.mean_, index=genes)
    source_scaled_input = log_cpm.reindex(columns=genes).fillna(source_center)
    baseline_mean = pd.Series(baseline[0].mean_, index=genes)
    baseline_input = log_cpm.reindex(columns=genes).fillna(baseline_mean)
    scores = pd.DataFrame(index=log_cpm.index)
    scores["source_standardized_kbridge"] = model.predict_proba(source_scaled_input.to_numpy())[:, 1]
    scores["within_tonsil_kbridge"] = model.project_z(target_scaled.to_numpy())
    scores["source_standardized_elastic_net"] = baseline.predict_proba(baseline_input.to_numpy())[:, 1]
    scores["within_tonsil_elastic_net"] = baseline[1].decision_function(target_scaled.to_numpy())
    return scores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--external", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("validation_results"))
    options = parser.parse_args()
    options.output.mkdir(parents=True, exist_ok=True)
    source.DATA = options.data
    x, y = source.load_human_training()
    model = source.KBridge(list(x.columns)).fit(x.to_numpy(), y.to_numpy())
    baseline = source.make_pipeline(
        source.StandardScaler(), source.LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.5,
            class_weight="balanced", max_iter=8000, random_state=42
        )
    ).fit(x.to_numpy(), y.to_numpy())
    outputs = []
    input_panels = []
    coverages = []
    for id in ["GSM4074307_Tonsil1Bcells", "GSM4074308_Tonsil2Bcells"]:
        meta, log_cpm = load_tonsil(options.external, id, x.columns)
        meta["accession"] = id.split("_")[0]
        outputs.append(meta.join(apply_models(log_cpm, model, baseline)))
        input_panels.append(meta.join(log_cpm))
        coverages.append({"accession": id.split("_")[0], "shared_genes": len(log_cpm.columns),
                          "missing_genes": ";".join(sorted(set(x.columns) - set(log_cpm.columns)))})
    full = pd.concat(outputs)
    assert full.index.is_unique
    groups = full.groupby("accession")
    rows = []
    for label, group in list(groups) + [("pooled", full)]:
        for col in group.columns[4:]:
            assert set(group["PhenoOriginal"]) == {"Mem", "GC", "PB/PC"}
            y = group["PhenoOriginal"].eq("PB/PC")
            rows.append({"sample": label, "model_scaling": col, "cells": len(y),
                         "B_cells": int((~y).sum()), "plasmablasts_and_plasma_cells": int(y.sum()),
                         "AUROC": roc_auc_score(y, group[col]),
                         "AUPRC": average_precision_score(y, group[col])})
    pd.concat(input_panels).to_csv(options.output / "tonsil_input_log_cpm.csv", index_label="cell")
    full.to_csv(options.output / "tonsil_predictions.csv", index_label="cell")
    pd.DataFrame(coverages).to_csv(options.output / "tonsil_gene_coverage.csv", index=False)
    summary = pd.DataFrame(rows)
    summary.to_csv(options.output / "tonsil_performance.csv", index=False)
    print(full.groupby(["accession", "Donor", "PhenoOriginal"]).size().to_string())
    print(pd.DataFrame(coverages).to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
