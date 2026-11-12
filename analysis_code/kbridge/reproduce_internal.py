"""Recompute five-fold internal classification, ablations and optional permutations."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import kbridge_model as source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name("data"))
    parser.add_argument("--output", type=Path, default=Path("internal_results"))
    parser.add_argument("--permutations", action="store_true", help="Also rerun 100 label permutations.")
    options = parser.parse_args()
    source.DATA = options.data
    options.output.mkdir(parents=True, exist_ok=True)
    x, y = source.load_human_training()
    x.join(y.rename("label")).to_csv(options.output / "training_panel_log_expression.csv", index_label="cell")
    metrics, predictions, splits = source.evaluate_models(x, y)
    metrics.to_csv(options.output / "model_metrics.csv", index=False)
    predictions.to_csv(options.output / "oof_predictions.csv", index_label="cell")
    source.evaluate_ablation(x, y, splits).to_csv(options.output / "ablation.csv", index=False)
    if options.permutations:
        repeated, permutation = source.seed_and_permutation_tests(x, y)
        repeated.to_csv(options.output / "seed_sensitivity.csv", index=False)
        pd.DataFrame({"AUROC": permutation}).to_csv(options.output / "permutation.csv", index=False)
        observed = metrics.loc[metrics.model.eq("K-BRIDGE"), "AUROC"].item()
        print("Permutation p:", (np.count_nonzero(permutation >= observed) + 1) / (len(permutation) + 1))


if __name__ == "__main__":
    main()
