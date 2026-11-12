"""Paired cross-validation estimates and training-fold gate stability for GSE214265."""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

import kbridge_model as source

METRICS = {"AUROC": roc_auc_score, "AUPRC": average_precision_score}


def validation_splits(x, y):
    """Reproduce the original partition and all ten original repeated partitions."""
    yield "primary", StratifiedKFold(5, shuffle=True, random_state=42).split(x, y), 42
    for repeat in range(10):
        yield (f"repeat_{repeat + 1}",
               StratifiedKFold(5, shuffle=True, random_state=100 + repeat).split(x, y),
               1000 + 10 * repeat)


def fit_validation_folds(x, y):
    full_primary = np.empty(len(y))
    no_graph = np.empty(len(y))
    fold_weights = []
    fold_metrics = []
    for label, splits, start_seed in validation_splits(x, y):
        repeat_predictions = np.empty(len(y))
        for fold, (train, test) in enumerate(splits, start=1):
            random_state = start_seed + fold - 1
            full = source.KBridge(list(x.columns), seed=random_state).fit(x.iloc[train].to_numpy(), y[train])
            repeat_predictions[test] = full.predict_proba(x.iloc[test].to_numpy())[:, 1]
            fold_weights.append(pd.Series(full.gate, index=x.columns, name=f"{label}_fold_{fold}"))
            if label == "primary":
                variant = source.KBridge(list(x.columns), seed=random_state, use_graph=False)
                no_graph[test] = variant.fit(x.iloc[train].to_numpy(), y[train]).predict_proba(x.iloc[test].to_numpy())[:, 1]
        if label == "primary":
            full_primary = repeat_predictions.copy()
        fold_metrics.append({"partition": label, **{key: score(y, repeat_predictions)
                                                    for key, score in METRICS.items()}})
    return full_primary, no_graph, pd.DataFrame(fold_weights), pd.DataFrame(fold_metrics)


def paired_differences(y, predictions, n_boot=5000):
    """Resample identical cells across models, stratified by label."""
    rng = np.random.default_rng(42)
    class_indices = [np.flatnonzero(y == label) for label in np.unique(y)]
    boot_indices = [
        np.concatenate([rng.choice(indices, size=len(indices), replace=True)
                        for indices in class_indices])
        for _ in range(n_boot)
    ]
    scores = []
    for metric, score in METRICS.items():
        estimates = {key: score(y, pred) for key, pred in predictions.items()}
        values = {key: np.array([score(y[index], pred[index]) for index in boot_indices])
                  for key, pred in predictions.items()}
        for key in predictions:
            if key == "K-BRIDGE":
                continue
            delta = values["K-BRIDGE"] - values[key]
            low, high = np.quantile(delta, [0.025, 0.975])
            scores.append({"metric": metric, "comparator": key, "kbridge": estimates["K-BRIDGE"],
                           "estimate": estimates[key], "delta": estimates["K-BRIDGE"] - estimates[key],
                           "CI95_low": low, "CI95_high": high})
    return pd.DataFrame(scores)


def gate_stability(weights):
    pairwise = []
    ranks = weights.rank(axis=1, ascending=False, method="min")
    for left, right in itertools.combinations(weights.index, 2):
        correlation = spearmanr(weights.loc[left], weights.loc[right]).statistic
        left_top = set(ranks.columns[ranks.loc[left] <= 10])
        right_top = set(ranks.columns[ranks.loc[right] <= 10])
        overlap = len(left_top & right_top) / len(left_top | right_top)
        pairwise.append((left, right, correlation, overlap))
    stats = weights.agg(["mean", "std", "min", "max"]).T
    stats["top10_count"] = (ranks <= 10).sum()
    return stats.sort_values("mean", ascending=False), pd.DataFrame(
        pairwise, columns=["left_fold", "right_fold", "spearman", "top10_jaccard"]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, default=Path(__file__).with_name("internal_results") / "oof_predictions.csv")
    parser.add_argument("--output", type=Path, default=Path("validation_results"))
    options = parser.parse_args()
    source.DATA = options.data
    options.output.mkdir(parents=True, exist_ok=True)
    x, labels = source.load_human_training()
    y = labels.to_numpy()
    submitted = pd.read_csv(options.predictions, index_col=0)
    assert list(x.index) == list(submitted.index)
    assert np.array_equal(y, submitted["label"].to_numpy())
    full, no_graph, weights, fold_metrics = fit_validation_folds(x, y)
    np.testing.assert_allclose(full, submitted["K-BRIDGE"], atol=1e-11, rtol=1e-10)
    predictions = {key: value.to_numpy() for key, value in submitted.drop(columns="label").items()}
    predictions["No graph"] = no_graph
    paired = paired_differences(y, predictions)
    stability, pairwise = gate_stability(weights)
    paired.to_csv(options.output / "paired_differences.csv", index=False)
    weights.to_csv(options.output / "fold_gate_values.csv", index_label="fold")
    stability.to_csv(options.output / "gate_stability_by_gene.csv", index_label="gene")
    pairwise.to_csv(options.output / "gate_stability_between_folds.csv", index=False)
    fold_metrics.to_csv(options.output / "repeated_validation.csv", index=False)
    pd.DataFrame({**predictions, "label": y}, index=x.index).to_csv(
        options.output / "validation_predictions.csv", index_label="cell"
    )
    print(paired.to_string(index=False))
    print(stability.head(15).to_string())
    print(pairwise[["spearman", "top10_jaccard"]].describe().to_string())


if __name__ == "__main__":
    main()
