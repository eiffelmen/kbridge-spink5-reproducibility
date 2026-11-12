"""Reproduce the patient-level effects and test details reported in the manuscript."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, wilcoxon

HERE = Path(__file__).resolve().parent
SEED = 20261112
N_BOOTSTRAP = 10000


def rank_correlation_bootstrap(values, indices):
    """Return Spearman coefficients for two scores against SPINK5 by resample."""
    chunks = []
    for chunk in np.array_split(indices, 10):
        rank = rankdata(values[chunk], axis=1)
        rank -= rank.mean(axis=1, keepdims=True)
        cross_product = (rank[:, :, :1] * rank[:, :, 1:]).sum(axis=1)
        sum_squares = np.square(rank).sum(axis=1)
        chunks.append(cross_product / np.sqrt(sum_squares[:, :1] * sum_squares[:, 1:]))
    return np.concatenate(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paired", type=Path,
                        default=HERE.parent / "public_hnsc/results/paired_tumour_normal.tsv")
    parser.add_argument("--composition", type=Path,
                        default=HERE / "results/tcga_composition_scores.tsv")
    parser.add_argument("--output", type=Path, default=HERE / "results/statistical_reporting.tsv")
    args = parser.parse_args()

    paired = pd.read_csv(args.paired, sep="\t")
    assert paired["patient"].is_unique
    difference = paired["tumour_minus_normal"].to_numpy()
    assert np.isfinite(difference).all()
    assert np.all(difference != 0) and len(np.unique(np.abs(difference))) == len(difference)
    np.testing.assert_allclose(difference, paired["SPINK5_log2_tumour"] - paired["SPINK5_log2_normal"])
    signed_rank = wilcoxon(difference, alternative="two-sided", method="exact")
    rng = np.random.default_rng(SEED)
    pair_draws = rng.integers(0, len(difference), size=(N_BOOTSTRAP, len(difference)))
    paired_ci = np.quantile(np.median(difference[pair_draws], axis=1), [0.025, 0.975])
    result = [{
        "comparison": "Paired primary tumour minus matched solid normal",
        "n": len(difference),
        "unit": "patient pair",
        "effect": "Median log2 expression difference",
        "estimate": np.median(difference),
        "CI95_low": paired_ci[0],
        "CI95_high": paired_ci[1],
        "test": "Wilcoxon signed rank (exact)",
        "statistic": "W",
        "statistic_value": signed_rank.statistic,
        "degrees_of_freedom": None,
        "P_two_sided": signed_rank.pvalue,
    }]

    scores = pd.read_csv(args.composition, sep="\t", index_col="sampleID")
    score_names = ["B lineage", "Fibroblasts"]
    scores = scores[["SPINK5_log2", *score_names]].dropna()
    assert scores.index.is_unique
    n = len(scores)
    rng = np.random.default_rng(SEED)
    patient_draws = rng.integers(0, n, size=(N_BOOTSTRAP, n))
    correlations = rank_correlation_bootstrap(scores.to_numpy(), patient_draws)
    ci = np.quantile(correlations, [0.025, 0.975], axis=0)
    for index, name in enumerate(score_names):
        estimate = spearmanr(scores["SPINK5_log2"], scores[name], alternative="two-sided")
        degrees = n - 2
        t_statistic = estimate.statistic * np.sqrt(degrees / (1 - estimate.statistic**2))
        result.append({
            "comparison": f"Bulk SPINK5 versus {name} score",
            "n": n,
            "unit": "primary tumour patient",
            "effect": "Spearman rho",
            "estimate": estimate.statistic,
            "CI95_low": ci[0, index],
            "CI95_high": ci[1, index],
            "test": "Spearman rank correlation (asymptotic)",
            "statistic": "t",
            "statistic_value": t_statistic,
            "degrees_of_freedom": degrees,
            "P_two_sided": estimate.pvalue,
        })

    table = pd.DataFrame(result)
    table["degrees_of_freedom"] = table["degrees_of_freedom"].astype("Int64")
    table["bootstrap_resamples"] = N_BOOTSTRAP
    table["bootstrap_seed"] = SEED
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
