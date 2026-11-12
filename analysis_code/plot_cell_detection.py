"""Plot patient-level SPINK5 detection fractions from the GSE103322 analysis."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

ORDER = ["Malignant epithelial", "B cells", "Other immune", "Stromal"]


def plot(input_path, output_path):
    data = pd.read_csv(input_path, sep="\t")
    assert not data.duplicated(["compartment", "patient"]).any()
    assert (data["size"] > 0).all()
    assert (data["sum"].between(0, data["size"])).all()
    data["fraction"] = data["sum"] / data["size"]
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42,
        "svg.fonttype": "none", "axes.spines.right": False,
        "axes.spines.top": False, "axes.spines.left": False,
        "axes.linewidth": 0.65, "figure.facecolor": "white",
    })
    fig, (ax, tally) = plt.subplots(
        1, 2, figsize=(6.7, 2.43), sharey=True,
        gridspec_kw={"width_ratios": [3.4, 1.15]},
    )
    fig.subplots_adjust(left=0.232, right=0.967, top=0.82, bottom=0.23, wspace=0.09)
    colors = ["#20617F", "#5A6F78", "#778791", "#778791"]
    for row, (compartment, color) in enumerate(zip(ORDER, colors)):
        group = data.loc[data["compartment"].eq(compartment)].sort_values("patient")
        assert len(group) > 0
        for fraction, tied in group.groupby("fraction", sort=True):
            offsets = np.linspace(-0.23, 0.23, len(tied)) if len(tied) > 1 else [0]
            ax.scatter(
                [fraction] * len(tied), row + np.asarray(offsets), s=13,
                color=color, edgecolor="white", linewidth=0.34, zorder=3,
            )
        tally.text(0.55, row, f"{group['sum'].sum():,} / {group['size'].sum():,}",
                   ha="center", va="center", color="#172A34")
    ax.set(yticks=list(range(len(ORDER))), yticklabels=ORDER, ylim=(3.5, -0.5),
           xlim=(-0.017, 0.67), xlabel="Cells with detectable SPINK5 per patient")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xticks([0, .2, .4, .6])
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.grid(axis="x", color="#DBE3E8", linewidth=0.7, zorder=0)
    tally.set_xlim(0, 1)
    tally.axis("off")
    tally.text(0.55, -0.68, "Detected / profiled\ncells", ha="center",
               va="bottom", fontweight="bold", color="#172A34")
    fig.text(0.23, 0.045, "Each point represents one patient within a compartment.",
             color="#45565F", fontsize=7)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="puram_patient_detection.tsv")
    parser.add_argument("output", type=Path, help="Vector PDF output")
    args = parser.parse_args()
    plot(args.input, args.output)
