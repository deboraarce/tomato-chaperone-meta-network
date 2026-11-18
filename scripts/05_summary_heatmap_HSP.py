"""Summary heatmap of HSP genes across contrasts.

Arce DP & Pratta GR (2026)

Input: HSP_forest_data.csv (from 04_forest_plots_HSP.py).
Cells are coloured by log2FC and marked with "*" when p < ALPHA; missing
values are grey. Genes are ordered by hierarchical clustering (average
linkage, Euclidean distance; missing values set to 0 for clustering only).

Requirements: pandas, numpy, scipy, matplotlib
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

INPUT_CSV = "HSP_forest_data.csv"
ALPHA = 0.05

STUDY_ORDER = [
    "GSE294808_heat-stress", "GSE294808_multi-stress", "GSE245512_et_HS",
    "GSE152620_STvsC", "GSE152620_SvsC",
    "GSE151277_D-5dvsD-1d", "GSE151277_heat24hvs0h",
    "GSE174607_Mic24hvs0h",
    "GSE245759_SHSvsCT", "GSE245759_HSvsCT",
    "own_cai_BRvsMG", "own_cai_RRvsMG",
    "own_PIMPI_LA0722_BRvsMG", "own_PIMPI_LA0722_RRvsMG",
    "own_F1_BRvsMG", "own_F1_RRvsMG",
]

OUTPUT_PNG = "HSP_summary_heatmap.png"
OUTPUT_PDF = "HSP_summary_heatmap.pdf"


def main():
    res = pd.read_csv(INPUT_CSV)

    pivot_fc = res.pivot(index="Gene", columns="Study", values="logFC").reindex(columns=STUDY_ORDER)
    pivot_p = res.pivot(index="Gene", columns="Study", values="pvalue").reindex(columns=STUDY_ORDER)

    Z = linkage(pivot_fc.fillna(0).values, method="average", metric="euclidean")
    gene_order = pivot_fc.index[leaves_list(Z)].tolist()

    pivot_fc = pivot_fc.loc[gene_order]
    pivot_p = pivot_p.loc[gene_order]

    fig, ax = plt.subplots(figsize=(9, 0.32 * len(gene_order) + 2))

    vmax = np.nanmax(np.abs(pivot_fc.values))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#e8e8e8")
    im = ax.imshow(np.ma.masked_invalid(pivot_fc.values), cmap=cmap,
                   vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(STUDY_ORDER)))
    ax.set_xticklabels(STUDY_ORDER, rotation=90, fontsize=8)
    ax.set_yticks(range(len(gene_order)))
    ax.set_yticklabels(gene_order, fontsize=7)

    for i in range(pivot_fc.shape[0]):
        for j in range(pivot_fc.shape[1]):
            p = pivot_p.values[i, j]
            if pd.notna(p) and p < ALPHA:
                ax.text(j, i, "*", ha="center", va="center",
                        color="black", fontsize=8, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("log2FC", fontsize=9)

    ax.set_title(f"Differential expression of HSP genes\n(* p < {ALPHA}; grey = no data)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    plt.savefig(OUTPUT_PDF, bbox_inches="tight")


if __name__ == "__main__":
    main()
