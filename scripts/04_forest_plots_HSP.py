"""Forest plots of HSP genes across differential expression contrasts.

Arce DP & Pratta GR (2026)

Input: one table per contrast with columns GeneID, logFC, pvalue (in this
order). The standard error of each log2FC is reconstructed from its
two-sided p value, assuming a Wald test:
    z = qnorm(1 - p/2);  SE = |log2FC| / z;  95% CI = log2FC +/- 1.96 SE
p values are bounded to [1e-300, 0.999999]. Genes absent from a contrast
are reported as missing.

Output:
    HSP_forest_data.csv   long table (gene, contrast, log2FC, p, SE, 95% CI)
    HSP_forest_plots.pdf  one forest plot per gene

Requirements: pandas, numpy, scipy, matplotlib
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# =========================== Settings ===========================

INPUT_DIR = "data/input_meta_tomato"

FILES = {
    "GSE294808_heat-stress": "#1_GSE294808_heat-stress_log2FC.csv",
    "GSE294808_multi-stress": "#2_GSE294808_multi-stress_log2FC.csv",
    "GSE245512_et_HS": "#3_GSE245512_et&HS_log2FC.csv",
    "own_cai_BRvsMG": "#4_cai_BRvsMG.csv",
    "own_cai_RRvsMG": "#5_cai_RRvsMG.csv",
    "own_PIMPI_LA0722_BRvsMG": "#6_PIMPI_LA0722_BRvsMG.csv",
    "own_PIMPI_LA0722_RRvsMG": "#7_PIMPI_LA0722_RRvsMG.csv",
    "own_F1_BRvsMG": "#8_F1_BRvsMG.csv",
    "own_F1_RRvsMG": "#9_F1_RRvsMG.csv",
    "GSE152620_STvsC": "#10_GSE152620_STvsC.csv",
    "GSE152620_SvsC": "#11_GSE152620_SvsC.csv",
    "GSE151277_D-5dvsD-1d": "#12_GSE151277_D-5d_vs_D-1d.csv",
    "GSE151277_heat24hvs0h": "#13_GSE151277_Tabla_Para_RRA_heat 24 hs_vs_heat 0 hs.csv",
    "GSE174607_Mic24hvs0h": "#14_GSE174607_Tabla_Para_RRA_Mic_Tom_24h_vs_Mic_Tom_0h.csv",
    "GSE245759_SHSvsCT": "#15_GSE245759_Tabla_Para_RRA_SHSvsCT.csv",
    "GSE245759_HSvsCT": "#16_GSE245759_Tabla_Para_RRA_HSvsCT.csv",
}

STUDY_ORDER = list(FILES.keys())

GENES = [
    "Solyc01g009200", "Solyc01g009220", "Solyc01g098790", "Solyc01g098810",
    "Solyc01g102960", "Solyc02g080410", "Solyc02g093600", "Solyc03g082420",
    "Solyc03g113930", "Solyc03g123540", "Solyc04g014480", "Solyc04g072250",
    "Solyc04g082720", "Solyc04g082740", "Solyc05g014280", "Solyc06g076520",
    "Solyc06g076540", "Solyc06g076560", "Solyc06g076570", "Solyc07g064020",
    "Solyc08g062340", "Solyc08g062450", "Solyc08g078700", "Solyc08g078710",
    "Solyc08g078720", "Solyc09g007140", "Solyc09g011710", "Solyc09g015000",
    "Solyc09g015020", "Solyc10g086680", "Solyc11g020330", "Solyc11g071560",
    "Solyc12g042830",
]

OUTPUT_CSV = "HSP_forest_data.csv"
OUTPUT_PDF = "HSP_forest_plots.pdf"
ALPHA = 0.05

# =========================== Data ===========================

def load_tables(folder, files_dict):
    dfs = {}
    for name, filename in files_dict.items():
        df = pd.read_csv(os.path.join(folder, filename))
        df.columns = [c.strip() for c in df.columns]
        df = df.rename(columns={df.columns[0]: "GeneID",
                                df.columns[1]: "logFC",
                                df.columns[2]: "pvalue"})
        df["GeneID"] = df["GeneID"].astype(str).str.strip()
        n_dup = df["GeneID"].duplicated().sum()
        if n_dup > 0:
            print(f"{name}: {n_dup} duplicated GeneID, first occurrence kept")
            df = df.drop_duplicates(subset="GeneID", keep="first")
        dfs[name] = df.set_index("GeneID")[["logFC", "pvalue"]]
        print(f"{name}: {df.shape[0]} genes")
    return dfs


def se_and_ci(logfc, pvalue):
    p_clip = min(max(pvalue, 1e-300), 0.999999)
    z = norm.ppf(1 - p_clip / 2)
    if z <= 0 or np.isnan(z):
        return np.nan, np.nan, np.nan
    se = abs(logfc) / z
    return se, logfc - 1.96 * se, logfc + 1.96 * se


def build_long_table(dfs, genes, study_order):
    rows = []
    for gene in genes:
        for st in study_order:
            d = dfs[st]
            if gene in d.index:
                logfc = float(d.loc[gene, "logFC"])
                p = float(d.loc[gene, "pvalue"])
                se, lo, hi = se_and_ci(logfc, p)
                rows.append([gene, st, logfc, p, se, lo, hi])
            else:
                rows.append([gene, st, np.nan, np.nan, np.nan, np.nan, np.nan])
    return pd.DataFrame(rows, columns=["Gene", "Study", "logFC", "pvalue",
                                       "SE", "CI_low", "CI_high"])

# =========================== Plots ===========================

def plot_forest(res, study_order, pdf_file, alpha=0.05):
    genes = res["Gene"].unique().tolist()
    with PdfPages(pdf_file) as pdf:
        for gene in genes:
            sub = res[res["Gene"] == gene].set_index("Study").reindex(study_order)
            n = len(study_order)
            fig, ax = plt.subplots(figsize=(7.5, 0.35 * n + 1.5))
            ypos = np.arange(n)[::-1]

            for i, st in zip(ypos, study_order):
                row = sub.loc[st]
                if pd.isna(row["logFC"]):
                    ax.text(0, i, "no data", va="center", ha="center",
                            fontsize=8, color="gray", style="italic")
                    continue
                logfc, lo, hi, p = row["logFC"], row["CI_low"], row["CI_high"], row["pvalue"]
                color = "#c0392b" if p < alpha else "#7f8c8d"
                facecolor = color if p < alpha else "white"
                ax.hlines(i, lo, hi, color=color, lw=1.5, zorder=2)
                ax.plot(logfc, i, marker="o", markersize=6,
                        markerfacecolor=facecolor, markeredgecolor=color, zorder=3)

            ax.axvline(0, color="black", lw=0.8, linestyle="--", zorder=1)
            ax.set_yticks(ypos)
            ax.set_yticklabels(study_order, fontsize=8)
            ax.set_xlabel("log2FC (95% CI reconstructed from p value)", fontsize=9)
            ax.set_title(gene, fontsize=12, fontweight="bold")
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(axis="x", labelsize=8)

            legend_elems = [
                Line2D([0], [0], marker="o", color="#c0392b",
                       markerfacecolor="#c0392b", label=f"p < {alpha}", lw=1.5),
                Line2D([0], [0], marker="o", color="#7f8c8d",
                       markerfacecolor="white", label=f"p >= {alpha}", lw=1.5),
            ]
            ax.legend(handles=legend_elems, loc="lower right", fontsize=7, frameon=False)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    print(f"{pdf_file}: {len(genes)} pages")

# =========================== Main ===========================

if __name__ == "__main__":
    dfs = load_tables(INPUT_DIR, FILES)
    res = build_long_table(dfs, GENES, STUDY_ORDER)
    res.to_csv(OUTPUT_CSV, index=False)
    print(f"{OUTPUT_CSV}: {res.shape[0]} rows")
    plot_forest(res, STUDY_ORDER, OUTPUT_PDF, alpha=ALPHA)
