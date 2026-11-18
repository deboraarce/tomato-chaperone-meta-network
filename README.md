# Meta-analysis and co-expression network analysis of tomato stress and ripening transcriptomes

Code accompanying:

> Arce DP & Pratta GR (2026). [Meta-analysis as a bioinformatic and statistical framework for identifying key biomarkers in tomato plant breeding programs]. *Frontiers in Systems Biology*. 

Archived version of this repository:

The pipeline combines a systematic search of public RNA-Seq datasets, two meta-analysis methods (robust rank aggregation and a random-effects model), a weighted gene co-expression network (WGCNA) built from fruit transcriptomes of *Solanum lycopersicum* cv. Caimanta, *S. pimpinellifolium* LA0722 and their F1 hybrid at three ripening stages, and summary plots of heat shock protein genes.

## Repository structure

```
scripts/                      analysis scripts, run in numerical order
data/input_meta_tomato/       16 contrast tables (GeneID, logFC, pvalue) used by scripts 01, 04
data/input_metafor/           contrast tables with standard errors and diseno.csv, used by script 02
data/own/                     count tables of Caimanta, LA0722 and F1, used by scripts 00b and 03
data/public/                  files downloaded from GEO, used by script 00b
results/                      created by the scripts
```

## Scripts

Run from the repository root.

| Script | Analysis | Main input | Main output |
| --- | --- | --- | --- |
| `00_geo_systematic_search.py` | GEO DataSets search (NCBI Entrez) for tomato abiotic stress RNA-Seq studies, with preliminary classification of available processed files; followed by manual PRISMA screening | query in the script; `NCBI_EMAIL` environment variable | `results/systematic_review/*.csv`, `search_log.txt` |
| `00b_DESeq2_contrast_tables.R` | DESeq2 contrast tables for datasets without author-provided differential expression: Caimanta, LA0722 and F1 (BR vs MG and RR vs MG, per genotype) and GSE174607 (24 h vs 0 h) | `data/own/*.csv`, `data/public/GSE174607_counts_anno.txt` | tables #4-#9 and #14 in `data/input_meta_tomato/` |
| `01_meta_analysis_RRA.R` | Robust rank aggregation of up-regulated genes across 16 contrasts | `data/input_meta_tomato/*.csv` | RRA results, significant genes, volcano plot, heatmap |
| `02_meta_analysis_random_effects.R` | Per-gene random-effects meta-analysis (metafor `rma.mv`) with biological process (stress vs ripening) as moderator; 15 contrasts (GSE245512 excluded) | `data/input_metafor/` | `resultados_moderador/*.csv` |
| `03_DESeq2_kmeans_WGCNA.R` | DESeq2, k-means clustering and WGCNA of the 27 fruit samples (3 genotypes x 3 stages x 3 replicates); export of modules to Cytoscape | `data/own/combined.csv` | DEG tables, PCA, modules, module-trait correlations, Cytoscape files |
| `04_forest_plots_HSP.py` | Forest plots of 33 HSP genes across the 16 contrasts; SE and 95% CI reconstructed from p values | `data/input_meta_tomato/*.csv` | `HSP_forest_data.csv`, `HSP_forest_plots.pdf` |
| `05_summary_heatmap_HSP.py` | Summary heatmap of the forest plots | `HSP_forest_data.csv` | `HSP_summary_heatmap.png`, `.pdf` |

## Data sources and normalization

| Contrast tables | Source | Differential expression |
| --- | --- | --- |
| GSE294808, GSE245512, GSE152620, GSE245759 | GEO, results of the original authors | DESeq2 (authors) |
| GSE151277 | GEO, results of the original authors | Cuffdiff, Cufflinks v2.2.1 (authors) |
| GSE174607 | GEO raw counts | DESeq2, this study (script 00b) |
| Caimanta, LA0722, F1 | Cacchiarelli et al. (2024), [BioProject accession] | DESeq2, this study (script 00b) |

Public datasets were not re-normalized when the original authors provided differential expression results. Because robust rank aggregation uses only the order of genes within each contrast, it does not require effect sizes on a common scale.

## Main parameters

- DESeq2 contrasts (script 00b): design `~ condition`, 3 vs 3 replicates, genes with >= 10 total counts, median-of-ratios normalization, Wald test.
- RRA: genes ranked by decreasing log2FC within each contrast (log2FC > 0); genes present in >= 3 contrasts; significance at RRA score < 0.05.
- Random effects: REML, heterogeneity per process at study and contrast level; genes present in >= 2 studies per process; Benjamini-Hochberg FDR < 0.05 and |log2FC| >= 1.
- WGCNA (script 03): design `~ genotipo + etapa`; genes with >= 10 counts in >= 3 samples; DEGs at padj < 0.05 and |log2FC| > 1.2; signed network, soft-thresholding power 12, deepSplit 4, minModuleSize 30, mergeCutHeight 0.20; Cytoscape export with TOM threshold 0.05.

## Requirements

- R >= 4.x with DESeq2, WGCNA, RobustRankAggreg, metafor, data.table, dplyr, ggplot2, pheatmap, tibble, RColorBrewer, cluster.
- Python >= 3.9 with biopython, pandas, numpy, scipy, matplotlib.

Each R script writes a `sessionInfo_*.txt` file with the R and package versions used.

## Citation

If you use this code, please cite the article above and the Zenodo archive (see `CITATION.cff`).

## License

MIT (see `LICENSE`).
