# Data

Place the input files in these folders before running the scripts.

| Folder | Contents |
| --- | --- |
| `input_meta_tomato/` | The 16 contrast tables (`#1_...csv` to `#16_...csv`) with columns `GeneID`, `logFC`, `pvalue`. Tables #4-#9 and #14 are produced by `scripts/00b_DESeq2_contrast_tables.R`. |
| `input_metafor/` | One cleaned table per contrast with columns `GeneID`, `logFC`, `seLOGFC`, and `diseno.csv` with columns `Archivo`, `Contrast`, `Study`, `Grupo` (`estres` or `maduracion`). |
| `own/` | Count tables of Caimanta, LA0722 and F1 per comparison (`toc_*_PvsV.csv`, `toc_*_RvsV.csv`), and `combined.csv` (27 samples) for script 03. |
| `public/` | `GSE174607_counts_anno.txt`, downloaded from GEO. |

Raw reads of Caimanta, LA0722 and F1: [BioProject accession]. Public datasets: GEO GSE294808, GSE245512, GSE152620, GSE151277, GSE174607 and GSE245759.
