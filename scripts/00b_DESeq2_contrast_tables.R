# ============================================================
# Differential expression tables for meta-analysis input (DESeq2)
# Arce DP & Pratta GR (2026)
#
# Used for contrasts without author-provided differential expression:
#   - Caimanta, LA0722 and F1: BR vs MG and RR vs MG, per genotype
#     (tables #4-#9)
#   - GSE174607 (Micro-Tom fruit, heat stress): 24 h vs 0 h, from the
#     raw counts deposited in GEO (table #14)
#
# Each contrast: 3 vs 3 biological replicates, design ~ condition,
# genes with >= 10 total counts, median-of-ratios normalization and
# Wald test (DESeq2 defaults).
#
# Output per contrast: GeneID, logFC, lfcSE, pvalue, padj
# ============================================================

library(DESeq2)

output_dir <- "data/input_meta_tomato"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

run_contrast <- function(counts_file, control_cols, treated_cols,
                         control_label, treated_label, output_file) {
  counts <- if (grepl("\\.csv$", counts_file)) {
    read.csv(counts_file, row.names = 1, check.names = FALSE)
  } else {
    read.delim(counts_file, row.names = 1, check.names = FALSE)
  }
  counts <- round(as.matrix(counts[, c(control_cols, treated_cols)]))

  coldata <- data.frame(
    row.names = colnames(counts),
    condition = factor(c(rep(control_label, length(control_cols)),
                         rep(treated_label, length(treated_cols))),
                       levels = c(control_label, treated_label))
  )

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = coldata,
                                design    = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", treated_label, control_label))

  out <- data.frame(
    GeneID = sub("\\.[0-9]+$", "", rownames(res)),
    logFC  = res$log2FoldChange,
    lfcSE  = res$lfcSE,
    pvalue = res$pvalue,
    padj   = res$padj
  )
  write.csv(out, file.path(output_dir, output_file), row.names = FALSE)
  invisible(out)
}

# ------------------------------------------------------------
# Own data: one count table per genotype and comparison
# (columns: 3 BR or RR samples, then 3 MG samples)
# ------------------------------------------------------------

own_contrasts <- data.frame(
  counts_file = c("data/own/toc_Cai_PvsV.csv", "data/own/toc_Cai_RvsV.csv",
                  "data/own/toc_P_PvsV.csv",   "data/own/toc_P_RvsV.csv",
                  "data/own/toc_F1_PvsV.csv",  "data/own/toc_F1_RvsV.csv"),
  treated     = c("BR", "RR", "BR", "RR", "BR", "RR"),
  output_file = c("#4_cai_BRvsMG.csv",          "#5_cai_RRvsMG.csv",
                  "#6_PIMPI_LA0722_BRvsMG.csv", "#7_PIMPI_LA0722_RRvsMG.csv",
                  "#8_F1_BRvsMG.csv",           "#9_F1_RRvsMG.csv")
)

for (i in seq_len(nrow(own_contrasts))) {
  run_contrast(
    counts_file   = own_contrasts$counts_file[i],
    control_cols  = 4:6,
    treated_cols  = 1:3,
    control_label = "MG",
    treated_label = own_contrasts$treated[i],
    output_file   = own_contrasts$output_file[i]
  )
}

# ------------------------------------------------------------
# GSE174607: raw counts from GEO (GSE174607_counts_anno.txt)
# ------------------------------------------------------------

run_contrast(
  counts_file   = "data/public/GSE174607_counts_anno.txt",
  control_cols  = c("Mic_Tom_0h_1",  "Mic_Tom_0h_2",  "Mic_Tom_0h_3"),
  treated_cols  = c("Mic_Tom_24h_1", "Mic_Tom_24h_2", "Mic_Tom_24h_3"),
  control_label = "h0",
  treated_label = "h24",
  output_file   = "#14_GSE174607_Tabla_Para_RRA_Mic_Tom_24h_vs_Mic_Tom_0h.csv"
)

writeLines(capture.output(sessionInfo()), "sessionInfo_00b_DESeq2.txt")
