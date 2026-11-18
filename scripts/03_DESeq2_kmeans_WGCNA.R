#------------------------------------------------------------------------------
# Differential expression, k-means clustering and WGCNA co-expression analysis
# of Caimanta (Cai), S. pimpinellifolium LA0722 (P) and F1 fruit at three
# ripening stages
# Arce DP & Pratta GR (2026)
#
# Design: 3 ripening stages (V = mature green, P = breaker, R = red ripe)
#         x 3 genotypes (Cai, F, P) x 3 biological replicates = 27 samples
# Sample names: <stage><replicate>_<genotype>, e.g. "V1_Cai", "P2_F", "R3_P"
#------------------------------------------------------------------------------

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tibble)
library(cluster)
library(WGCNA)

options(stringsAsFactors = FALSE)

counts_file <- "data/own/combined.csv"

#------------------------------------------------------------------------------
# 1. Count matrix and sample information
#------------------------------------------------------------------------------
counts <- read.csv(counts_file, row.names = 1)
counts <- round(as.matrix(counts))

sample_names <- colnames(counts)
etapa    <- substr(sample_names, 1, 1)
genotipo <- sub("^[VPR][0-9]_", "", sample_names)

colData <- data.frame(
  row.names = sample_names,
  etapa    = factor(etapa,    levels = c("V", "P", "R")),
  genotipo = factor(genotipo, levels = c("Cai", "F", "P"))
)
print(table(colData$etapa, colData$genotipo))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = colData,
  design    = ~ genotipo + etapa
)

# Keep genes with >= 10 counts in >= 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
resultsNames(dds)

#------------------------------------------------------------------------------
# 2. Differential expression by ripening stage, controlling for genotype
#------------------------------------------------------------------------------
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ genotipo)
res_etapa_global <- results(dds_lrt)
res_etapa_global <- res_etapa_global[order(res_etapa_global$padj), ]
write.csv(as.data.frame(res_etapa_global), "resultados_efecto_global_etapa.csv")

res_P_vs_V <- results(dds, contrast = c("etapa", "P", "V"))
res_R_vs_V <- results(dds, contrast = c("etapa", "R", "V"))
res_R_vs_P <- results(dds, contrast = c("etapa", "R", "P"))

res_P_vs_V <- res_P_vs_V[order(res_P_vs_V$padj), ]
res_R_vs_V <- res_R_vs_V[order(res_R_vs_V$padj), ]
res_R_vs_P <- res_R_vs_P[order(res_R_vs_P$padj), ]

write.csv(as.data.frame(res_P_vs_V), "resultados_Pinton_vs_Verde.csv")
write.csv(as.data.frame(res_R_vs_V), "resultados_Rojo_vs_Verde.csv")
write.csv(as.data.frame(res_R_vs_P), "resultados_Rojo_vs_Pinton.csv")

summary(res_P_vs_V, alpha = 0.05)
summary(res_R_vs_V, alpha = 0.05)
summary(res_R_vs_P, alpha = 0.05)

#------------------------------------------------------------------------------
# 3. Variance-stabilized data and PCA
#------------------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("etapa", "genotipo"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = etapa, shape = genotipo)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("V" = "forestgreen", "P" = "orange", "R" = "red")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA - ripening stage and genotype")
ggsave("PCA_etapa_genotipo.png", plot = p, width = 7, height = 5, dpi = 300)

counts_norm <- counts(dds, normalized = TRUE)
write.csv(counts_norm, "cuentas_normalizadas.csv")

#------------------------------------------------------------------------------
# 4. DEGs: padj < 0.05 and |log2FC| > 1.2 in BR vs MG or RR vs MG
#------------------------------------------------------------------------------
sig_genes_R_vs_V <- res_R_vs_V[which(res_R_vs_V$padj < 0.05 & abs(res_R_vs_V$log2FoldChange) > 1.2), ]
sig_genes_P_vs_V <- res_P_vs_V[which(res_P_vs_V$padj < 0.05 & abs(res_P_vs_V$log2FoldChange) > 1.2), ]

degs_union <- unique(c(rownames(sig_genes_R_vs_V), rownames(sig_genes_P_vs_V)))
cat("DEGs used for clustering and WGCNA:", length(degs_union), "\n")

vst_mat    <- assay(vsd)[degs_union, ]
vst_scaled <- t(scale(t(vst_mat)))

#------------------------------------------------------------------------------
# 5. K-means clustering of DEGs
#------------------------------------------------------------------------------
wss <- vector()
sil <- vector()
k_range <- 2:10
dist_vst <- dist(vst_scaled)

set.seed(123)
for (i in seq_along(k_range)) {
  k <- k_range[i]
  km_tmp <- kmeans(vst_scaled, centers = k, nstart = 25, iter.max = 50)
  wss[i] <- km_tmp$tot.withinss
  sil[i] <- summary(silhouette(km_tmp$cluster, dist_vst))$avg.width
}

png("Kmeans_evaluation.png", width = 800, height = 400)
par(mfrow = c(1, 2))
plot(k_range, wss, type = "b", pch = 19, xlab = "Number of clusters (k)", ylab = "WSS", main = "Elbow method")
plot(k_range, sil, type = "b", pch = 19, xlab = "Number of clusters (k)", ylab = "Mean silhouette", main = "Silhouette method")
dev.off()

k_optimo <- 4
set.seed(123)
km_res <- kmeans(vst_scaled, centers = k_optimo, nstart = 25, iter.max = 50)

df_clusters <- data.frame(
  GeneID  = names(km_res$cluster),
  Cluster = factor(km_res$cluster)
)
write.csv(df_clusters, "kmeans_gene_clusters.csv", row.names = FALSE)

annot_cols <- as.data.frame(colData(dds)[, c("etapa", "genotipo")])
rownames(annot_cols) <- colnames(vst_scaled)

annot_rows <- data.frame(Cluster = factor(km_res$cluster, levels = 1:k_optimo))
rownames(annot_rows) <- names(km_res$cluster)

vst_ordered <- vst_scaled[order(km_res$cluster), ]

ann_colors <- list(
  etapa    = c("V" = "forestgreen", "P" = "orange", "R" = "red"),
  genotipo = c("Cai" = "grey85", "F" = "grey55", "P" = "grey25"),
  Cluster  = c("1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3", "4" = "#FF7F00")
)

pheatmap(
  vst_ordered,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_row    = annot_rows,
  annotation_col    = annot_cols,
  annotation_colors = ann_colors,
  filename          = "Heatmap_Kmeans_DEGs.png",
  main              = "K-means clustering of DEGs"
)

#------------------------------------------------------------------------------
# 6. WGCNA
#------------------------------------------------------------------------------
enableWGCNAThreads()

datExpr <- t(assay(vsd)[degs_union, ])

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Scale-free topology fit across soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
write.csv(sft$fitIndices, "WGCNA_soft_threshold_fit.csv", row.names = FALSE)

png("WGCNA_soft_thresholding.png", width = 800, height = 400)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "blue")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# Network construction and module detection
soft_power <- 12

net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  deepSplit         = 4,
  minModuleSize     = 30,
  mergeCutHeight    = 0.20,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  verbose           = 3
)

moduleColors <- labels2colors(net$colors)

png("WGCNA_gene_dendrogram.png", width = 800, height = 600)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

#------------------------------------------------------------------------------
# 7. Module-trait relationships
# Trait columns are renamed explicitly because stage "P" (breaker) and
# genotype "P" (S. pimpinellifolium) would otherwise share the same name.
#------------------------------------------------------------------------------
traits <- model.matrix(~ 0 + etapa + genotipo, data = as.data.frame(colData(dds)))

nombres_esperados <- c("etapaV", "etapaP", "etapaR", "genotipoF", "genotipoP")
if (!identical(colnames(traits), nombres_esperados)) {
  stop("Unexpected trait column order: ", paste(colnames(traits), collapse = ", "))
}
colnames(traits) <- c("Verde", "Pinton", "Rojo", "F", "Pimpinellifolium")

MEs <- orderMEs(net$MEs)
moduleTraitCor    <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

write.csv(moduleTraitCor,    "WGCNA_module_trait_cor.csv")
write.csv(moduleTraitPvalue, "WGCNA_module_trait_pvalue.csv")

png("WGCNA_module_trait_relationships.png", width = 900, height = 650)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(8, 8.5, 4, 3))
labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = colnames(traits),
  yLabels       = names(MEs),
  ySymbols      = names(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.7,
  zlim          = c(-1, 1),
  main          = "Module - stage/genotype correlation"
)

n_etapa <- 3
n_total <- ncol(traits)
abline(v = n_etapa + 0.5, lwd = 2, col = "black")
mtext("RIPENING STAGE", side = 3, line = 2.2, at = (1 + n_etapa) / 2, cex = 1, font = 2)
mtext("GENOTYPE (vs Cai)", side = 3, line = 2.2, at = (n_etapa + 1 + n_total) / 2, cex = 1, font = 2)
dev.off()

res_wgcna <- data.frame(
  GeneID       = colnames(datExpr),
  WGCNA_Module = moduleColors
)
df_final_clusters <- merge(df_clusters, res_wgcna, by = "GeneID", all = TRUE)
write.csv(df_final_clusters, "resumen_clustering_kmeans_wgcna.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# 8. Export of selected modules to Cytoscape
#------------------------------------------------------------------------------
df_log2fc <- data.frame(
  GeneID        = rownames(res_R_vs_V),
  log2FC_R_vs_V = res_R_vs_V$log2FoldChange,
  padj_R_vs_V   = res_R_vs_V$padj,
  log2FC_P_vs_V = res_P_vs_V[rownames(res_R_vs_V), ]$log2FoldChange,
  padj_P_vs_V   = res_P_vs_V[rownames(res_R_vs_V), ]$padj
)

node_attributes_master <- df_final_clusters %>%
  left_join(df_log2fc, by = "GeneID")
write.csv(node_attributes_master, "Cytoscape_NODES_master_attributes.csv", row.names = FALSE)

df_equivalencia <- data.frame(
  Numero_Modulo = paste0("M", net$colors),
  Color_Modulo  = moduleColors
) %>% unique() %>% arrange(Numero_Modulo)
write.csv(df_equivalencia, "equivalencia_modulos_numero_color.csv", row.names = FALSE)

modulos_de_interes <- c("blue", "yellow", "green", "salmon", "cyan")
umbral_tom <- 0.05

dir.create("Cytoscape_Input_Files", showWarnings = FALSE)

# TOM computed within each module of interest
for (mod in modulos_de_interes) {
  genes_mod <- colnames(datExpr)[moduleColors == mod]
  if (length(genes_mod) < 2) next

  subTOM <- TOMsimilarityFromExpr(datExpr[, genes_mod], power = soft_power, networkType = "signed")
  colnames(subTOM) <- genes_mod
  rownames(subTOM) <- genes_mod

  exportNetworkToCytoscape(
    subTOM,
    edgeFile     = file.path("Cytoscape_Input_Files", paste0("Cytoscape_EDGES_Module_", mod, ".txt")),
    nodeFile     = file.path("Cytoscape_Input_Files", paste0("Cytoscape_NODES_Module_", mod, ".txt")),
    weighted     = TRUE,
    threshold    = umbral_tom,
    nodeNames    = genes_mod,
    altNodeNames = genes_mod,
    nodeAttr     = moduleColors[moduleColors == mod]
  )
}

writeLines(capture.output(sessionInfo()), "sessionInfo_03_WGCNA.txt")
