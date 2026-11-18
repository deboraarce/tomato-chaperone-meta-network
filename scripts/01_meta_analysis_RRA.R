# ============================================================
# Meta-analysis of up-regulated genes by Robust Rank Aggregation (RRA)
# Arce DP & Pratta GR (2026)
#
# Input: 16 tables (one per control-vs-treatment contrast) with
#        columns GeneID, logFC and pvalue.
# Output: RRA scores, significant genes, volcano plot and heatmap
#         of the top 30 genes.
# ============================================================

library(RobustRankAggreg)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)

input_dir <- "data/input_meta_tomato"

# ------------------------------------------------------------
# 1. Load contrast tables
# ------------------------------------------------------------

read_contrast <- function(file) {
  read.csv(file.path(input_dir, file)) %>% rename(log2FC = logFC)
}

exp1  <- read_contrast("#1_GSE294808_heat-stress_log2FC.csv")
exp2  <- read_contrast("#2_GSE294808_multi-stress_log2FC.csv")
exp3  <- read_contrast("#3_GSE245512_et&HS_log2FC.csv")
exp4  <- read_contrast("#4_cai_BRvsMG.csv")
exp5  <- read_contrast("#5_cai_RRvsMG.csv")
exp6  <- read_contrast("#6_PIMPI_LA0722_BRvsMG.csv")
exp7  <- read_contrast("#7_PIMPI_LA0722_RRvsMG.csv")
exp8  <- read_contrast("#8_F1_BRvsMG.csv")
exp9  <- read_contrast("#9_F1_RRvsMG.csv")
exp10 <- read_contrast("#10_GSE152620_STvsC.csv")
exp11 <- read_contrast("#11_GSE152620_SvsC.csv")
exp12 <- read_contrast("#12_GSE151277_D-5d_vs_D-1d.csv")
exp13 <- read_contrast("#13_GSE151277_Tabla_Para_RRA_heat 24 hs_vs_heat 0 hs.csv")
exp14 <- read_contrast("#14_GSE174607_Tabla_Para_RRA_Mic_Tom_24h_vs_Mic_Tom_0h.csv")
exp15 <- read_contrast("#15_GSE245759_Tabla_Para_RRA_SHSvsCT.csv")
exp16 <- read_contrast("#16_GSE245759_Tabla_Para_RRA_HSvsCT.csv")

# ------------------------------------------------------------
# 2. Ranked lists of up-regulated genes (log2FC > 0, decreasing log2FC)
# ------------------------------------------------------------

up1  <- exp1  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up2  <- exp2  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up3  <- exp3  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up4  <- exp4  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up5  <- exp5  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up6  <- exp6  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up7  <- exp7  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up8  <- exp8  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up9  <- exp9  %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up10 <- exp10 %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up11 <- exp11 %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up12 <- exp12 %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up13 <- exp13 %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up14 <- exp14 %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up15 <- exp15 %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)
up16 <- exp16 %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(GeneID)

# ------------------------------------------------------------
# 3. Rank aggregation
# ------------------------------------------------------------

result_up <- aggregateRanks(
  glist = list(up1, up2, up3, up4, up5, up6, up7, up8,
               up9, up10, up11, up12, up13, up14, up15, up16),
  method = "RRA"
)
result_up$direction <- "up"
colnames(result_up) <- c("GeneID", "RRA_Score", "direction")

# ------------------------------------------------------------
# 4. Add log2FC per contrast, number of contrasts and mean log2FC
# ------------------------------------------------------------

fc_exp1  <- exp1[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp1  = log2FC)
fc_exp2  <- exp2[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp2  = log2FC)
fc_exp3  <- exp3[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp3  = log2FC)
fc_exp4  <- exp4[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp4  = log2FC)
fc_exp5  <- exp5[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp5  = log2FC)
fc_exp6  <- exp6[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp6  = log2FC)
fc_exp7  <- exp7[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp7  = log2FC)
fc_exp8  <- exp8[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp8  = log2FC)
fc_exp9  <- exp9[,  c("GeneID", "log2FC")] %>% rename(log2FC_exp9  = log2FC)
fc_exp10 <- exp10[, c("GeneID", "log2FC")] %>% rename(log2FC_exp10 = log2FC)
fc_exp11 <- exp11[, c("GeneID", "log2FC")] %>% rename(log2FC_exp11 = log2FC)
fc_exp12 <- exp12[, c("GeneID", "log2FC")] %>% rename(log2FC_exp12 = log2FC)
fc_exp13 <- exp13[, c("GeneID", "log2FC")] %>% rename(log2FC_exp13 = log2FC)
fc_exp14 <- exp14[, c("GeneID", "log2FC")] %>% rename(log2FC_exp14 = log2FC)
fc_exp15 <- exp15[, c("GeneID", "log2FC")] %>% rename(log2FC_exp15 = log2FC)
fc_exp16 <- exp16[, c("GeneID", "log2FC")] %>% rename(log2FC_exp16 = log2FC)

meta_data <- result_up %>%
  left_join(fc_exp1,  by = "GeneID") %>%
  left_join(fc_exp2,  by = "GeneID") %>%
  left_join(fc_exp3,  by = "GeneID") %>%
  left_join(fc_exp4,  by = "GeneID") %>%
  left_join(fc_exp5,  by = "GeneID") %>%
  left_join(fc_exp6,  by = "GeneID") %>%
  left_join(fc_exp7,  by = "GeneID") %>%
  left_join(fc_exp8,  by = "GeneID") %>%
  left_join(fc_exp9,  by = "GeneID") %>%
  left_join(fc_exp10, by = "GeneID") %>%
  left_join(fc_exp11, by = "GeneID") %>%
  left_join(fc_exp12, by = "GeneID") %>%
  left_join(fc_exp13, by = "GeneID") %>%
  left_join(fc_exp14, by = "GeneID") %>%
  left_join(fc_exp15, by = "GeneID") %>%
  left_join(fc_exp16, by = "GeneID") %>%
  mutate(
    exp_count = rowSums(!is.na(data.frame(
      log2FC_exp1, log2FC_exp2, log2FC_exp3, log2FC_exp4,
      log2FC_exp5, log2FC_exp6, log2FC_exp7, log2FC_exp8,
      log2FC_exp9, log2FC_exp10, log2FC_exp11, log2FC_exp12,
      log2FC_exp13, log2FC_exp14, log2FC_exp15, log2FC_exp16))),
    log2FC_mean = rowMeans(data.frame(
      log2FC_exp1, log2FC_exp2, log2FC_exp3, log2FC_exp4,
      log2FC_exp5, log2FC_exp6, log2FC_exp7, log2FC_exp8,
      log2FC_exp9, log2FC_exp10, log2FC_exp11, log2FC_exp12,
      log2FC_exp13, log2FC_exp14, log2FC_exp15, log2FC_exp16), na.rm = TRUE)
  ) %>%
  filter(exp_count >= 3) %>%
  select(GeneID, direction, RRA_Score, exp_count, log2FC_mean,
         log2FC_exp1, log2FC_exp2, log2FC_exp3, log2FC_exp4,
         log2FC_exp5, log2FC_exp6, log2FC_exp7, log2FC_exp8,
         log2FC_exp9, log2FC_exp10, log2FC_exp11, log2FC_exp12,
         log2FC_exp13, log2FC_exp14, log2FC_exp15, log2FC_exp16)

meta_data <- meta_data[order(meta_data$RRA_Score), ]

# ------------------------------------------------------------
# 5. Export results
# ------------------------------------------------------------

write.csv(meta_data, file = "meta_RRA_UP_results_tomato.csv", row.names = FALSE)

sig_up <- subset(meta_data, RRA_Score < 0.05)
write.csv(sig_up, file = "meta_RRA_UP_significant_tomato.csv", row.names = FALSE)

cat("Up-regulated genes present in >= 3 contrasts:", nrow(meta_data), "\n")
cat("Significant genes (RRA score < 0.05):", nrow(sig_up), "\n")

# ------------------------------------------------------------
# 6. Volcano plot
# ------------------------------------------------------------

volcano_data <- meta_data %>%
  mutate(
    neg_log10_score = -log10(RRA_Score),
    significance = case_when(
      RRA_Score < 0.05 ~ "Up-regulated (RRA < 0.05)",
      TRUE             ~ "Not significant"
    )
  )

top_labels <- volcano_data %>%
  filter(RRA_Score < 0.05) %>%
  arrange(RRA_Score) %>%
  slice_head(n = 10)

p_volcano <- ggplot(
  volcano_data,
  aes(x = log2FC_mean, y = neg_log10_score, color = significance)
) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up-regulated (RRA < 0.05)" = "#E74C3C",
               "Not significant"           = "grey70"),
    name = NULL
  ) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(
    data = top_labels,
    aes(label = GeneID),
    size = 2.8, color = "black",
    hjust = -0.1, vjust = 0.5,
    check_overlap = TRUE
  ) +
  labs(
    title    = "RRA - up-regulated genes",
    subtitle = paste0("Significant genes: ", nrow(sig_up)),
    x = "Mean log2FC",
    y = "-log10(RRA score)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey40"),
    legend.position = "top"
  )

ggsave("meta_RRA_UP_volcano_tomato.png", plot = p_volcano,
       width = 8, height = 6, dpi = 600)

# ------------------------------------------------------------
# 7. Heatmap of the top 30 significant genes
# Missing values (gene absent from a contrast) are set to 0
# for clustering only.
# ------------------------------------------------------------

if (nrow(sig_up) > 0) {
  heatmap_colors <- colorRampPalette(c("white", "#FADBD8", "#E74C3C", "#943126"))(100)

  heatmap_ready <- sig_up %>%
    arrange(RRA_Score) %>%
    head(30)

  heatmap_data <- heatmap_ready %>%
    select(GeneID,
           log2FC_exp1, log2FC_exp2, log2FC_exp3, log2FC_exp4,
           log2FC_exp5, log2FC_exp6, log2FC_exp7, log2FC_exp8,
           log2FC_exp9, log2FC_exp10, log2FC_exp11, log2FC_exp12,
           log2FC_exp13, log2FC_exp14, log2FC_exp15, log2FC_exp16) %>%
    column_to_rownames("GeneID") %>%
    as.matrix()

  colnames(heatmap_data) <- paste0("exp", 1:16)

  heatmap_data_cluster <- heatmap_data
  heatmap_data_cluster[is.na(heatmap_data_cluster)] <- 0

  pheatmap(
    heatmap_data_cluster,
    color         = heatmap_colors,
    cluster_cols  = FALSE,
    cluster_rows  = TRUE,
    show_rownames = TRUE,
    fontsize_row  = 7.5,
    fontsize_col  = 10,
    main          = "Top 30 up-regulated meta-genes (present in >= 3 contrasts)",
    filename      = "meta_RRA_UP_heatmap_tomato.png",
    width         = 8,
    height        = 9
  )
}

writeLines(capture.output(sessionInfo()), "sessionInfo_01_RRA.txt")
