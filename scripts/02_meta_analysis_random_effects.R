###############################################################################
# Random-effects meta-analysis per gene with biological process as moderator
# Arce DP & Pratta GR (2026)
#
# For each gene, one multilevel model is fitted to all contrasts:
#   rma.mv(yi = logFC, V = SE^2, mods = ~ 0 + Grupo, random = ...)
# giving the pooled log2FC for stress, the pooled log2FC for ripening,
# and the ripening - stress difference with its p value.
#
# Random structure: separate heterogeneity per process at study and
# contrast level (struct = "DIAG"). Contrasts sharing a control or series
# (e.g. BR and RR vs MG of the same genotype) are grouped in "Study".
# If this model fails, a model with common heterogeneity
# (~ 1 | Study/Contrast) is used, and then a DerSimonian-Laird model.
#
# Input (in data_dir): one cleaned table per contrast with columns
#   GeneID, logFC, seLOGFC; and diseno.csv with columns
#   Archivo, Contrast, Study, Grupo ("estres" or "maduracion").
###############################################################################

suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
  library(parallel)
})

# ---------------------------------------------------------------- Settings
data_dir   <- "data/input_metafor"
dir_salida <- file.path(data_dir, "resultados_moderador")
dir.create(dir_salida, showWarnings = FALSE, recursive = TRUE)

n_cores                  <- max(1, detectCores() - 1)
fdr_cutoff               <- 0.05
lfc_cutoff               <- 1
min_studies              <- 2      # minimum distinct studies per process
incluir_GSE245512        <- FALSE  # pre-filtered list of 191 genes; excluded
heterogeneidad_por_grupo <- TRUE
solo_prueba              <- NA     # integer = random subset of genes for testing

# ------------------------------------------------------------ Data
diseno <- fread(file = file.path(data_dir, "diseno.csv"))
if (!incluir_GSE245512) diseno <- diseno[Study != "GSE245512"]

datos <- rbindlist(lapply(seq_len(nrow(diseno)), function(i) {
  d <- fread(file = file.path(data_dir, diseno$Archivo[i]))
  cat(sprintf("%-20s %-36s %6d genes\n", diseno$Contrast[i], diseno$Archivo[i], nrow(d)))
  d[, `:=`(Contrast = diseno$Contrast[i], Study = diseno$Study[i], Grupo = diseno$Grupo[i])]
}))
datos[, Grupo := factor(Grupo, levels = c("estres", "maduracion"))]

cobertura <- datos[, .(n_est = uniqueN(Study[Grupo == "estres"]),
                       n_mad = uniqueN(Study[Grupo == "maduracion"])), by = GeneID]
genes_ok <- cobertura[n_est >= min_studies & n_mad >= min_studies, GeneID]
cat(sprintf("\nTotal genes: %d | with >= %d studies in both processes: %d\n",
            nrow(cobertura), min_studies, length(genes_ok)))
if (!is.na(solo_prueba)) {
  set.seed(1)
  genes_ok <- sample(genes_ok, min(solo_prueba, length(genes_ok)))
}
datos <- datos[GeneID %in% genes_ok]

# ------------------------------------------------------ Per-gene models
ajustar <- function(g) {
  ctrl <- list(optimizer = "nlminb", iter.max = 1000)
  if (heterogeneidad_por_grupo) {
    res <- tryCatch(
      rma.mv(yi = logFC, V = seLOGFC^2, mods = ~ 0 + Grupo,
             random = list(~ Grupo | Study, ~ Grupo | Contrast), struct = c("DIAG", "DIAG"),
             data = g, method = "REML", control = ctrl),
      error = function(e) NULL)
    if (!is.null(res)) return(list(res = res, modelo = "rma.mv_tau2_por_grupo"))
  }
  res <- tryCatch(
    rma.mv(yi = logFC, V = seLOGFC^2, mods = ~ 0 + Grupo,
           random = ~ 1 | Study/Contrast, data = g, method = "REML", control = ctrl),
    error = function(e) NULL)
  if (!is.null(res)) return(list(res = res, modelo = "rma.mv_tau2_comun"))
  res <- tryCatch(rma(yi = logFC, vi = seLOGFC^2, mods = ~ 0 + Grupo, data = g, method = "DL"),
                  error = function(e) NULL)
  if (!is.null(res)) return(list(res = res, modelo = "rma_DL"))
  NULL
}

meta_gen <- function(g) {
  fit <- ajustar(g)
  if (is.null(fit)) return(NULL)
  res <- fit$res
  b <- as.numeric(res$b); V <- vcov(res)   # order: estres, maduracion
  dif <- b[2] - b[1]
  se_dif <- sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
  if (fit$modelo == "rma.mv_tau2_por_grupo") {
    tau2_e <- res$tau2[1] + res$gamma2[1]; tau2_m <- res$tau2[2] + res$gamma2[2]
  } else if (inherits(res, "rma.mv")) {
    tau2_e <- tau2_m <- sum(res$sigma2)
  } else tau2_e <- tau2_m <- res$tau2
  ge <- g[Grupo == "estres"]; gm <- g[Grupo == "maduracion"]
  data.table(
    GeneID = g$GeneID[1],
    logFC_estres = b[1], se_estres = res$se[1], p_estres = res$pval[1],
    logFC_maduracion = b[2], se_maduracion = res$se[2], p_maduracion = res$pval[2],
    Dif_mad_menos_estres = dif, se_dif = se_dif,
    p_dif = 2 * pnorm(-abs(dif / se_dif)),
    Tau2_estres = tau2_e, Tau2_maduracion = tau2_m,
    QE_pvalor = res$QEp,
    k_estres = nrow(ge), estudios_estres = uniqueN(ge$Study),
    k_maduracion = nrow(gm), estudios_maduracion = uniqueN(gm$Study),
    Concordancia_estres = mean(sign(ge$logFC) == sign(b[1])),
    Concordancia_maduracion = mean(sign(gm$logFC) == sign(b[2])),
    Modelo = fit$modelo
  )
}

por_gen <- split(datos, by = "GeneID", keep.by = TRUE)
cat(sprintf("Fitting %d genes on %d core(s)...\n", length(por_gen), n_cores))
t0 <- Sys.time()
res <- if (.Platform$OS.type == "windows" || n_cores == 1) {
  lapply(por_gen, function(g) suppressWarnings(meta_gen(g)))
} else {
  mclapply(por_gen, function(g) suppressWarnings(meta_gen(g)), mc.cores = n_cores)
}
res <- rbindlist(res)
cat(sprintf("Done: %d genes in %.1f min\n", nrow(res),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ------------------------------------------------------ FDR and classification
res[, `:=`(FDR_estres     = p.adjust(p_estres, "BH"),
           FDR_maduracion = p.adjust(p_maduracion, "BH"),
           FDR_dif        = p.adjust(p_dif, "BH"))]

sig_e <- res$FDR_estres < fdr_cutoff & abs(res$logFC_estres) >= lfc_cutoff
sig_m <- res$FDR_maduracion < fdr_cutoff & abs(res$logFC_maduracion) >= lfc_cutoff
res[, Clase := fcase(
  sig_e & sig_m & logFC_estres > 0 & logFC_maduracion > 0, "UP_comun",
  sig_e & sig_m & logFC_estres < 0 & logFC_maduracion < 0, "DOWN_comun",
  sig_e & sig_m,                                          "Discordante",
  sig_e & !sig_m,                                         "Solo_estres",
  !sig_e & sig_m,                                         "Solo_maduracion",
  default = "NS")]

# Genes ranked by the larger of the two FDR values
res[, rank_fdr := pmax(FDR_estres, FDR_maduracion)]
setorder(res, rank_fdr)

setcolorder(res, c("GeneID", "Clase", "logFC_estres", "FDR_estres", "logFC_maduracion",
                   "FDR_maduracion", "Dif_mad_menos_estres", "FDR_dif"))
fwrite(res, file.path(dir_salida, "resultados_metaanalisis_moderador.csv"))
fwrite(res[Clase == "UP_comun"],   file.path(dir_salida, "genes_UP_comunes_moderador.csv"))
fwrite(res[Clase == "DOWN_comun"], file.path(dir_salida, "genes_DOWN_comunes_moderador.csv"))

cat(sprintf("\nClassification (FDR < %.2f and |log2FC| >= %g):\n", fdr_cutoff, lfc_cutoff))
print(res[, .N, by = Clase][order(-N)])
cat(sprintf("\nGenes with different effect in stress and ripening (FDR_dif < %.2f): %d\n",
            fdr_cutoff, sum(res$FDR_dif < fdr_cutoff, na.rm = TRUE)))
cat("Models used:\n"); print(res[, .N, by = Modelo])

writeLines(capture.output(sessionInfo()), file.path(dir_salida, "sessionInfo_02_metafor.txt"))
