# Required for reproducible Kenward-Roger inference; do not silently fall back.
if (!requireNamespace("pbkrtest", quietly = TRUE))
  stop("Install pbkrtest before running this script; Kenward-Roger inference is required.")
emmeans::emm_options(lmer.df = "kenward-roger", disable.pbkrtest = FALSE)

# 04_phase3_inflammatory_panel_FDR.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

panel <- list(
  list(d1 = "PCRD1", d5 = "PCRD5", label = "PCR", group = "Inflamatorio/hematologico"),
  list(d1 = "VHSD1", d5 = "VHSD5", label = "VHS", group = "Inflamatorio/hematologico"),
  list(d1 = "FERD1", d5 = "FERD5", label = "Ferritina", group = "Inflamatorio/hematologico"),
  list(d1 = "NLRD1", d5 = "NLRD5", label = "NLR", group = "Inflamatorio/hematologico"),
  list(d1 = "PLRD1", d5 = "PLRD5", label = "PLR", group = "Inflamatorio/hematologico"),
  list(d1 = "CRP/ABD1", d5 = "CRP/ABD5", label = "PCR/albumina", group = "Inflamatorio/hematologico"),
  list(d1 = "IL-2 (pg/mL)D1", d5 = "il-2 delirium d5", label = "IL-2", group = "Citocina"),
  list(d1 = "IL-4 (pg/mL)D1", d5 = "IL-4 (pg/mL)D5", label = "IL-4", group = "Citocina"),
  list(d1 = "IL-6 (pg/mL)D1", d5 = "il6 deliriumd5", label = "IL-6", group = "Citocina"),
  list(d1 = "IL-10 (pg/mL)D1", d5 = "IL-10 (pg/mL)D5", label = "IL-10", group = "Citocina"),
  list(d1 = "IL-17A (pg/mL)D1", d5 = "il-17Ad5", label = "IL-17A", group = "Citocina"),
  list(d1 = "IFN-γ (pg/mL)D1", d5 = "ifn-gamad5", label = "IFN-gamma", group = "Citocina"),
  list(d1 = "TNF (pg/mL)D1", d5 = "tnfdeliriumd5", label = "TNF-alfa", group = "Citocina"),
  list(d1 = "TEAC D1", d5 = "TEACD5", label = "TEAC", group = "Estresse oxidativo"),
  list(d1 = "TBARSD1", d5 = "TBARSD5", label = "TBARS", group = "Estresse oxidativo"),
  list(d1 = "GSHD1", d5 = "GSHD5", label = "GSH", group = "Estresse oxidativo")
)


results <- list()
contrast_rows <- list()

for (bm in panel) {
  d1v <- raw[[bm$d1]]; d5v <- raw[[bm$d5]]
  if (is.null(d1v) || is.null(d5v)) {
    cat(sprintf("  AVISO: coluna nao encontrada para %s (d1=%s d5=%s) - PULADO\n", bm$label, bm$d1, bm$d5))
    next
  }
  long <- data.frame(
    id = rep(raw$ID, 2),
    time = factor(rep(c("D1","D5"), each = nrow(raw)), levels = c("D1","D5")),
    phenotype = factor(c(raw$fenotipo_D1, raw$fenotipo_D5), levels = c("G1_free","G2_delirium","G3_coma")),
    y = c(d1v, d5v)
  )
  long <- long[!is.na(long$y) & !is.na(long$phenotype), ]
  long$logy <- log(long$y)

  fit <- tryCatch(
    lmer(logy ~ phenotype * time + (1|id), data = long, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("  ERRO ao ajustar %s: %s\n", bm$label, conditionMessage(fit)))
    next
  }
  singular <- isSingular(fit, tol = 1e-4)
  an <- tryCatch(anova(fit, type = 3), error = function(e) NULL)
  p_phenotype <- if (!is.null(an) && "phenotype" %in% rownames(an)) an["phenotype","Pr(>F)"] else NA
  p_interaction <- if (!is.null(an) && any(grepl(":", rownames(an)))) an[grep(":", rownames(an)),"Pr(>F)"] else NA

  results[[bm$label]] <- data.frame(
    biomarker = bm$label, grupo = bm$group, n_obs = nrow(long),
    singular = singular, p_phenotype_raw = p_phenotype, p_interaction_raw = p_interaction
  )

  emm <- tryCatch(emmeans(fit, ~ phenotype | time), error = function(e) NULL)
  if (!is.null(emm)) {
    pw <- as.data.frame(pairs(emm, adjust = "tukey"))
    pw$biomarker <- bm$label
    contrast_rows[[bm$label]] <- pw
  }
  cat(sprintf("  %-14s n=%3d singular=%-5s p(fenotipo)=%.4f p(interacao)=%.4f\n",
              bm$label, nrow(long), singular, p_phenotype, p_interaction))
}

res_df <- do.call(rbind, results)
rownames(res_df) <- NULL

res_df$q_FDR_phenotype <- p.adjust(res_df$p_phenotype_raw, method = "BH")
res_df$q_FDR_interaction <- p.adjust(res_df$p_interaction_raw, method = "BH")
res_df <- res_df[order(res_df$q_FDR_phenotype), ]

print(res_df[, c("biomarker","grupo","n_obs","singular","p_phenotype_raw","q_FDR_phenotype",
                  "p_interaction_raw","q_FDR_interaction")], row.names = FALSE, digits = 4)

n_sig_raw <- sum(res_df$p_phenotype_raw < 0.05, na.rm = TRUE)
n_sig_fdr <- sum(res_df$q_FDR_phenotype < 0.05, na.rm = TRUE)
cat(sprintf("\nBiomarcadores com p bruto <0.05 (efeito fenotipo): %d de %d\n", n_sig_raw, nrow(res_df)))
cat(sprintf("Biomarcadores com q-FDR <0.05 (efeito fenotipo):    %d de %d\n", n_sig_fdr, nrow(res_df)))

write.csv(res_df, file.path(out_dir, "table13_14_inflammatory_panel_FDR.csv"), row.names = FALSE)

contrasts_df <- do.call(rbind, contrast_rows)
write.csv(contrasts_df, file.path(out_dir, "table13_14_inflammatory_panel_contrasts.csv"), row.names = FALSE)

cat("\n=== FASE 3 (painel inflamatorio/oxidativo + FDR) CONCLUIDA ===\n")
cat("Saidas: table13_14_inflammatory_panel_FDR.csv, table13_14_inflammatory_panel_contrasts.csv\n")
