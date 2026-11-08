# Required for reproducible Kenward-Roger inference; do not silently fall back.
if (!requireNamespace("pbkrtest", quietly = TRUE))
  stop("Install pbkrtest before running this script; Kenward-Roger inference is required.")
emmeans::emm_options(lmer.df = "kenward-roger", disable.pbkrtest = FALSE)

# 40_TABELA7_estabilidade_contemporaneo.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(emmeans)
})
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

biomarkers <- list(
  NfL        = list(d1 = "NFLD1,2",        d5 = "NFLD5,2"),
  S100B      = list(d1 = "SB100D1,2",      d5 = "SB100D5,2"),
  `IL-6`     = list(d1 = "IL-6 (pg/mL)D1", d5 = "il6 deliriumd5"),
  `TNF-alfa` = list(d1 = "TNF (pg/mL)D1",  d5 = "tnfdeliriumd5")
)

build_long <- function(raw, d1col, d5col) {
  data.frame(id = rep(raw$ID, 2),
             time = factor(rep(c("D1","D5"), each = nrow(raw)), levels = c("D1","D5")),
             phenotype = factor(c(raw$fenotipo_D1, raw$fenotipo_D5),
                                 levels = c("G1_free","G2_delirium","G3_coma")),
             PRISMIV  = rep(raw[["PRISMIV%"]], 2),
             sedation = rep(raw[["sedation time"]], 2),
             y = c(raw[[d1col]], raw[[d5col]]))
}

model_specs <- list(
  M0_raw            = list(label = "Raw",               covars = character(0)),
  M1_PRISM          = list(label = "PRISM IV",           covars = "PRISMIV"),
  M2_sedation       = list(label = "Sedation",           covars = "sedation"),
  M3_PRISM_sedation = list(label = "PRISM IV + Sedation", covars = c("PRISMIV","sedation"))
)

extract_G2vG1 <- function(fit, tpoint) {
  emm <- emmeans(fit, ~ phenotype | time)
  con <- pairs(emm, adjust = "bonferroni")
  cdf <- as.data.frame(con)
  cdf <- cdf[cdf$time == tpoint, ]
  row <- cdf[cdf$contrast == "G1_free - G2_delirium", ]
  if (nrow(row) == 0) stop(sprintf("Contrast G1-G2 not found at %s", tpoint))
  beta <- -row$estimate[1]
  se   <- row$SE[1]
  p    <- row$p.value[1]
  ci   <- -confint(con)[confint(con)$time == tpoint & confint(con)$contrast == "G1_free - G2_delirium",
                          c("upper.CL","lower.CL")]
  data.frame(beta = round(beta,3), SE = round(se,3),
             CI_lo = round(as.numeric(ci[1,1]),3), CI_hi = round(as.numeric(ci[1,2]),3),
             p_bonferroni = round(p,4))
}

rows <- list()
for (bm in names(biomarkers)) {
  long <- build_long(raw, biomarkers[[bm]]$d1, biomarkers[[bm]]$d5)
  long <- long[!is.na(long$y) & !is.na(long$phenotype), ]
  long$logy <- log(long$y)

  for (mnm in names(model_specs)) {
    spec <- model_specs[[mnm]]
    covars <- if (length(spec$covars) > 0) paste("+", paste(spec$covars, collapse = " + ")) else ""
    form <- as.formula(sprintf("logy ~ phenotype * time %s + (1|id)", covars))
    fit <- lmer(form, data = long, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
    an <- anova(fit, type = 3, ddf = "Satterthwaite")

    d1c <- extract_G2vG1(fit, "D1")
    d5c <- extract_G2vG1(fit, "D5")

    rows[[length(rows)+1]] <- data.frame(
      Biomarker = bm, Model = spec$label,
      G2vG1_D1_beta = d1c$beta, G2vG1_D1_CI = sprintf("%.3f to %.3f", d1c$CI_lo, d1c$CI_hi),
      G2vG1_D1_p_bonf = d1c$p_bonferroni,
      G2vG1_D5_beta = d5c$beta, G2vG1_D5_CI = sprintf("%.3f to %.3f", d5c$CI_lo, d5c$CI_hi),
      G2vG1_D5_p_bonf = d5c$p_bonferroni,
      Phenotype_F = round(an["phenotype","F value"],2),
      Phenotype_df = sprintf("%.0f,%.1f", an["phenotype","NumDF"], an["phenotype","DenDF"]),
      Phenotype_p = round(an["phenotype","Pr(>F)"],4),
      Interaction_F = round(an["phenotype:time","F value"],2),
      Interaction_df = sprintf("%.0f,%.1f", an["phenotype:time","NumDF"], an["phenotype:time","DenDF"]),
      Interaction_p = round(an["phenotype:time","Pr(>F)"],4),
      N_observations = nobs(fit), N_patients = length(unique(long$id)),
      singular = isSingular(fit, tol = 1e-4)
    )
  }
}

tab7 <- do.call(rbind, rows)
write.csv(tab7, file.path(out_dir, "TABELA7_stability_CONTEMPORANEO.csv"), row.names = FALSE)
cat("=== Supplementary Table 6 - fenotipo contemporaneo, G2-vs-G1 por dia ===\n\n")
print(tab7, row.names = FALSE)
cat("\nNenhuma correcao adicional foi aplicada ENTRE os 4 modelos (robustez do mesmo contraste).\n")
cat("Bonferroni aplicado DENTRO de cada tempo x modelo, sobre as 3 comparacoes fenotipicas.\n")
