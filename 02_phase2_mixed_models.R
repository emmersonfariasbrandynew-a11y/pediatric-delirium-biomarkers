# Required for reproducible Kenward-Roger inference; do not silently fall back.
if (!requireNamespace("pbkrtest", quietly = TRUE))
  stop("Install pbkrtest before running this script; Kenward-Roger inference is required.")
emmeans::emm_options(lmer.df = "kenward-roger", disable.pbkrtest = FALSE)

# 02_phase2_mixed_models.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(performance)
})

options(stringsAsFactors = FALSE)
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

long <- data.frame(
  id        = rep(raw$ID, 2),
  time      = factor(rep(c("D1","D5"), each = nrow(raw)), levels = c("D1","D5")),
  phenotype = c(raw$fenotipo_D1, raw$fenotipo_D5),
  NfL       = c(raw[["NFLD1,2"]], raw[["NFLD5,2"]]),
  S100B     = c(raw[["SB100D1,2"]], raw[["SB100D5,2"]]),
  PRISMIV   = rep(raw[["PRISMIV%"]], 2),
  sedation  = rep(raw[["sedation time"]], 2),
  sepsis    = rep(raw$sepsis, 2)
)
long$phenotype <- factor(long$phenotype, levels = c("G1_free","G2_delirium","G3_coma"))
long$logNfL   <- log(long$NfL)
long$logS100B <- log(long$S100B)

cat(sprintf("Long dataset: %d rows (%d patients x 2 timepoints)\n", nrow(long), nrow(raw)))
cat(sprintf("Missing phenotype: %d | Missing NfL: %d | Missing S100B: %d\n\n",
            sum(is.na(long$phenotype)), sum(is.na(long$NfL)), sum(is.na(long$S100B))))


describe_biomarker <- function(df, var) {
  agg <- aggregate(df[[var]], by = list(time = df$time, phenotype = df$phenotype),
                    FUN = function(x) {
                      x <- x[!is.na(x)]
                      c(n = length(x), median = median(x), IQR_lo = quantile(x,.25),
                        IQR_hi = quantile(x,.75), min = min(x), max = max(x))
                    })
  res <- do.call(data.frame, agg)
  names(res) <- c("time","phenotype","n","median","IQR_lo","IQR_hi","min","max")
  res
}
nfl_desc <- describe_biomarker(long, "NfL")
s100b_desc <- describe_biomarker(long, "S100B")
cat("NfL (ng/mL) por fenotipo contemporaneo e momento:\n"); print(nfl_desc)
cat("\nS100B (ng/mL) por fenotipo contemporaneo e momento:\n"); print(s100b_desc)
write.csv(nfl_desc, file.path(out_dir, "table_NfL_distribution.csv"), row.names = FALSE)
write.csv(s100b_desc, file.path(out_dir, "table_S100B_distribution.csv"), row.names = FALSE)

cat("\nShapiro-Wilk (log-escala, checagem de normalidade dos residuos brutos):\n")
cat(sprintf("  log(NfL):   W=%.4f p=%.4f\n", shapiro.test(long$logNfL)$statistic, shapiro.test(long$logNfL)$p.value))
cat(sprintf("  log(S100B): W=%.4f p=%.4f\n\n", shapiro.test(long$logS100B)$statistic, shapiro.test(long$logS100B)$p.value))


fit_models <- function(df, yvar, label) {
  cat(sprintf("---- %s ----\n", label))
  f0 <- as.formula(sprintf("%s ~ phenotype * time + (1|id)", yvar))
  f1 <- as.formula(sprintf("%s ~ phenotype * time + PRISMIV + (1|id)", yvar))
  f2 <- as.formula(sprintf("%s ~ phenotype * time + sedation + (1|id)", yvar))
  f3 <- as.formula(sprintf("%s ~ phenotype * time + PRISMIV + sedation + (1|id)", yvar))

  models <- list(M0_raw = f0, M1_PRISM = f1, M2_sedation = f2, M3_PRISM_sedation = f3)
  fits <- list()
  for (nm in names(models)) {
    fit <- tryCatch(
      lmer(models[[nm]], data = df, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) e, warning = function(w) w)
    fits[[nm]] <- fit
    conv_ok  <- !inherits(fit, "error") && !inherits(fit, "warning")
    singular <- if (!inherits(fit, "error")) isSingular(fit, tol = 1e-4) else NA
    cat(sprintf("  %-20s convergiu=%s  singular=%s\n", nm,
                ifelse(conv_ok || inherits(fit,"lmerModLmerTest"), "sim", "verificar"),
                as.character(singular)))
  }
  fits
}

nfl_fits   <- fit_models(long, "logNfL", "NfL (log)")
cat("\n")
s100b_fits <- fit_models(long, "logS100B", "S100B (log)")
cat("\n")

saveRDS(list(nfl_fits = nfl_fits, s100b_fits = s100b_fits, long = long),
        file.path(out_dir, "phase2_models.rds"))

cat("=== RESUMO MODELO PARCIMONIOSO (M0) - NfL ===\n")
print(summary(nfl_fits$M0_raw))
cat("\n=== RESUMO MODELO PARCIMONIOSO (M0) - S100B ===\n")
print(summary(s100b_fits$M0_raw))

for (bm in c("NfL","S100B")) {
  fit <- if (bm == "NfL") nfl_fits$M0_raw else s100b_fits$M0_raw
  cat(sprintf("\n-- %s --\n", bm))
  vc <- as.data.frame(VarCorr(fit))
  cat("Variancia do intercepto aleatorio:\n"); print(vc)
  cat(sprintf("Singular: %s\n", isSingular(fit, tol = 1e-4)))
  png(file.path(out_dir, sprintf("diag_residuals_%s.png", bm)), width = 900, height = 450)
  par(mfrow = c(1,2))
  plot(fitted(fit), resid(fit), main = paste(bm, "- Residuos vs ajustados"),
       xlab = "Ajustado", ylab = "Residuo"); abline(h = 0, lty = 2)
  qqnorm(resid(fit), main = paste(bm, "- QQ residuos")); qqline(resid(fit))
  dev.off()
}

for (bm in c("NfL","S100B")) {
  fit <- if (bm == "NfL") nfl_fits$M0_raw else s100b_fits$M0_raw
  cat(sprintf("\n-- %s : contrastes entre fenotipos, dentro de cada tempo --\n", bm))
  emm1 <- emmeans(fit, ~ phenotype | time)
  print(pairs(emm1, adjust = "tukey"))
  cat(sprintf("\n-- %s : mudanca D5-D1 dentro de cada fenotipo --\n", bm))
  emm2 <- emmeans(fit, ~ time | phenotype)
  print(pairs(emm2, adjust = "tukey"))
  write.csv(as.data.frame(pairs(emm1, adjust = "tukey")),
            file.path(out_dir, sprintf("contrasts_%s_phenotype_within_time.csv", bm)), row.names = FALSE)
  write.csv(as.data.frame(pairs(emm2, adjust = "tukey")),
            file.path(out_dir, sprintf("contrasts_%s_time_within_phenotype.csv", bm)), row.names = FALSE)
}

cat("\n=== FASE 2 (modelos NfL/S100B) CONCLUIDA ===\n")
cat("Saidas em output/: table_NfL_distribution.csv, table_S100B_distribution.csv,\n")
cat("diag_residuals_NfL.png, diag_residuals_S100B.png, contrasts_*.csv,\n")
cat("phase2_models.rds (todos os 8 modelos ajustados)\n")
