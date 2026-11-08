# Required for reproducible Kenward-Roger inference; do not silently fall back.
if (!requireNamespace("pbkrtest", quietly = TRUE))
  stop("Install pbkrtest before running this script; Kenward-Roger inference is required.")
emmeans::emm_options(lmer.df = "kenward-roger", disable.pbkrtest = FALSE)

# 15_FIGURA2_PAINELB_forest_plot_cohens_d.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(emmeans)
})
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw  <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
mods <- readRDS(file.path(out_dir, "phase2_models.rds"))   # M0_raw/M1_PRISM/M2_sedation/M3_PRISM_sedation

sigma_total_lmer <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  tau00 <- sum(vc$vcov[vc$grp == "id" & is.na(vc$var2)])
  sqrt(sigma(fit)^2 + tau00)
}

std_from_fit <- function(fit, label, biomarker) {
  emm <- emmeans(fit, ~ phenotype)
  ed  <- tryCatch(df.residual(fit), error = function(e) NA)
  if (is.na(ed) || !is.finite(ed)) ed <- Inf
  eff <- eff_size(emm, sigma = sigma_total_lmer(fit), edf = ed)
  d <- as.data.frame(eff)
  data.frame(model = label, biomarker = biomarker, contrast = d$contrast,
             cohens_d = d$effect.size, lower95 = d$lower.CL, upper95 = d$upper.CL,
             SE = d$SE, edf_used = ed)
}

rows_14 <- list()
for (bm in c("NfL", "S100B")) {
  fits <- if (bm == "NfL") mods$nfl_fits else mods$s100b_fits
  rows_14[[paste0(bm, "_Raw")]]       <- std_from_fit(fits$M0_raw,            "Raw",                    bm)
  rows_14[[paste0(bm, "_PRISM")]]     <- std_from_fit(fits$M1_PRISM,          "Adj. PRISM IV",          bm)
  rows_14[[paste0(bm, "_Sedation")]]  <- std_from_fit(fits$M2_sedation,       "Adj. Sedation",          bm)
  rows_14[[paste0(bm, "_PRISMSed")]]  <- std_from_fit(fits$M3_PRISM_sedation, "Adj. PRISMIV+Sedation",  bm)
}

build_long <- function(raw, d1col, d5col) {
  data.frame(id = rep(raw$ID, 2),
             time = factor(rep(c("D1","D5"), each = nrow(raw)), levels = c("D1","D5")),
             phenotype = factor(c(raw$fenotipo_D1, raw$fenotipo_D5), levels = c("G1_free","G2_delirium","G3_coma")),
             PRISMIV = rep(raw[["PRISMIV%"]], 2), sedation = rep(raw[["sedation time"]], 2),
             sepsis = rep(raw$sepsis, 2), y = c(raw[[d1col]], raw[[d5col]]))
}
bm_cols <- list(NfL = list(d1="NFLD1,2", d5="NFLD5,2"), S100B = list(d1="SB100D1,2", d5="SB100D5,2"))

rows_56 <- list()
for (bm in names(bm_cols)) {
  long <- build_long(raw, bm_cols[[bm]]$d1, bm_cols[[bm]]$d5)
  long <- long[!is.na(long$y) & !is.na(long$phenotype), ]
  long$logy <- log(long$y)

  fit_plus_sepsis <- lmer(logy ~ phenotype * time + sepsis + (1|id), data = long,
                           REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  rows_56[[paste0(bm, "_PlusSepsis")]] <- std_from_fit(fit_plus_sepsis, "Sepsis-adjusted", bm)

}

std_from_lm <- function(fit, data_for_emm, label, biomarker) {
  emm <- emmeans(fit, ~ phenotype_D1, data = data_for_emm)
  eff <- eff_size(emm, sigma = sigma(fit), edf = df.residual(fit))
  d <- as.data.frame(eff)
  data.frame(model = label, biomarker = biomarker, contrast = d$contrast,
             cohens_d = d$effect.size, lower95 = d$lower.CL, upper95 = d$upper.CL,
             SE = d$SE, edf_used = df.residual(fit))
}

rows_78 <- list()
for (bm in names(bm_cols)) {
  d1v <- raw[[bm_cols[[bm]]$d1]]; d5v <- raw[[bm_cols[[bm]]$d5]]
  delta <- log(d5v) - log(d1v)
  df_delta <- data.frame(id = raw$ID, phenotype_D1 = raw$fenotipo_D1, sepsis = raw$sepsis, delta = delta)
  df_delta <- df_delta[!is.na(df_delta$delta) & !is.na(df_delta$phenotype_D1), ]
  df_delta$phenotype_D1 <- factor(df_delta$phenotype_D1, levels = c("G1_free","G2_delirium","G3_coma"))

  fit_delta <- lm(delta ~ phenotype_D1, data = df_delta)
  rows_78[[paste0(bm, "_DeltaCohort")]] <- std_from_lm(fit_delta, df_delta, "Delta (D5-D1) Cohort", bm)

}

forest_df <- do.call(rbind, c(rows_14, rows_56, rows_78))
rownames(forest_df) <- NULL

recode_contrast <- function(x) {
  x <- gsub("phenotype_D1", "", x); x <- gsub("phenotype", "", x)
  x <- gsub("G1_free", "G1", x); x <- gsub("G2_delirium", "G2", x); x <- gsub("G3_coma", "G3", x)
  x
}
forest_df$contrast_label <- recode_contrast(forest_df$contrast)

model_order <- c("Raw","Adj. PRISM IV","Adj. PRISMIV+Sedation","Adj. Sedation",
                  "Sepsis-adjusted",
                  "Delta (D5-D1) Cohort")
forest_df$model <- factor(forest_df$model, levels = rev(model_order))

write.csv(forest_df, file.path(out_dir, "FIGURA2_PAINELB_forest_data_cohens_d.csv"), row.names = FALSE)
cat("=== FIGURA 2 - PAINEL B: Cohen's d (standardized mean difference via emmeans::eff_size) ===\n\n")
print(forest_df[order(forest_df$biomarker, forest_df$model), c("biomarker","model","contrast_label","cohens_d","lower95","upper95")],
      row.names = FALSE)
cat("\nSalvo em output/FIGURA2_PAINELB_forest_data_cohens_d.csv\n")
