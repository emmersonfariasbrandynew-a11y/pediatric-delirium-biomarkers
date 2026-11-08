# Required for reproducible Kenward-Roger inference; do not silently fall back.
if (!requireNamespace("pbkrtest", quietly = TRUE))
  stop("Install pbkrtest before running this script; Kenward-Roger inference is required.")
emmeans::emm_options(lmer.df = "kenward-roger", disable.pbkrtest = FALSE)

# 32_TABELA23_mortality_LMM_full.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw$PRISMIV <- raw[["PRISMIV%"]]
cat(sprintf("Loaded %d patients\n", nrow(raw)))

biomarkers <- list(
  "CRP (mg/L)"      = list(d1 = "PCRD1",             d5 = "PCRD5"),
  "ESR (mm/h)"      = list(d1 = "VHSD1",             d5 = "VHSD5"),
  "Ferritin (ng/mL)"= list(d1 = "FERD1",             d5 = "FERD5"),
  "NLR"             = list(d1 = "NLRD1",             d5 = "NLRD5"),
  "PLR"             = list(d1 = "PLRD1",             d5 = "PLRD5"),
  "CRP/Albumin"     = list(d1 = "CRP/ABD1",          d5 = "CRP/ABD5"),
  "VIS"             = list(d1 = "VISd1",             d5 = "VISD5"),
  "IL-2 (pg/mL)"    = list(d1 = "IL-2 (pg/mL)D1",    d5 = "il-2 delirium d5"),
  "IL-4 (pg/mL)"    = list(d1 = "IL-4 (pg/mL)D1",    d5 = "IL-4 (pg/mL)D5"),
  "IL-6 (pg/mL)"    = list(d1 = "IL-6 (pg/mL)D1",    d5 = "il6 deliriumd5"),
  "IL-10 (pg/mL)"   = list(d1 = "IL-10 (pg/mL)D1",   d5 = "IL-10 (pg/mL)D5"),
  "IL-17A (pg/mL)"  = list(d1 = "IL-17A (pg/mL)D1",  d5 = "il-17Ad5"),
  "IFN-gamma (pg/mL)" = list(d1 = "IFN-γ (pg/mL)D1", d5 = "ifn-gamad5"),
  "TNF-alfa (pg/mL)"= list(d1 = "TNF (pg/mL)D1",     d5 = "tnfdeliriumd5"),
  "TEAC (mmol/L)"   = list(d1 = "TEAC D1",           d5 = "TEACD5"),
  "TBARS (umol/L)"  = list(d1 = "TBARSD1",           d5 = "TBARSD5"),
  "GSH (umol/mL)"   = list(d1 = "GSHD1",             d5 = "GSHD5"),
  "S100B (ng/mL)"   = list(d1 = "SB100D1,2",         d5 = "SB100D5,2"),
  "NfL (ng/mL)"     = list(d1 = "NFLD1,2",           d5 = "NFLD5,2")
)

covariate_sets <- list(
  "Unadjusted"          = character(0),
  "PRISM-IV"             = c("PRISMIV"),
  "PRISM-IV + Sepsis"    = c("PRISMIV", "sepsis")
)

build_long <- function(raw, d1col, d5col) {
  data.frame(
    id = rep(raw$ID, 2),
    logy = c(log(raw[[d1col]]), log(raw[[d5col]])),
    time = factor(rep(c("D1","D5"), each = nrow(raw)), levels = c("D1","D5")),
    mortality = factor(rep(ifelse(raw$obitooverall == 1, "NonSurvivor", "Survivor"), 2),
                        levels = c("Survivor","NonSurvivor")),
    PRISMIV = rep(raw$PRISMIV, 2),
    sepsis = rep(raw$sepsis, 2)
  )
}

fmt_p <- function(p) if (is.na(p)) NA else if (p < 0.001) "<0.001" else sprintf("%.3f", p)

results <- list()
for (bm in names(biomarkers)) {
  d1col <- biomarkers[[bm]]$d1; d5col <- biomarkers[[bm]]$d5
  if (!(d1col %in% names(raw)) || !(d5col %in% names(raw))) {
    cat(sprintf("SKIP %s: column not found (%s / %s)\n", bm, d1col, d5col))
    next
  }
  long <- build_long(raw, d1col, d5col)
  long <- long[!is.na(long$logy) & is.finite(long$logy) & !is.na(long$mortality), ]

  for (model_name in names(covariate_sets)) {
    covs <- covariate_sets[[model_name]]
    fstr <- paste("logy ~ mortality*time", if (length(covs)) paste("+", paste(covs, collapse=" + ")) else "", "+ (1|id)")
    f <- as.formula(fstr)
    fit_warnings <- character(0)
    fit <- withCallingHandlers(
      tryCatch(
        lmer(f, data = long, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
        error = function(e) NULL),
      warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (is.null(fit)) { cat(sprintf("FIT FAILED: %s / %s\n", bm, model_name)); next }
    a <- tryCatch(anova(fit, type = 3), error = function(e) NULL)
    if (is.null(a)) { cat(sprintf("ANOVA FAILED: %s / %s\n", bm, model_name)); next }

    n_obs <- nobs(fit)
    n_patients <- length(unique(long$id))
    singular <- isSingular(fit, tol = 1e-4)
    convergence_flag <- if (singular) "singular fit" else if (length(fit_warnings) > 0) paste(fit_warnings, collapse = "; ") else "OK"
    if (singular || length(fit_warnings) > 0) {
      cat(sprintf("FLAG %s / %s: singular=%s warnings=[%s]\n", bm, model_name, singular, paste(fit_warnings, collapse="; ")))
    }

    mort_p <- a["mortality", "Pr(>F)"]
    time_p <- a["time", "Pr(>F)"]
    inter_p <- a["mortality:time", "Pr(>F)"]

    emm_mort <- tryCatch(emmeans(fit, ~ mortality), error = function(e) NULL)
    pw_mort_p <- NA; pw_mort_dir <- "NA"
    if (!is.null(emm_mort)) {
      pw <- summary(pairs(emm_mort, adjust = "bonferroni"))
      pw_mort_p <- pw$p.value[1]
      est <- pw$estimate[1]
      pw_mort_dir <- if (pw_mort_p < 0.05) {
        if (est < 0) "NonSurvivor > Survivor" else "Survivor > NonSurvivor"
      } else "Survivor ~= NonSurvivor"
    }

    emm_time <- tryCatch(emmeans(fit, ~ time), error = function(e) NULL)
    pw_time_p <- NA; pw_time_dir <- "NA"
    if (!is.null(emm_time)) {
      pwt <- summary(pairs(emm_time, adjust = "bonferroni"))
      pw_time_p <- pwt$p.value[1]
      estt <- pwt$estimate[1]
      pw_time_dir <- if (pw_time_p < 0.05) {
        if (estt < 0) "D1 < D5" else "D1 > D5"
      } else "D1 ~= D5"
    }

    results[[paste(bm, model_name)]] <- data.frame(
      Biomarker = bm, Model = model_name,
      Mortality_p = round(mort_p, 3), Time_p = round(time_p, 3), Interaction_p = round(inter_p, 3),
      Pairwise_Mortality = sprintf("%s (p=%s)", pw_mort_dir, fmt_p(pw_mort_p)),
      Pairwise_Time = sprintf("%s (p=%s)", pw_time_dir, fmt_p(pw_time_p)),
      N_observations = n_obs, N_patients = n_patients, Convergence = convergence_flag
    )
    cat(sprintf("%s / %s: mort_p=%.3f time_p=%.3f inter_p=%.3f N_obs=%d N_pat=%d %s\n",
                bm, model_name, mort_p, time_p, inter_p, n_obs, n_patients, convergence_flag))
  }
}

tab23 <- do.call(rbind, results)
rownames(tab23) <- NULL
write.csv(tab23, file.path(out_dir, "TABELA23_mortality_LMM_full.csv"), row.names = FALSE)
cat("\n=== Supplementary Table 14 - full reproduction ===\n")
print(tab23, row.names = FALSE)
