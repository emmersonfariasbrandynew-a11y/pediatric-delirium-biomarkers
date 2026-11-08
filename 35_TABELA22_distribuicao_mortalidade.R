# 35_TABELA22_distribuicao_mortalidade.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw <- raw[!is.na(raw$ID) & raw$ID != "", ]
raw$obitooverall <- suppressWarnings(as.numeric(raw$obitooverall))
cat(sprintf("Loaded analysis-ready dataset: %d patients\n", nrow(raw)))

biomarkers <- list(
  "CRP (mg/L)"        = list(d1 = "PCRD1",             d5 = "PCRD5",             dp = 1),
  "ESR (mm/h)"        = list(d1 = "VHSD1",             d5 = "VHSD5",             dp = 0),
  "Ferritin (ng/mL)"  = list(d1 = "FERD1",             d5 = "FERD5",             dp = 1),
  "NLR"               = list(d1 = "NLRD1",             d5 = "NLRD5",             dp = 2),
  "PLR"               = list(d1 = "PLRD1",             d5 = "PLRD5",             dp = 2),
  "CRP/Albumin"       = list(d1 = "CRP/ABD1",          d5 = "CRP/ABD5",          dp = 2),
  "VIS (score)"       = list(d1 = "VISd1",             d5 = "VISD5",             dp = 0),
  "IL-2 (pg/mL)"      = list(d1 = "IL-2 (pg/mL)D1",    d5 = "il-2 delirium d5",  dp = 2),
  "IL-4 (pg/mL)"      = list(d1 = "IL-4 (pg/mL)D1",    d5 = "IL-4 (pg/mL)D5",    dp = 2),
  "IL-6 (pg/mL)"      = list(d1 = "IL-6 (pg/mL)D1",    d5 = "il6 deliriumd5",    dp = 2),
  "IL-10 (pg/mL)"     = list(d1 = "IL-10 (pg/mL)D1",   d5 = "IL-10 (pg/mL)D5",   dp = 2),
  "IL-17A (pg/mL)"    = list(d1 = "IL-17A (pg/mL)D1",  d5 = "il-17Ad5",          dp = 2),
  "IFN-gamma (pg/mL)" = list(d1 = "IFN-γ (pg/mL)D1",   d5 = "ifn-gamad5",        dp = 2),
  "TNF-alfa (pg/mL)"  = list(d1 = "TNF (pg/mL)D1",     d5 = "tnfdeliriumd5",     dp = 2),
  "TEAC (mmol/L)"     = list(d1 = "TEAC D1",           d5 = "TEACD5",            dp = 2),
  "TBARS (umol/L)"    = list(d1 = "TBARSD1",           d5 = "TBARSD5",           dp = 2),
  "GSH (umol/L)"      = list(d1 = "GSHD1",             d5 = "GSHD5",             dp = 3),
  "S100B (ng/mL)"     = list(d1 = "SB100D1,2",         d5 = "SB100D5,2",         dp = 3),
  "NfL (ng/mL)"       = list(d1 = "NFLD1,2",           d5 = "NFLD5,2",           dp = 3)
)

fmt_med_iqr <- function(x, dp) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE, type = 6)
  fmt <- sprintf("%%.%df (%%.%df%s%%.%df)", dp, dp, "–", dp)
  sprintf(fmt, q[1], q[2], q[3])
}
fmt_mean_sd <- function(x) sprintf("%.2f ± %.2f", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))

results <- data.frame(Biomarker=character(), Time=character(), Scale=character(),
                       Survivors=character(), NonSurvivors=character())
for (bm in names(biomarkers)) {
  dp <- biomarkers[[bm]]$dp
  for (tp in c("d1","d5")) {
    col <- biomarkers[[bm]][[tp]]
    v <- suppressWarnings(as.numeric(raw[[col]]))
    surv <- v[raw$obitooverall == 0 & !is.na(v)]
    nsurv <- v[raw$obitooverall == 1 & !is.na(v)]
    results <- rbind(results, data.frame(
      Biomarker = bm, Time = toupper(tp), Scale = "Raw",
      Survivors = fmt_med_iqr(surv, dp), NonSurvivors = fmt_med_iqr(nsurv, dp)
    ))
    results <- rbind(results, data.frame(
      Biomarker = bm, Time = toupper(tp), Scale = "Log",
      Survivors = fmt_mean_sd(log(surv)), NonSurvivors = fmt_mean_sd(log(nsurv))
    ))
  }
  cat(sprintf("%s done\n", bm))
}

write.csv(results, file.path(out_dir, "TABELA22_distribuicao_mortalidade.csv"), row.names = FALSE)
cat("\n=== Supplementary Table 13 - full reproduction ===\n")
print(results, row.names = FALSE)
