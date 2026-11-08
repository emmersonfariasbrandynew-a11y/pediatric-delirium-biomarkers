# 34_TABELA2_wilcoxon_pareado_completo.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw <- raw[!is.na(raw$ID) & raw$ID != "", ]
cat(sprintf("Loaded analysis-ready dataset: %d patients\n", nrow(raw)))

biomarkers <- list(
  "CRP, mg/L"        = list(d1 = "PCRD1",             d5 = "PCRD5"),
  "ESR, mm/h"        = list(d1 = "VHSD1",             d5 = "VHSD5"),
  "Ferritin, ng/mL"  = list(d1 = "FERD1",             d5 = "FERD5"),
  "NLR"              = list(d1 = "NLRD1",             d5 = "NLRD5"),
  "PLR"              = list(d1 = "PLRD1",             d5 = "PLRD5"),
  "CRP/albumin ratio"= list(d1 = "CRP/ABD1",          d5 = "CRP/ABD5"),
  "IL-2, pg/mL"      = list(d1 = "IL-2 (pg/mL)D1",    d5 = "il-2 delirium d5"),
  "IL-4, pg/mL"      = list(d1 = "IL-4 (pg/mL)D1",    d5 = "IL-4 (pg/mL)D5"),
  "IL-6, pg/mL"      = list(d1 = "IL-6 (pg/mL)D1",    d5 = "il6 deliriumd5"),
  "IL-10, pg/mL"     = list(d1 = "IL-10 (pg/mL)D1",   d5 = "IL-10 (pg/mL)D5"),
  "IL-17A, pg/mL"    = list(d1 = "IL-17A (pg/mL)D1",  d5 = "il-17Ad5"),
  "IFN-gamma, pg/mL" = list(d1 = "IFN-γ (pg/mL)D1",   d5 = "ifn-gamad5"),
  "TNF-alfa, pg/mL"  = list(d1 = "TNF (pg/mL)D1",     d5 = "tnfdeliriumd5"),
  "TEAC, mmol/L"     = list(d1 = "TEAC D1",           d5 = "TEACD5"),
  "TBARS, umol/L"    = list(d1 = "TBARSD1",           d5 = "TBARSD5"),
  "GSH, umol/mL"     = list(d1 = "GSHD1",             d5 = "GSHD5"),
  "S100B, ng/mL"     = list(d1 = "SB100D1,2",         d5 = "SB100D5,2"),
  "NfL, ng/mL"       = list(d1 = "NFLD1,2",           d5 = "NFLD5,2")
)

fmt_med_iqr <- function(x, dp) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE, type = 6)
  fmt <- sprintf("%%.%df (%%.%df%s%%.%df)", dp, dp, "–", dp)
  sprintf(fmt, q[1], q[2], q[3])
}
fmt_p <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)
decimals <- c("CRP, mg/L"=1, "ESR, mm/h"=1, "Ferritin, ng/mL"=1, "NLR"=1, "PLR"=1,
              "CRP/albumin ratio"=2, "IL-2, pg/mL"=2, "IL-4, pg/mL"=2, "IL-6, pg/mL"=2,
              "IL-10, pg/mL"=2, "IL-17A, pg/mL"=2, "IFN-gamma, pg/mL"=2, "TNF-alfa, pg/mL"=2,
              "TEAC, mmol/L"=2, "TBARS, umol/L"=2, "GSH, umol/mL"=3, "S100B, ng/mL"=3, "NfL, ng/mL"=3)

results <- data.frame(Biomarker = character(), N = integer(), D1 = character(),
                       D5 = character(), p_value = character())
for (bm in names(biomarkers)) {
  d1col <- biomarkers[[bm]]$d1; d5col <- biomarkers[[bm]]$d5
  d1v <- suppressWarnings(as.numeric(raw[[d1col]]))
  d5v <- suppressWarnings(as.numeric(raw[[d5col]]))
  keep <- !is.na(d1v) & !is.na(d5v)
  d1v <- d1v[keep]; d5v <- d5v[keep]
  test <- wilcox.test(d1v, d5v, paired = TRUE)
  dp <- decimals[[bm]]
  results <- rbind(results, data.frame(
    Biomarker = bm, N = length(d1v),
    D1 = fmt_med_iqr(d1v, dp), D5 = fmt_med_iqr(d5v, dp), p_value = fmt_p(test$p.value)
  ))
  cat(sprintf("%-20s n=%d  D1=%s  D5=%s  p=%s\n", bm, length(d1v), fmt_med_iqr(d1v, dp), fmt_med_iqr(d5v, dp), fmt_p(test$p.value)))
}

write.csv(results, file.path(out_dir, "TABELA2_wilcoxon_pareado_18biomarcadores.csv"), row.names = FALSE)
cat("\n=== Supplementary Table 2 - full reproduction (18 biomarkers) ===\n")
print(results, row.names = FALSE)
