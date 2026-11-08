# 39_TABELA4_distribuicao_fenotipo_contemporaneo.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw <- raw[!is.na(raw$ID) & raw$ID != "", ]

biomarkers <- list(
  "CRP (mg/L)"        = list(d1 = "PCRD1",            d5 = "PCRD5",            dp = 1, lp = 2),
  "ESR (mm/h)"        = list(d1 = "VHSD1",            d5 = "VHSD5",            dp = 0, lp = 2),
  "Ferritin (ng/mL)"  = list(d1 = "FERD1",            d5 = "FERD5",            dp = 1, lp = 2),
  "NLR"               = list(d1 = "NLRD1",            d5 = "NLRD5",            dp = 2, lp = 2),
  "PLR"               = list(d1 = "PLRD1",            d5 = "PLRD5",            dp = 2, lp = 2),
  "CRP/Albumin ratio" = list(d1 = "CRP/ABD1",         d5 = "CRP/ABD5",         dp = 2, lp = 2),
  "VIS"               = list(d1 = "VISd1",            d5 = "VISD5",            dp = 0, lp = 2),
  "IL-2 (pg/mL)"      = list(d1 = "IL-2 (pg/mL)D1",   d5 = "il-2 delirium d5", dp = 2, lp = 2),
  "IL-4 (pg/mL)"      = list(d1 = "IL-4 (pg/mL)D1",   d5 = "IL-4 (pg/mL)D5",   dp = 2, lp = 2),
  "IL-6 (pg/mL)"      = list(d1 = "IL-6 (pg/mL)D1",   d5 = "il6 deliriumd5",   dp = 2, lp = 2),
  "IL-10 (pg/mL)"     = list(d1 = "IL-10 (pg/mL)D1",  d5 = "IL-10 (pg/mL)D5",  dp = 2, lp = 2),
  "IL-17A (pg/mL)"    = list(d1 = "IL-17A (pg/mL)D1", d5 = "il-17Ad5",         dp = 2, lp = 2),
  "IFN-gamma (pg/mL)" = list(d1 = "IFN-γ (pg/mL)D1",  d5 = "ifn-gamad5",       dp = 2, lp = 2),
  "TNF-alfa (pg/mL)"  = list(d1 = "TNF (pg/mL)D1",    d5 = "tnfdeliriumd5",    dp = 2, lp = 2),
  "TEAC (mmol/L)"     = list(d1 = "TEAC D1",          d5 = "TEACD5",           dp = 2, lp = 3),
  "TBARS (umol/L)"    = list(d1 = "TBARSD1",          d5 = "TBARSD5",          dp = 2, lp = 3),
  "GSH (umol/mL)"     = list(d1 = "GSHD1",            d5 = "GSHD5",            dp = 3, lp = 3),
  "S100B (ng/mL)"     = list(d1 = "SB100D1,2",        d5 = "SB100D5,2",        dp = 3, lp = 3),
  "NfL (ng/mL)"       = list(d1 = "NFLD1,2",          d5 = "NFLD5,2",          dp = 3, lp = 3)
)
groups <- c(G1 = "G1_free", G2 = "G2_delirium", G3 = "G3_coma")

fmt_med_iqr <- function(x, dp) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE, type = 6)
  sprintf(sprintf("%%.%df (%%.%df%s%%.%df)", dp, dp, "–", dp), q[1], q[2], q[3])
}
fmt_mean_sd <- function(x, lp) sprintf(sprintf("%%.%df ± %%.%df", lp, lp), mean(x), sd(x))

results <- NULL
counts <- NULL
for (bm in names(biomarkers)) {
  b <- biomarkers[[bm]]
  for (tp in c("d1", "d5")) {
    v <- suppressWarnings(as.numeric(raw[[b[[tp]]]]))
    ph <- if (tp == "d1") raw$fenotipo_D1 else raw$fenotipo_D5
    raw_row <- data.frame(Biomarker = bm, Time = toupper(tp), Scale = "Raw")
    log_row <- data.frame(Biomarker = bm, Time = toupper(tp), Scale = "Log")
    for (g in names(groups)) {
      x <- v[ph == groups[[g]] & !is.na(ph) & !is.na(v)]
      raw_row[[g]] <- fmt_med_iqr(x, b$dp)
      log_row[[g]] <- fmt_mean_sd(log(x), b$lp)
      counts <- rbind(counts, data.frame(Biomarker = bm, Time = toupper(tp), Group = g, n = length(x)))
    }
    results <- rbind(results, raw_row, log_row)
  }
}
write.csv(results, file.path(out_dir, "TABELA4_distribuicao_fenotipo_contemporaneo.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "TABELA4_n_por_grupo.csv"), row.names = FALSE)
cat("Group sizes (contemporaneous phenotype):\n")
print(table(D1 = raw$fenotipo_D1)); print(table(D5 = raw$fenotipo_D5))
cat(sprintf("Biomarker-level n range: %d-%d\n", min(counts$n), max(counts$n)))
print(results, row.names = FALSE)
