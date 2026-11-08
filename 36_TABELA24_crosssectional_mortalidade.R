# 36_TABELA24_crosssectional_mortalidade.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw <- raw[!is.na(raw$ID) & raw$ID != "", ]
raw$obitooverall <- suppressWarnings(as.numeric(raw$obitooverall))
cat(sprintf("Loaded analysis-ready dataset: %d patients\n", nrow(raw)))

biomarkers <- list(
  "CRP (mg/L)"        = "PCRD1|PCRD5",
  "ESR (mm/h)"        = "VHSD1|VHSD5",
  "Ferritin (ng/mL)"  = "FERD1|FERD5",
  "NLR"               = "NLRD1|NLRD5",
  "PLR"               = "PLRD1|PLRD5",
  "CRP/Albumin"       = "CRP/ABD1|CRP/ABD5",
  "VIS (score)"       = "VISd1|VISD5",
  "IL-2 (pg/mL)"      = "IL-2 (pg/mL)D1|il-2 delirium d5",
  "IL-4 (pg/mL)"      = "IL-4 (pg/mL)D1|IL-4 (pg/mL)D5",
  "IL-6 (pg/mL)"      = "IL-6 (pg/mL)D1|il6 deliriumd5",
  "IL-10 (pg/mL)"     = "IL-10 (pg/mL)D1|IL-10 (pg/mL)D5",
  "IL-17A (pg/mL)"    = "IL-17A (pg/mL)D1|il-17Ad5",
  "IFN-gamma (pg/mL)" = "IFN-γ (pg/mL)D1|ifn-gamad5",
  "TNF-alfa (pg/mL)"  = "TNF (pg/mL)D1|tnfdeliriumd5",
  "TEAC (mmol/L)"     = "TEAC D1|TEACD5",
  "TBARS (umol/L)"    = "TBARSD1|TBARSD5",
  "GSH (umol/L)"      = "GSHD1|GSHD5",
  "S100B (ng/mL)"     = "SB100D1,2|SB100D5,2",
  "NfL (ng/mL)"       = "NFLD1,2|NFLD5,2"
)

wilcox_z <- function(nsurv, surv) {
  n1 <- length(nsurv); n2 <- length(surv)
  r <- rank(c(surv, nsurv))  # survivors first as reference group x
  r1 <- r[1:n2]
  W <- sum(r1) - n2*(n2+1)/2
  mu <- n2*n1/2
  ties <- table(c(surv, nsurv))
  tie_correction <- sum(ties^3 - ties) / ((n1+n2)*(n1+n2-1))
  sigma2 <- (n1*n2/12) * ((n1+n2+1) - tie_correction)
  correction <- sign(W - mu) * 0.5
  if (W == mu) correction <- 0
  Z <- (W - mu - correction) / sqrt(sigma2)
  p <- 2 * pnorm(-abs(Z))
  list(Z = Z, p = p)
}

direction_label <- function(Z, p, group="") {
  if (p >= 0.05) return(sprintf("No difference (p = %.3f)", p))
  if (Z < 0) return(sprintf("Higher in non-survivors (p = %.3f)", p))
  return(sprintf("Higher in survivors (p = %.3f)", p))
}

cohend <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  sp <- sqrt(((n1-1)*var(x) + (n2-1)*var(y)) / (n1+n2-2))
  (mean(x) - mean(y)) / sp
}

results <- data.frame(Biomarker=character(), D1_Z=numeric(), D1_p=numeric(), D1_dir=character(),
                       D5_Z=numeric(), D5_p=numeric(), D5_dir=character(),
                       DeltaLog_MeanDiff=numeric(), t=numeric(), df=numeric(), p=numeric(),
                       CI_lo=numeric(), CI_hi=numeric(), Cohens_d=numeric(), Delta_dir=character())

for (bm in names(biomarkers)) {
  cols <- strsplit(biomarkers[[bm]], "\\|")[[1]]
  d1col <- cols[1]; d5col <- cols[2]
  d1v <- suppressWarnings(as.numeric(raw[[d1col]]))
  d5v <- suppressWarnings(as.numeric(raw[[d5col]]))

  keep1 <- !is.na(d1v)
  surv1 <- d1v[keep1 & raw$obitooverall==0]; nsurv1 <- d1v[keep1 & raw$obitooverall==1]
  r1 <- wilcox_z(nsurv1, surv1)

  keep5 <- !is.na(d5v)
  surv5 <- d5v[keep5 & raw$obitooverall==0]; nsurv5 <- d5v[keep5 & raw$obitooverall==1]
  r5 <- wilcox_z(nsurv5, surv5)

  keepD <- !is.na(d1v) & !is.na(d5v)
  dlog <- log(d5v[keepD]) - log(d1v[keepD])
  obito_d <- raw$obitooverall[keepD]
  dsurv <- dlog[obito_d==0]; dnsurv <- dlog[obito_d==1]
  tt <- t.test(dsurv, dnsurv)
  d <- cohend(dsurv, dnsurv)
  dir_delta <- if (tt$p.value >= 0.05) sprintf("No difference (p = %.3f)", tt$p.value) else
    if (mean(dnsurv) > mean(dsurv)) sprintf("Higher in non-survivors (p = %.3f)", tt$p.value) else
      sprintf("Higher in survivors (p = %.3f)", tt$p.value)

  results <- rbind(results, data.frame(
    Biomarker = bm,
    D1_Z = round(r1$Z,3), D1_p = round(r1$p,3), D1_dir = direction_label(r1$Z, r1$p),
    D5_Z = round(r5$Z,3), D5_p = round(r5$p,3), D5_dir = direction_label(r5$Z, r5$p),
    DeltaLog_MeanDiff = round(mean(dsurv)-mean(dnsurv),3),
    t = round(unname(tt$statistic),3), df = round(unname(tt$parameter),1), p = round(tt$p.value,3),
    CI_lo = round(tt$conf.int[1],2), CI_hi = round(tt$conf.int[2],2),
    Cohens_d = round(d,2), Delta_dir = dir_delta
  ))
  cat(sprintf("%s: D1 Z=%.3f p=%.3f | D5 Z=%.3f p=%.3f | dlog t=%.3f p=%.3f d=%.2f\n",
              bm, r1$Z, r1$p, r5$Z, r5$p, unname(tt$statistic), tt$p.value, d))
}

write.csv(results, file.path(out_dir, "TABELA24_crosssectional_mortalidade.csv"), row.names = FALSE)
cat("\n=== Supplementary Table 15 - full reproduction ===\n")
print(results, row.names = FALSE)
