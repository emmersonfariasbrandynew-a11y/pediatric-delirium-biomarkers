# 43_TNF_D1_influence_check.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages(library(survival))

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw$TNF_D1 <- raw[["TNF (pg/mL)D1"]]

d <- raw[!is.na(raw$TNF_D1) & !is.na(raw$time90) & !is.na(raw$obitooverall), ]
cat(sprintf("N patients in crude Day-1 TNF-alfa Cox model: %d\n", nrow(d)))

fit <- coxph(Surv(time90, obitooverall) ~ TNF_D1, data = d)
cat("\n=== Full-sample crude model ===\n")
print(summary(fit)$coefficients)
hr_full <- exp(coef(fit)["TNF_D1"])
p_full  <- summary(fit)$coefficients["TNF_D1", "Pr(>|z|)"]
cat(sprintf("HR=%.4f, p=%.4f\n", hr_full, p_full))

db <- residuals(fit, type = "dfbeta")
d$dfbeta <- as.numeric(db)
d_sorted <- d[order(-abs(d$dfbeta)), c("ID", "TNF_D1", "obitooverall", "time90", "dfbeta")]
highest_dfbeta_id <- d_sorted$ID[1]

d_loo <- d[d$ID != highest_dfbeta_id, ]
fit_loo <- coxph(Surv(time90, obitooverall) ~ TNF_D1, data = d_loo)
hr_loo <- exp(coef(fit_loo)["TNF_D1"])
p_loo  <- summary(fit_loo)$coefficients["TNF_D1", "Pr(>|z|)"]
cat("\n=== Leave-one-out (excluding the patient with the highest |dfbeta|) ===\n")
cat(sprintf("HR=%.4f (full sample: %.4f), p=%.4f (full sample: %.4f)\n", hr_loo, hr_full, p_loo, p_full))

cat("\n=== Full leave-one-out loop (does any single exclusion change p<0.05 status?) ===\n")
flips <- 0
flip_id <- NA
for (id in d$ID) {
  dsub <- d[d$ID != id, ]
  fsub <- tryCatch(coxph(Surv(time90, obitooverall) ~ TNF_D1, data = dsub), error = function(e) NULL)
  if (is.null(fsub)) next
  psub <- summary(fsub)$coefficients["TNF_D1", "Pr(>|z|)"]
  if ((p_full < 0.05) != (psub < 0.05)) {
    flips <- flips + 1
    flip_id <- id
  }
}
cat(sprintf("Total patients whose individual exclusion flips significance status: %d / %d\n", flips, nrow(d)))

result <- data.frame(
  N = nrow(d), HR_full = round(hr_full,3), p_full = round(p_full,4),
  highest_dfbeta_ID = "withheld",
  HR_leave1out_dfbeta = round(hr_loo,3), p_leave1out_dfbeta = round(p_loo,4),
  significance_flip_ID = ifelse(is.na(flip_id), NA, "withheld"),
  n_flips_significance = flips
)
write.csv(result, file.path(out_dir, "TABELA27_TNF_D1_influence_check.csv"), row.names = FALSE)
cat("\nSalvo: output/TABELA27_TNF_D1_influence_check.csv\n")
