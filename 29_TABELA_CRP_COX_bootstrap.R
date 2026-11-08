# 29_TABELA_CRP_COX_bootstrap.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages(library(survival))
set.seed(20260818)

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw$PRISMIV <- raw[["PRISMIV%"]]
raw$PCR_D1  <- raw[["PCRD1"]]

d <- raw[!is.na(raw$PCR_D1) & !is.na(raw$time90) & !is.na(raw$obitooverall), ]
n <- nrow(d)
cat(sprintf("N com PCR D1 valido: %d\n", n))
cat(sprintf("Obitos (horizonte aberto): %d/%d (%.1f%%)\n", sum(d$obitooverall==1), n, 100*mean(d$obitooverall==1)))

fit_crude <- coxph(Surv(time90, obitooverall) ~ PCR_D1, data = d)
fit_adj   <- coxph(Surv(time90, obitooverall) ~ PCR_D1 + PRISMIV, data = d)

hr_crude <- exp(coef(fit_crude)["PCR_D1"])
hr_adj   <- exp(coef(fit_adj)["PCR_D1"])
p_crude  <- summary(fit_crude)$coefficients["PCR_D1","Pr(>|z|)"]
p_adj    <- summary(fit_adj)$coefficients["PCR_D1","Pr(>|z|)"]
ci_crude <- exp(confint(fit_crude))["PCR_D1",]
ci_adj   <- exp(confint(fit_adj))["PCR_D1",]

R_BOOT <- 5000
boot_hr_crude <- numeric(R_BOOT)
boot_hr_adj   <- numeric(R_BOOT)
for (i in 1:R_BOOT) {
  idx <- sample(seq_len(n), n, replace = TRUE)
  db <- d[idx, ]
  bc <- tryCatch(coxph(Surv(time90, obitooverall) ~ PCR_D1, data = db), error = function(e) NULL)
  boot_hr_crude[i] <- if (is.null(bc)) NA else exp(coef(bc)["PCR_D1"])
  ba <- tryCatch(coxph(Surv(time90, obitooverall) ~ PCR_D1 + PRISMIV, data = db), error = function(e) NULL)
  boot_hr_adj[i] <- if (is.null(ba)) NA else exp(coef(ba)["PCR_D1"])
}
boot_ci_crude <- quantile(boot_hr_crude, c(0.025, 0.975), na.rm = TRUE)
boot_ci_adj   <- quantile(boot_hr_adj, c(0.025, 0.975), na.rm = TRUE)
boot_p_crude  <- mean(boot_hr_crude <= 1, na.rm = TRUE) * 2
boot_p_crude  <- min(boot_p_crude, 2 - boot_p_crude)
boot_p_adj    <- mean(boot_hr_adj <= 1, na.rm = TRUE) * 2
boot_p_adj    <- min(boot_p_adj, 2 - boot_p_adj)

result <- data.frame(
  Model = c("Crude", "Adjusted (PRISM-IV)"),
  HR = round(c(hr_crude, hr_adj), 4),
  CI95_lo = round(c(ci_crude[1], ci_adj[1]), 4),
  CI95_hi = round(c(ci_crude[2], ci_adj[2]), 4),
  p_value = round(c(p_crude, p_adj), 3),
  Bootstrap_p = round(c(boot_p_crude, boot_p_adj), 3),
  Bootstrap_CI_lo = round(c(boot_ci_crude[1], boot_ci_adj[1]), 4),
  Bootstrap_CI_hi = round(c(boot_ci_crude[2], boot_ci_adj[2]), 4)
)
write.csv(result, file.path(out_dir, "CRP_Cox_openhorizon_bootstrap.csv"), row.names = FALSE)
cat("\n=== CRP (PCR) Day-1 Cox, open horizon, with bootstrap ===\n")
print(result, row.names = FALSE)
