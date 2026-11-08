# 07_phase4_mortality.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(survival)
})

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw$PRISMIV <- raw[["PRISMIV%"]]

n_death <- sum(raw$obitooverall == 1)
n_alive <- sum(raw$obitooverall == 0)
cat(sprintf("Obitos: %d/%d | Sobreviventes (censurados): %d/%d\n", n_death, nrow(raw), n_alive, nrow(raw)))
cat(sprintf("Tempo de seguimento (time90): mediana=%.0f dias, IIQ=%.0f-%.0f, min=%.0f, max=%.0f\n",
            median(raw$time90), quantile(raw$time90,.25), quantile(raw$time90,.75), min(raw$time90), max(raw$time90)))
cat("Definicao de censura: sobreviventes censurados em time90 (dias desde a\n")
cat("inclusao/D0 ate obito ou fim do seguimento). Momento inicial: inclusao\n")
cat("no estudo (D0/admissao), conforme reconstrucao da Fase 1.\n\n")


raw$TNFa_D1 <- raw[["TNF (pg/mL)D1"]]
q1_tnfa <- quantile(raw$TNFa_D1, 0.25, type = 6, na.rm = TRUE)
q3_tnfa <- quantile(raw$TNFa_D1, 0.75, type = 6, na.rm = TRUE)
iqr_tnfa <- as.numeric(q3_tnfa - q1_tnfa)
cat(sprintf("IIQ do TNF-alfa D1 (type=6, estilo SPSS) = %.4f pg/mL (Q1=%.3f, Q3=%.3f)\n\n",
            iqr_tnfa, q1_tnfa, q3_tnfa))

cox_tnfa_crude <- coxph(Surv(time90, obitooverall) ~ TNFa_D1, data = raw)
b1  <- coef(cox_tnfa_crude)[1]; se1 <- sqrt(vcov(cox_tnfa_crude)[1,1])
cat(sprintf("Modelo BRUTO: beta=%.8f  EP=%.8f  p=%.6f\n", b1, se1, summary(cox_tnfa_crude)$coefficients[,"Pr(>|z|)"]))
cat(sprintf("HR por 1 pg/mL = %.6f (95%% CI %.6f-%.6f)\n",
            exp(b1), exp(b1 - 1.96*se1), exp(b1 + 1.96*se1)))
hr_iqr <- exp(b1 * iqr_tnfa)
ci_iqr_lo <- exp((b1 - 1.96*se1) * iqr_tnfa)
ci_iqr_hi <- exp((b1 + 1.96*se1) * iqr_tnfa)
cat(sprintf("HR por incremento de 1 IIQ = %.4f (95%% CI %.4f-%.4f)\n\n", hr_iqr, ci_iqr_lo, ci_iqr_hi))

cox_tnfa_adj <- coxph(Surv(time90, obitooverall) ~ TNFa_D1 + PRISMIV, data = raw)
b2  <- coef(cox_tnfa_adj)[1]; se2 <- sqrt(vcov(cox_tnfa_adj)[1,1])
cat(sprintf("Modelo + PRISM IV: beta(TNFa)=%.8f  EP=%.8f  p=%.6f\n",
            b2, se2, summary(cox_tnfa_adj)$coefficients[1,"Pr(>|z|)"]))
cat(sprintf("HR por 1 pg/mL (ajustado) = %.6f (95%% CI %.6f-%.6f)\n",
            exp(b2), exp(b2 - 1.96*se2), exp(b2 + 1.96*se2)))
hr_iqr_adj <- exp(b2 * iqr_tnfa)
cat(sprintf("HR por incremento de 1 IIQ (ajustado) = %.4f (95%% CI %.4f-%.4f)\n\n",
            hr_iqr_adj, exp((b2-1.96*se2)*iqr_tnfa), exp((b2+1.96*se2)*iqr_tnfa)))


for (mname in c("cox_tnfa_crude","cox_tnfa_adj")) {
  fit <- get(mname)
  zp <- cox.zph(fit)
  cat(sprintf("-- %s --\n", mname)); print(zp$table)
  png(file.path(out_dir, sprintf("schoenfeld_%s.png", mname)), width = 800, height = 400 * ceiling(nrow(zp$table)-1)/1)
  plot(zp)
  dev.off()
  cat("\n")
}


set.seed(20260818)
R <- 2000
n <- nrow(raw)
boot_betas <- numeric(R); boot_cindex <- numeric(R)
apparent_cindex <- summary(cox_tnfa_crude)$concordance[1]

for (i in seq_len(R)) {
  idx <- sample(seq_len(n), n, replace = TRUE)
  bs <- raw[idx, ]
  fit_bs <- tryCatch(coxph(Surv(time90, obitooverall) ~ TNFa_D1, data = bs), error = function(e) NULL)
  if (is.null(fit_bs)) { boot_betas[i] <- NA; boot_cindex[i] <- NA; next }
  boot_betas[i] <- coef(fit_bs)[1]
  c_orig <- tryCatch(summary(coxph(Surv(time90, obitooverall) ~ TNFa_D1, data = raw,
                                    init = coef(fit_bs), control = coxph.control(iter.max = 0)))$concordance[1],
                      error = function(e) NA)
  boot_cindex[i] <- summary(fit_bs)$concordance[1] - c_orig  # optimism for this resample
}
boot_betas <- boot_betas[!is.na(boot_betas)]
boot_cindex <- boot_cindex[!is.na(boot_cindex)]

cat(sprintf("Reamostragens validas: %d de %d\n", length(boot_betas), R))
cat(sprintf("Beta original=%.6f | Beta bootstrap medio=%.6f | IC95%% percentil beta: %.6f a %.6f\n",
            b1, mean(boot_betas), quantile(boot_betas,.025), quantile(boot_betas,.975)))
cat(sprintf("HR original=%.4f | HR bootstrap (percentil 2.5-97.5): %.4f - %.4f\n",
            exp(b1), exp(quantile(boot_betas,.025)), exp(quantile(boot_betas,.975))))
cat(sprintf("C-index aparente=%.4f | otimismo medio do C-index=%.4f | C-index corrigido=%.4f\n\n",
            apparent_cindex, mean(boot_cindex), apparent_cindex - mean(boot_cindex)))

png(file.path(out_dir, "bootstrap_TNFa_HR_distribution.png"), width = 700, height = 500)
hist(exp(boot_betas), breaks = 40, main = "Distribuicao bootstrap do HR (TNF-alfa D1)",
     xlab = "HR por 1 pg/mL", col = "lightblue")
abline(v = exp(b1), col = "red", lwd = 2, lty = 2)
dev.off()

saveRDS(list(cox_tnfa_crude = cox_tnfa_crude,
             cox_tnfa_adj = cox_tnfa_adj, boot_betas = boot_betas),
        file.path(out_dir, "phase4_mortality_models.rds"))

cat("\n=== FASE 4 (mortalidade) CONCLUIDA ===\n")
cat("Saidas: schoenfeld_*.png (2 modelos TNF, horizonte aberto),\n")
cat("bootstrap_TNFa_HR_distribution.png, phase4_mortality_models.rds\n")
