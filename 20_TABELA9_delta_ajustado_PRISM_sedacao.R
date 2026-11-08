# 20_TABELA9_delta_ajustado_PRISM_sedacao.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

biomarkers <- list(
  NfL      = list(d1 = "NFLD1,2",          d5 = "NFLD5,2"),
  S100B    = list(d1 = "SB100D1,2",        d5 = "SB100D5,2"),
  `IL-6`   = list(d1 = "IL-6 (pg/mL)D1",   d5 = "il6 deliriumd5"),
  `TNF-alfa` = list(d1 = "TNF (pg/mL)D1",  d5 = "tnfdeliriumd5")
)

results <- list()
for (bm in names(biomarkers)) {
  d1col <- biomarkers[[bm]]$d1; d5col <- biomarkers[[bm]]$d5
  delta <- log(raw[[d5col]]) - log(raw[[d1col]])
  df <- data.frame(id = raw$ID, phenotype = raw$fenotipo_D1, delta = delta,
                    PRISMIV = raw[["PRISMIV%"]], sedation = raw[["sedation time"]])
  df <- df[!is.na(df$delta) & !is.na(df$phenotype), ]
  df$phenotype <- factor(df$phenotype, levels = c("G1_free","G2_delirium","G3_coma"))

  fit <- lm(delta ~ phenotype + PRISMIV + sedation, data = df)
  a <- anova(fit)
  sm <- summary(fit)
  model_F <- sm$fstatistic[1]; model_df1 <- sm$fstatistic[2]; model_df2 <- sm$fstatistic[3]
  model_p <- pf(model_F, model_df1, model_df2, lower.tail = FALSE)

  phen_F <- a["phenotype","F value"]; phen_df <- a["phenotype","Df"]
  phen_p <- a["phenotype","Pr(>F)"]
  resid_df <- a["Residuals","Df"]

  ss_phen <- a["phenotype","Sum Sq"]; ss_resid <- a["Residuals","Sum Sq"]
  partial_eta2 <- ss_phen / (ss_phen + ss_resid)

  coefs <- sm$coefficients
  prism_beta <- coefs["PRISMIV","Estimate"]; prism_se <- coefs["PRISMIV","Std. Error"]; prism_p <- coefs["PRISMIV","Pr(>|t|)"]
  sed_beta <- coefs["sedation","Estimate"]; sed_se <- coefs["sedation","Std. Error"]; sed_p <- coefs["sedation","Pr(>|t|)"]

  results[[bm]] <- data.frame(
    Biomarker = bm,
    Model_F = round(model_F,2), Model_df1 = model_df1, Model_df2 = model_df2, Model_p = round(model_p,3),
    Phenotype_F = round(phen_F,2), Phenotype_df = phen_df, Phenotype_p = round(phen_p,3),
    PRISM_beta = round(prism_beta,3), PRISM_SE = round(prism_se,3), PRISM_p = round(prism_p,3),
    Sedation_beta = round(sed_beta,3), Sedation_SE = round(sed_se,3), Sedation_p = round(sed_p,3),
    Partial_eta2 = round(partial_eta2,3)
  )
}

tab9 <- do.call(rbind, results)
write.csv(tab9, file.path(out_dir, "TABELA9_delta_ajustado_PRISM_sedacao_CONTEMPORANEO.csv"), row.names = FALSE)
cat("=== TABELA 9 SDC - delta ajustado PRISM+sedacao, fenotipo contemporaneo ===\n\n")
print(tab9, row.names = FALSE)
