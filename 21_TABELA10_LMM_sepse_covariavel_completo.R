# 21_TABELA10_LMM_sepse_covariavel_completo.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4); library(lmerTest)
})
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

biomarkers <- list(
  NfL      = list(d1 = "NFLD1,2",          d5 = "NFLD5,2"),
  S100B    = list(d1 = "SB100D1,2",        d5 = "SB100D5,2"),
  `IL-6`   = list(d1 = "IL-6 (pg/mL)D1",   d5 = "il6 deliriumd5"),
  `TNF-alfa` = list(d1 = "TNF (pg/mL)D1",  d5 = "tnfdeliriumd5")
)

build_long <- function(raw, d1col, d5col) {
  data.frame(id = rep(raw$ID, 2),
             time = factor(rep(c("D1","D5"), each = nrow(raw)), levels = c("D1","D5")),
             phenotype = factor(c(raw$fenotipo_D1, raw$fenotipo_D5), levels = c("G1_free","G2_delirium","G3_coma")),
             sepsis = rep(raw$sepsis, 2),
             PRISMIV = rep(raw[["PRISMIV%"]], 2), sedation = rep(raw[["sedation time"]], 2),
             y = c(raw[[d1col]], raw[[d5col]]))
}

results <- list()
for (bm in names(biomarkers)) {
  long <- build_long(raw, biomarkers[[bm]]$d1, biomarkers[[bm]]$d5)
  long <- long[!is.na(long$y) & !is.na(long$phenotype), ]
  long$logy <- log(long$y)

  fit <- lmer(logy ~ phenotype * time + sepsis + PRISMIV + sedation + (1|id), data = long,
              REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  an <- anova(fit, type = 3, ddf = "Satterthwaite")
  cf <- summary(fit)$coefficients

  results[[bm]] <- data.frame(
    Biomarker = bm,
    Phenotype_F = round(an["phenotype","F value"],2), Phenotype_p = round(an["phenotype","Pr(>F)"],3),
    Time_F = round(an["time","F value"],2), Time_p = round(an["time","Pr(>F)"],3),
    Interaction_F = round(an["phenotype:time","F value"],2), Interaction_p = round(an["phenotype:time","Pr(>F)"],3),
    Sepsis_beta = round(cf["sepsis","Estimate"],3), Sepsis_SE = round(cf["sepsis","Std. Error"],3), Sepsis_p = round(cf["sepsis","Pr(>|t|)"],3),
    PRISM_beta = round(cf["PRISMIV","Estimate"],3), PRISM_SE = round(cf["PRISMIV","Std. Error"],3), PRISM_p = round(cf["PRISMIV","Pr(>|t|)"],3),
    Sedation_beta = round(cf["sedation","Estimate"],3), Sedation_SE = round(cf["sedation","Std. Error"],3), Sedation_p = round(cf["sedation","Pr(>|t|)"],3),
    singular = isSingular(fit, tol = 1e-4)
  )
}

tab10 <- do.call(rbind, results)
write.csv(tab10, file.path(out_dir, "TABELA10_LMM_sepse_covariavel_CONTEMPORANEO.csv"), row.names = FALSE)
cat("=== TABELA 10 SDC - LMM +sepse+PRISM+sedacao, fenotipo contemporaneo, 4 biomarcadores ===\n\n")
print(tab10, row.names = FALSE)
