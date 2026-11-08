# Reproduces the original ROC selection scores; subgroup components are used only
# to document selection, not as separately reported subgroup inference.
# 33_TABELA14_consistency_score.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({library(lme4); library(lmerTest)})

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw$PRISMIV <- raw[["PRISMIV%"]]

biomarkers <- list(
  NfL      = list(d1 = "NFLD1,2",        d5 = "NFLD5,2"),
  S100B    = list(d1 = "SB100D1,2",      d5 = "SB100D5,2"),
  `IL-6`   = list(d1 = "IL-6 (pg/mL)D1", d5 = "il6 deliriumd5"),
  `TNF-alfa` = list(d1 = "TNF (pg/mL)D1", d5 = "tnfdeliriumd5")
)

long_all <- do.call(rbind, lapply(names(biomarkers), function(bm) {
  d1col <- biomarkers[[bm]]$d1; d5col <- biomarkers[[bm]]$d5
  data.frame(id = raw$ID, biomarker = bm,
             logy = c(log(raw[[d1col]]), log(raw[[d5col]])),
             time = rep(c("D1","D5"), each = nrow(raw)),
             phenotype = c(raw$fenotipo_D1, raw$fenotipo_D5),
             sepsis = rep(raw$sepsis, 2),
             PRISMIV = rep(raw[["PRISMIV%"]], 2),
             sedation = rep(raw[["sedation time"]], 2))
}))

crit4 <- sapply(names(biomarkers), function(bm) {
  d <- long_all[long_all$biomarker == bm & !is.na(long_all$logy) & !is.na(long_all$phenotype), ]
  fit <- lmer(logy ~ phenotype * time + sepsis + PRISMIV + sedation + (1|id), data = d, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))
  a <- anova(fit, type = 3, ddf = "Satterthwaite")
  a["phenotype", "Pr(>F)"]
})

run_delta <- function(bm, septic_only) {
  d1col <- biomarkers[[bm]]$d1; d5col <- biomarkers[[bm]]$d5
  delta <- log(raw[[d5col]]) - log(raw[[d1col]])
  df <- data.frame(id = raw$ID, phenotype = raw$fenotipo_D1, delta = delta,
                    PRISMIV = raw[["PRISMIV%"]], sedation = raw[["sedation time"]],
                    sepsis = raw$sepsis)
  if (septic_only) df <- df[df$sepsis == 1, ]
  df <- df[!is.na(df$delta) & !is.na(df$phenotype), ]
  df$phenotype <- factor(df$phenotype, levels = c("G1_free","G2_delirium","G3_coma"))
  fit <- lm(delta ~ phenotype + PRISMIV + sedation, data = df)
  a <- anova(fit)
  a["phenotype", "Pr(>F)"]
}

run_lmm_septic <- function(bm) {
  d <- long_all[long_all$biomarker == bm & long_all$sepsis == 1 &
                !is.na(long_all$logy) & !is.na(long_all$phenotype), ]
  fit <- lmer(logy ~ phenotype * time + PRISMIV + sedation + (1|id), data = d, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))
  a <- anova(fit, type = 3, ddf = "Satterthwaite")
  a["phenotype", "Pr(>F)"]
}

results <- data.frame(Biomarker = character(), LMM_Septic = numeric(),
                       Delta_Adjusted = numeric(), Delta_Septic = numeric(),
                       LMM_Sepsis_Adjusted_Overall = numeric(),
                       Consistency_Score = character())
for (bm in names(biomarkers)) {
  p1 <- run_lmm_septic(bm)
  p2 <- run_delta(bm, septic_only = FALSE)
  p3 <- run_delta(bm, septic_only = TRUE)
  p4 <- crit4[[bm]]
  score <- sum(c(p1, p2, p3, p4) < 0.05)
  results <- rbind(results, data.frame(
    Biomarker = bm, LMM_Septic = round(p1, 3), Delta_Adjusted = round(p2, 3),
    Delta_Septic = round(p3, 3), LMM_Sepsis_Adjusted_Overall = round(p4, 3),
    Consistency_Score = sprintf("%d / 4", score)
  ))
}
selection <- data.frame(Biomarker=results$Biomarker,
  Original_Consistency_Score=results$Consistency_Score,
  Included_in_Figure3=ifelse(results$Biomarker=="TNF-alfa", "No", "Yes"))
write.csv(selection, file.path(out_dir, "TABELA14_consistency_score.csv"), row.names=FALSE)
cat("Original ROC selection scores (phenotype main effects only):\n")
print(selection, row.names=FALSE)
