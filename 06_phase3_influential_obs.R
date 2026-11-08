# 06_phase3_influential_obs.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(influence.ME)
})

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
phase2 <- readRDS(file.path(out_dir, "phase2_models.rds"))
long <- phase2$long  # already has logNfL, logS100B, phenotype, time, id

long$TNFa    <- c(raw[["TNF (pg/mL)D1"]], raw[["tnfdeliriumd5"]])
long$logTNFa <- log(long$TNFa)

fit_nfl   <- lmer(logNfL   ~ phenotype * time + (1|id), data = long, REML = TRUE,
                   control = lmerControl(optimizer = "bobyqa"))
fit_s100b <- lmer(logS100B ~ phenotype * time + (1|id), data = long, REML = TRUE,
                   control = lmerControl(optimizer = "bobyqa"))
fit_tnfa  <- lmer(logTNFa  ~ phenotype * time + (1|id), data = long, REML = TRUE,
                   control = lmerControl(optimizer = "bobyqa"))

models <- list(NfL = fit_nfl, S100B = fit_s100b, `TNF-alfa` = fit_tnfa)
n_groups <- length(unique(long$id))
cutoff <- 4 / n_groups

summary_rows <- list()

for (bm in names(models)) {
  fit <- models[[bm]]
  cat(sprintf("---- %s ----\n", bm))

  std_resid <- resid(fit, type = "pearson", scaled = TRUE)
  n_extreme_resid <- sum(abs(std_resid) > 2.5)
  cat(sprintf("Residuos padronizados: %d de %d observacoes com |resid|>2.5\n",
              n_extreme_resid, length(std_resid)))

  infl <- tryCatch(influence(fit, group = "id"), error = function(e) e)
  if (inherits(infl, "error")) {
    cat(sprintf("  ERRO ao calcular influencia: %s\n", conditionMessage(infl)))
    next
  }
  cooksd <- cooks.distance(infl)
  cooksd_df <- data.frame(id = rownames(as.matrix(cooksd)), cooksD = as.numeric(cooksd))
  cooksd_df <- cooksd_df[order(-cooksd_df$cooksD), ]
  top5 <- head(cooksd_df, 5)
  influential_ids <- cooksd_df$id[cooksd_df$cooksD > cutoff]


  png(file.path(out_dir, sprintf("influence_%s.png", gsub("-","",bm))), width = 1000, height = 450)
  par(mfrow = c(1,2))
  plot(fitted(fit), std_resid, main = paste(bm, "- Residuos padronizados"),
       xlab = "Ajustado", ylab = "Residuo padronizado", pch = 19, col = rgb(0,0,0,.5))
  abline(h = c(-2.5,0,2.5), lty = c(2,1,2), col = c("red","black","red"))
  plot(seq_along(cooksd_df$cooksD), sort(cooksd_df$cooksD, decreasing = TRUE),
       type = "h", main = paste(bm, "- Distancia de Cook (por paciente)"),
       xlab = "Paciente (ordenado)", ylab = "Cook's D")
  abline(h = cutoff, col = "red", lty = 2)
  dev.off()

  if (length(influential_ids) > 0) {
    long_excl <- long[!(long$id %in% as.numeric(influential_ids)), ]
    yvar <- switch(bm, NfL = "logNfL", S100B = "logS100B", `TNF-alfa` = "logTNFa")
    f2 <- as.formula(sprintf("%s ~ phenotype * time + (1|id)", yvar))
    fit_excl <- tryCatch(lmer(f2, data = long_excl, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa")),
                          error = function(e) e)
    if (!inherits(fit_excl, "error")) {
      cmp <- data.frame(
        coef = names(fixef(fit)),
        com_todos = round(fixef(fit), 4),
        sem_influentes = round(fixef(fit_excl)[names(fixef(fit))], 4)
      )
      cmp$diferenca_abs <- round(abs(cmp$com_todos - cmp$sem_influentes), 4)
      print(cmp, row.names = FALSE)
      an_full  <- anova(fit)
      an_excl  <- anova(fit_excl)
      cat(sprintf("p(fenotipo) com todos = %.4f | sem influentes = %.4f\n",
                   an_full["phenotype","Pr(>F)"], an_excl["phenotype","Pr(>F)"]))
    }
  } else {
    cat("\nNenhum paciente excede o corte - nao ha necessidade de reajuste sem influentes.\n")
  }
  cat("\n")

  summary_rows[[bm]] <- data.frame(
    biomarker = bm, n_extreme_resid = n_extreme_resid,
    n_influential_patients = length(influential_ids),
    influential_ids = ifelse(length(influential_ids) == 0, "", "withheld"),
    max_cooksD = round(max(cooksd_df$cooksD), 4)
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "table23_influential_observations_summary.csv"), row.names = FALSE)

