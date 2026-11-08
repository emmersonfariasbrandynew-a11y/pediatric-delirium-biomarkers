# 10_pacote_AICBIC_NfL_S100B.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages(library(lme4))
out_dir <- "output"
mods <- readRDS(file.path(out_dir, "phase2_models.rds"))

extract_row <- function(fit, biomarker, spec) {
  if (inherits(fit, "error") || inherits(fit, "simpleError")) {
    return(data.frame(biomarker = biomarker, especificacao = spec,
                       N_obs = NA, N_pacientes = NA, AIC = NA, BIC = NA, logLik = NA,
                       n_parametros = NA, singular = NA, convergiu = FALSE,
                       nota = paste("erro:", conditionMessage(fit))))
  }
  s <- summary(fit)
  ll <- logLik(fit)
  data.frame(
    biomarker      = biomarker,
    especificacao  = spec,
    N_obs          = nobs(fit),
    N_pacientes    = ngrps(fit)[[1]],
    AIC            = round(AIC(fit), 2),
    BIC            = round(BIC(fit), 2),
    logLik         = round(as.numeric(ll), 3),
    n_parametros   = attr(ll, "df"),
    singular       = isSingular(fit, tol = 1e-4),
    convergiu      = TRUE,
    nota           = ""
  )
}

spec_labels <- c(
  M0_raw = "M0: fenotipo*tempo (bruto)",
  M1_PRISM = "M1: + PRISM IV",
  M2_sedation = "M2: + duracao sedacao",
  M3_PRISM_sedation = "M3: + PRISM IV + duracao sedacao"
)

rows <- list()
for (nm in names(spec_labels)) {
  rows[[paste0("NfL_", nm)]]   <- extract_row(mods$nfl_fits[[nm]],   "NfL",   spec_labels[[nm]])
  rows[[paste0("S100B_", nm)]] <- extract_row(mods$s100b_fits[[nm]], "S100B", spec_labels[[nm]])
}
tab <- do.call(rbind, rows)
rownames(tab) <- NULL

cat("=== TABELA AIC/BIC - MODELOS NfL E S100B (M0-M3) ===\n\n")
print(tab, row.names = FALSE)

cat("\n=== MELHOR ESPECIFICACAO POR AIC (menor = melhor ajuste) ===\n")
for (bm in c("NfL","S100B")) {
  sub <- tab[tab$biomarker == bm & tab$convergiu, ]
  best <- sub[which.min(sub$AIC), ]
  cat(sprintf("%s: %s (AIC=%.2f, BIC=%.2f)\n", bm, best$especificacao, best$AIC, best$BIC))
}

write.csv(tab, file.path(out_dir, "table_AICBIC_NfL_S100B.csv"), row.names = FALSE)


cat("\n=== CONCLUIDO ===\n")
