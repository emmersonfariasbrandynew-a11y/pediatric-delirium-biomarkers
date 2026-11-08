# 08_pacote2_4_sedacao_completa.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

q6 <- function(x, p) quantile(x, p, type = 6, na.rm = TRUE)  # SPSS-style, para consistencia com o manuscrito (SPSS v30)

drugs <- list(
  list(label = "Midazolam (benzodiazepinico)", use_col = "benzo", dose_col = "dbenzo EV"),
  list(label = "Morfina",                       use_col = "morphine use", dose_col = "morphinedose EV"),
  list(label = "Fentanil",                      use_col = "fentanyl use", dose_col = "fentanyl EV dose"),
  list(label = "Cetamina",                      use_col = "cetamina", dose_col = "dceta EV"),
  list(label = "Dexmedetomidina",                use_col = "precedx", dose_col = "dprecedex"),
  list(label = "Propofol",                      use_col = "Propofol", dose_col = "Propofol dose")
)

cat(sprintf("N total da coorte = %d\n\n", nrow(raw)))

rows <- list()
for (d in drugs) {
  use <- raw[[d$use_col]]
  dose <- raw[[d$dose_col]]
  n_exp <- sum(use == 1, na.rm = TRUE)
  n_miss_use <- sum(is.na(use))
  dose_exposed <- dose[use == 1 & !is.na(use)]
  n_dose_valid <- sum(!is.na(dose_exposed))
  med_dose <- if (n_dose_valid > 0) median(dose_exposed, na.rm = TRUE) else NA
  q1 <- if (n_dose_valid > 0) q6(dose_exposed, .25) else NA
  q3 <- if (n_dose_valid > 0) q6(dose_exposed, .75) else NA

  cat(sprintf("-- %s --\n", d$label))
  cat(sprintf("Exposicao: %d/%d (%.1f%%) | dados de uso ausentes: %d\n",
              n_exp, nrow(raw), 100*n_exp/nrow(raw), n_miss_use))
  if (n_dose_valid > 0) {
    cat(sprintf("Dose entre expostos (n=%d): mediana=%.4f, IIQ=%.4f-%.4f\n",
                n_dose_valid, med_dose, q1, q3))
  } else {
    cat("Dose: sem dados validos entre expostos\n")
  }
  cat("\n")
  rows[[d$label]] <- data.frame(agente = d$label, n_exposto = n_exp, N_total = nrow(raw),
                                 pct_exposto = round(100*n_exp/nrow(raw),1),
                                 n_com_dose = n_dose_valid, dose_mediana = round(med_dose,4),
                                 dose_IIQ_lo = round(q1,4), dose_IIQ_hi = round(q3,4))
}
sed_df <- do.call(rbind, rows)
write.csv(sed_df, file.path(out_dir, "table_2_4_sedacao_completa.csv"), row.names = FALSE)
print(sed_df, row.names = FALSE)


