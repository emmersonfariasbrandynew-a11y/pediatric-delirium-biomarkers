# 14_pacote_wilcoxon_pareado_D1D5.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

bms <- list(
  "TNF-alfa" = list(d1 = "TNF (pg/mL)D1", d5 = "tnfdeliriumd5"),
  "IL-6"     = list(d1 = "IL-6 (pg/mL)D1", d5 = "il6 deliriumd5"),
  "S100B"    = list(d1 = "SB100D1,2", d5 = "SB100D5,2"),
  "NfL"      = list(d1 = "NFLD1,2", d5 = "NFLD5,2")
)

cat("=== WILCOXON PAREADO (SIGNED-RANK) D1 vs D5 ===\n\n")
results <- list()
for (nm in names(bms)) {
  b <- bms[[nm]]
  d1v <- raw[[b$d1]]; d5v <- raw[[b$d5]]
  ok <- !is.na(d1v) & !is.na(d5v)
  n_ok <- sum(ok)
  w <- wilcox.test(d1v[ok], d5v[ok], paired = TRUE)
  direcao <- ifelse(median(d5v[ok]-d1v[ok]) > 0, "D5 > D1", "D5 < D1")
  cat(sprintf("%-10s n=%d | mediana D1=%.3f D5=%.3f | direcao=%s | Wilcoxon V=%.0f p=%.6f\n",
              nm, n_ok, median(d1v[ok]), median(d5v[ok]), direcao, w$statistic, w$p.value))
  results[[nm]] <- data.frame(biomarker = nm, n = n_ok, mediana_D1 = median(d1v[ok]), mediana_D5 = median(d5v[ok]),
                               direcao = direcao, V = w$statistic, p_valor = w$p.value)
}
res_df <- do.call(rbind, results)
rownames(res_df) <- NULL
write.csv(res_df, file.path(out_dir, "CAMADA1_wilcoxon_pareado_D1D5.csv"), row.names = FALSE)
print(res_df, row.names = FALSE)

cat("\n=== CONCLUIDO ===\n")
