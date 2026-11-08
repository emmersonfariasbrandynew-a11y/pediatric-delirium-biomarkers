# 09_pacote_wat1_capd_e_checks.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))


cat("Definicao: CAPD positivo = CAPDTOTAL >=9 em pelo menos um momento (D1 ou D5).\n")
cat("WAT-1 positivo = pelo menos um WAT-1 (D1 ou D5) > 3 (ponto de corte do instrumento).\n\n")

any_positive <- function(a, b) {
  ifelse(is.na(a) & is.na(b), NA,
         ifelse(!is.na(a) & a, TRUE,
                ifelse(!is.na(b) & b, TRUE, FALSE)))
}
capd_pos <- any_positive(raw$CAPDTOTALD1 >= 9, raw$CAPDTOTALD5 >= 9)
wat1_pos <- any_positive(raw[["wat-1 d1"]] > 3, raw[["wat-1 d5"]] > 3)
n_wat1_avaliados <- sum(!is.na(raw[["wat-1 d1"]]) | !is.na(raw[["wat-1 d5"]]))

n_capd_nao_avaliavel <- sum(is.na(capd_pos))
cat(sprintf("Pacientes com CAPD nao avaliavel (RASS<=-4 em D1 E D5, coma nos dois momentos): %d/%d\n",
            n_capd_nao_avaliavel, nrow(raw)))
cat("  -> Estrutural, nao e dado faltante: CAPD nao e aplicavel a paciente comatoso.\n")
cat(sprintf("Pacientes avaliados com WAT-1 (D1 e/ou D5 disponivel): %d/%d\n", n_wat1_avaliados, nrow(raw)))
cat(sprintf("Pacientes com pelo menos 1 WAT-1 acima do corte (>3): %d/%d\n", sum(wat1_pos, na.rm=TRUE), nrow(raw)))
cat(sprintf("Pacientes com CAPD>=9 em pelo menos 1 momento (entre os %d avaliaveis): %d/%d\n\n",
            nrow(raw) - n_capd_nao_avaliavel, sum(capd_pos, na.rm=TRUE), nrow(raw) - n_capd_nao_avaliavel))

capd_pos_labeled <- factor(
  ifelse(is.na(capd_pos), "Nao avaliavel (coma, RASS<=-4)",
         ifelse(capd_pos, "CAPD positivo", "CAPD negativo")),
  levels = c("CAPD positivo", "CAPD negativo", "Nao avaliavel (coma, RASS<=-4)")
)
tab <- table(CAPD = capd_pos_labeled, WAT1_positivo = wat1_pos)
cat("Tabela cruzada CAPD x WAT-1 (N=45, com categoria explicita para nao-avaliavel):\n")
print(addmargins(tab))
n_both <- sum(capd_pos & wat1_pos, na.rm = TRUE)
n_avaliaveis <- nrow(raw) - n_capd_nao_avaliavel
cat(sprintf("\nPositivos simultaneamente para CAPD e WAT-1: %d/%d avaliaveis para CAPD (%.1f%%)\n",
            n_both, n_avaliaveis, 100*n_both/n_avaliaveis))

write.csv(as.data.frame.matrix(addmargins(tab)), file.path(out_dir, "table_2_5_CAPDxWAT1.csv"))


cat("\n=== CONCLUIDO ===\n")
