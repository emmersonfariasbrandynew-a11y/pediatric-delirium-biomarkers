# 01_phase1_data_reconstruction.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)

data_path <- "data/delirium_op.csv"
out_dir   <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

raw <- read.csv(data_path, header = TRUE, sep = ",", dec = ".",
                 na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM",
                 check.names = FALSE)

raw <- raw[!is.na(raw$ID) & raw$ID != "", ]
raw$ID <- as.integer(raw$ID)

empty_name_cols <- which(names(raw) == "" | is.na(names(raw)))
if (length(empty_name_cols) > 0) {
  cat(sprintf("Removendo %d coluna(s) sem nome (indices: %s)\n",
              length(empty_name_cols), paste(empty_name_cols, collapse = ",")))
  raw <- raw[, -empty_name_cols]
}
dup_names <- names(raw)[duplicated(names(raw))]
if (length(dup_names) > 0) {
  cat(sprintf("Colunas com nome duplicado (mantendo 1a ocorrencia, removendo repeticoes): %s\n",
              paste(unique(dup_names), collapse = ", ")))
  raw <- raw[, !duplicated(names(raw))]
}
cat(sprintf("Loaded %d patient rows (expected 45)\n\n", nrow(raw)))
stopifnot(nrow(raw) == 45)



classify_phenotype <- function(rass, capd) {
  ifelse(is.na(rass), NA,
    ifelse(rass <= -4, "G3_coma",
      ifelse(is.na(capd), NA,
        ifelse(capd >= 9, "G2_delirium", "G1_free"))))
}

raw$fenotipo_D1 <- classify_phenotype(raw$rassd1, raw$CAPDTOTALD1)
raw$fenotipo_D5 <- classify_phenotype(raw$rassd5, raw$CAPDTOTALD5)

cat("Distribuicao fenotipo_D1 (RASS/CAPD do Dia 1):\n")
print(table(raw$fenotipo_D1, useNA = "ifany"))
cat("\nDistribuicao fenotipo_D5 (RASS/CAPD do Dia 5):\n")
print(table(raw$fenotipo_D5, useNA = "ifany"))


feno_map <- c("1" = "G1_free", "2" = "G2_delirium", "3" = "G3_coma")
raw$fenotiponeuro_label <- feno_map[as.character(raw$fenotiponeuro)]


lvl <- c("G1_free","G2_delirium","G3_coma")
raw$fenotipo_D1_f <- factor(raw$fenotipo_D1, levels = lvl)
raw$fenotipo_D5_f <- factor(raw$fenotipo_D5, levels = lvl)

trans_tab <- table(D1 = raw$fenotipo_D1_f, D5 = raw$fenotipo_D5_f, useNA = "ifany")
print(addmargins(trans_tab))

n_missing_d1 <- sum(is.na(raw$fenotipo_D1))
n_missing_d5 <- sum(is.na(raw$fenotipo_D5))
cat(sprintf("\nParticipantes sem avaliacao valida no D1: %d\n", n_missing_d1))
cat(sprintf("Participantes sem avaliacao valida no D5: %d\n", n_missing_d5))

write.csv(as.data.frame.matrix(addmargins(trans_tab)),
          file.path(out_dir, "table3_transition_matrix_D1xD5.csv"))

saveRDS(raw, file.path(out_dir, "phase1_reconstructed_data.rds"))
write.csv(raw, file.path(out_dir, "phase1_reconstructed_data.csv"), row.names = FALSE, na = "")

cat("\n=== FASE 1 CONCLUIDA ===\n")
cat("Dataset reconstruido salvo em output/phase1_reconstructed_data.rds e .csv\n")
cat("Matriz de transicao salva em output/table3_transition_matrix_D1xD5.csv\n")
