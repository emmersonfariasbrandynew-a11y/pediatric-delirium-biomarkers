# 13_pacote_tabela_p_raw_q_FDR.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages(library(lme4))
out_dir <- "output"

mods <- readRDS(file.path(out_dir, "phase2_models.rds"))

get_omnibus <- function(fit, label) {
  an <- anova(fit, type = 3)
  p_pheno <- if ("phenotype" %in% rownames(an)) an["phenotype","Pr(>F)"] else NA
  p_inter <- if (any(grepl(":", rownames(an)))) an[grep(":", rownames(an)),"Pr(>F)"] else NA
  data.frame(biomarker = label, grupo = "Primario (pre-especificado)",
             n_obs = nobs(fit), singular = isSingular(fit, tol = 1e-4),
             p_phenotype_raw = p_pheno, p_interaction_raw = p_inter,
             familia = "1 - Primaria (NfL/S100B, sem FDR)")
}

fam1 <- rbind(
  get_omnibus(mods$nfl_fits$M0_raw,   "NfL"),
  get_omnibus(mods$s100b_fits$M0_raw, "S100B")
)
fam1$q_FDR_phenotype   <- NA_real_
fam1$q_FDR_interaction <- NA_real_

cat("=== FAMILIA 1 - PRIMARIA (NfL, S100B) - SEM CORRECAO DE FDR ===\n")
print(fam1[, c("biomarker","p_phenotype_raw","p_interaction_raw")], row.names = FALSE)

fam2 <- read.csv(file.path(out_dir, "table13_14_inflammatory_panel_FDR.csv"))
fam2$familia <- "2 - Exploratoria (16 biomarcadores, BH-FDR)"

cat("\n=== FAMILIA 2 - EXPLORATORIA (16 biomarcadores) - COM BH-FDR ===\n")
cat(sprintf("(reaproveitado de table13_14_inflammatory_panel_FDR.csv, %d biomarcadores)\n\n", nrow(fam2)))

cols <- c("biomarker","grupo","familia","n_obs","singular",
          "p_phenotype_raw","q_FDR_phenotype","p_interaction_raw","q_FDR_interaction")
master <- rbind(fam1[, cols], fam2[, cols])

cat("\n=== TABELA CONSOLIDADA - P-RAW / Q-FDR (TODOS OS BIOMARCADORES) ===\n\n")
print(master, row.names = FALSE, digits = 4)

write.csv(master, file.path(out_dir, "table_MASTER_p_raw_q_FDR.csv"), row.names = FALSE)

n_fam1_sig <- sum(fam1$p_phenotype_raw < 0.05, na.rm = TRUE)
n_fam2_sig_raw <- sum(fam2$p_phenotype_raw < 0.05, na.rm = TRUE)
n_fam2_sig_fdr <- sum(fam2$q_FDR_phenotype < 0.05, na.rm = TRUE)

cat(sprintf("\nResumo: Familia 1 (n=2) com p bruto<0.05 (efeito fenotipo): %d/2\n", n_fam1_sig))
cat(sprintf("Familia 2 (n=%d) com p bruto<0.05: %d | com q-FDR<0.05: %d\n",
            nrow(fam2), n_fam2_sig_raw, n_fam2_sig_fdr))


cat("\n=== CONCLUIDO ===\n")
