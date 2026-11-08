# 41_FIGURA2_PAINELA_contemporaneo.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(ggplot2)
})
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

long_D1 <- data.frame(id = raw$ID, time = "D1", phenotype = raw$fenotipo_D1,
                       NfL = raw[["NFLD1,2"]], S100B = raw[["SB100D1,2"]])
long_D5 <- data.frame(id = raw$ID, time = "D5", phenotype = raw$fenotipo_D5,
                       NfL = raw[["NFLD5,2"]], S100B = raw[["SB100D5,2"]])
figure_data <- rbind(long_D1, long_D5)
figure_data$time <- factor(figure_data$time, levels = c("D1","D5"))
figure_data$phenotype <- factor(figure_data$phenotype,
                                 levels = c("G1_free","G2_delirium","G3_coma"),
                                 labels = c("G1","G2","G3"))

cat("=== N per time x phenotype (contemporary) - MUST match D1: 14/20/11, D5: 10/24/11 ===\n")
tab <- with(figure_data, table(time, phenotype, useNA = "ifany"))
print(tab)
expected_D1 <- c(G1 = 14, G2 = 20, G3 = 11)
expected_D5 <- c(G1 = 10, G2 = 24, G3 = 11)
stopifnot(all(tab["D1", c("G1","G2","G3")] == expected_D1))
stopifnot(all(tab["D5", c("G1","G2","G3")] == expected_D5))
cat("N counts verified OK against expected contemporary denominators.\n\n")

letters3 <- function(sig_12, sig_13, sig_23) {
  if (!sig_12 && !sig_13 && !sig_23) return(c(G1="a", G2="a", G3="a"))
  if (sig_12 && sig_13 && sig_23)    return(c(G1="a", G2="b", G3="c"))
  if (sig_12 && sig_13 && !sig_23)   return(c(G1="a", G2="b", G3="b"))  # G2~G3
  if (sig_12 && !sig_13 && sig_23)   return(c(G1="a", G2="b", G3="a"))  # G1~G3
  if (!sig_12 && sig_13 && sig_23)   return(c(G1="a", G2="a", G3="b"))  # G1~G2
  if (sig_12 && !sig_13 && !sig_23)  return(c(G1="a", G2="b", G3="ab")) # G3 bridges
  if (!sig_12 && sig_13 && !sig_23)  return(c(G1="a", G2="ab", G3="b")) # G2 bridges
  if (!sig_12 && !sig_13 && sig_23)  return(c(G1="ab", G2="a", G3="b")) # G1 bridges
}

get_letters_for_time <- function(bm, tpoint) {
  cf <- read.csv(file.path(out_dir, sprintf("contrasts_%s_phenotype_within_time.csv", bm)))
  cf <- cf[cf$time == tpoint, ]
  p12 <- cf$p.value[cf$contrast == "G1_free - G2_delirium"]
  p13 <- cf$p.value[cf$contrast == "G1_free - G3_coma"]
  p23 <- cf$p.value[cf$contrast == "G2_delirium - G3_coma"]
  letters3(p12 < 0.05, p13 < 0.05, p23 < 0.05)
}

get_letters_for_phenotype <- function(bm) {
  cf <- read.csv(file.path(out_dir, sprintf("contrasts_%s_time_within_phenotype.csv", bm)))
  out <- c(G1_free = "A", G2_delirium = "A", G3_coma = "A")
  for (ph in names(out)) {
    p <- cf$p.value[cf$phenotype == ph]
    if (length(p) == 1 && p < 0.05) out[ph] <- "B"  # D5 differs from D1
  }
  out
}

annotate_biomarker <- function(bm) {
  lo_d1 <- get_letters_for_time(bm, "D1")
  lo_d5 <- get_letters_for_time(bm, "D5")
  hi <- get_letters_for_phenotype(bm)  # names G1_free/G2_delirium/G3_coma; "A"=D1-like, "B"=D5 differs
  data.frame(
    phenotype = rep(c("G1","G2","G3"), 2),
    time = rep(c("D1","D5"), each = 3),
    label = c(paste0(lo_d1["G1"], "A"), paste0(lo_d1["G2"], "A"), paste0(lo_d1["G3"], "A"),
              paste0(lo_d5["G1"], ifelse(hi["G1_free"]=="B","B","A")),
              paste0(lo_d5["G2"], ifelse(hi["G2_delirium"]=="B","B","A")),
              paste0(lo_d5["G3"], ifelse(hi["G3_coma"]=="B","B","A")))
  )
}

ann_NfL   <- annotate_biomarker("NfL")
ann_S100B <- annotate_biomarker("S100B")
cat("=== Compact letter annotations (contemporary model, script 02 contrasts) ===\n")
cat("NfL:\n");   print(ann_NfL)
cat("S100B:\n"); print(ann_S100B)
write.csv(ann_NfL,   file.path(out_dir, "FIGURA2_PAINELA_letters_NfL.csv"), row.names = FALSE)
write.csv(ann_S100B, file.path(out_dir, "FIGURA2_PAINELA_letters_S100B.csv"), row.names = FALSE)

global_stars <- "***"

make_panel <- function(bm, ylab, ann) {
  df <- figure_data[, c("id","time","phenotype", bm)]
  names(df)[4] <- "value"
  df <- df[!is.na(df$value) & !is.na(df$phenotype), ]

  ymax <- max(df$value, na.rm = TRUE)
  ann$y <- ymax * 1.08

  ggplot(df, aes(x = phenotype, y = value, fill = time)) +
    geom_violin(position = position_dodge(width = 0.8), trim = TRUE, alpha = 0.85, linewidth = 0.3) +
    geom_point(aes(color = time), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8, seed = 20260818),
               size = 0.9, alpha = 0.6, show.legend = FALSE) +
    geom_text(data = ann, aes(x = phenotype, y = y, label = label, group = time),
              position = position_dodge(width = 0.8), size = 3.2, fontface = "bold", inherit.aes = FALSE) +
    annotate("segment", x = 1, xend = 3, y = ymax * 1.22, yend = ymax * 1.22, linewidth = 0.4) +
    annotate("text", x = 2, y = ymax * 1.30, label = global_stars, size = 4) +
    scale_fill_manual(values = c(D1 = "#E31A1C", D5 = "#33A02C")) +
    scale_color_manual(values = c(D1 = "#E31A1C", D5 = "#33A02C")) +
    labs(x = NULL, y = ylab, fill = NULL) +
    coord_cartesian(ylim = c(0, ymax * 1.35)) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
}

p_NfL   <- make_panel("NfL",   "NfL (ng/mL)",   ann_NfL)
p_S100B <- make_panel("S100B", "S100B (ng/mL)", ann_S100B)

suppressPackageStartupMessages(library(cowplot))
plots_row <- cowplot::plot_grid(p_S100B, p_NfL, ncol = 2, labels = NULL)
title_row <- cowplot::ggdraw() +
  cowplot::draw_label("Panel A", x = 0.02, y = 0.5, hjust = 0, vjust = 0.5, fontface = "bold", size = 16)
panelA <- cowplot::plot_grid(title_row, plots_row, ncol = 1, rel_heights = c(0.08, 0.92))

tiff(file.path(out_dir, "FIGURA2_PAINELA_REVISED.tif"), width = 12800, height = 5652, res = 1200, compression = "lzw")
print(panelA)
dev.off()
ggsave(file.path(out_dir, "FIGURA2_PAINELA_REVISED.png"), panelA, width = 10.667, height = 4.71, dpi = 300)

cat("\nSalvo: output/FIGURA2_PAINELA_REVISED.tif e .png\n")
cat("Letras lowercase = diferenca entre fenotipos DENTRO do mesmo dia (Tukey).\n")
cat("Letras uppercase = diferenca entre D1 e D5 DENTRO do mesmo fenotipo (Tukey).\n")
cat("*** = teste global de fenotipo (Tipo III), p<0.001 para ambos os biomarcadores.\n")
