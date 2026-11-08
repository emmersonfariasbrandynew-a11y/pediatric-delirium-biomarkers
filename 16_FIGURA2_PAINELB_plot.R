# 16_FIGURA2_PAINELB_plot.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(grid)})
out_dir <- "output"
df <- read.csv(file.path(out_dir, "FIGURA2_PAINELB_forest_data_cohens_d.csv"), stringsAsFactors = FALSE)

model_order <- c("Raw","Adj. PRISM IV","Adj. PRISMIV+Sedation","Adj. Sedation",
                  "Sepsis-adjusted",
                  "Delta (D5-D1) Cohort")
model_labels <- c("Raw","Adj. PRISM IV","Adj. PRISM IV+Sedation","Adj. Sedation",
                   "Sepsis-adjusted",
                   "Delta (D5-D1) Cohort")
df$model <- factor(df$model, levels = rev(model_order), labels = rev(model_labels))
df$contrast_label <- factor(df$contrast_label, levels = c("G1 - G2","G2 - G3","G1 - G3"),
                             labels = c("G1 vs G2","G2 vs G3","G1 vs G3"))
df$biomarker <- factor(df$biomarker, levels = c("S100B","NfL"),
                        labels = c("S100B (log-transformed)","NfL (log-transformed)"))

pd <- position_dodge(width = 0.6)
p <- ggplot(df, aes(x = cohens_d, y = model, color = contrast_label, shape = contrast_label)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0, position = pd, linewidth = 0.6) +
  geom_point(position = pd, size = 2.6) +
  scale_color_manual(values = c("G1 vs G2" = "red", "G2 vs G3" = "blue", "G1 vs G3" = "darkgreen")) +
  scale_shape_manual(values = c("G1 vs G2" = 16, "G2 vs G3" = 17, "G1 vs G3" = 15)) +
  facet_wrap(~biomarker, scales = "free_x") +
  labs(x = "Standardized mean difference (Cohen's d, 95% CI)", y = NULL, color = NULL, shape = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold"), legend.position = "right",
        panel.grid.minor = element_blank(), plot.margin = margin(30, 10, 10, 10))

p <- cowplot::ggdraw(p) +
  cowplot::draw_label("Panel B", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
                       fontface = "bold", size = 16)

tiff(file.path(out_dir, "FIGURA2_PAINELB_forest_REVISED.tif"), width = 12800, height = 6400, res = 1200, compression = "lzw")
print(p)
dev.off()

ggsave(file.path(out_dir, "FIGURA2_PAINELB_forest_REVISED.png"), p, width = 10.667, height = 5.333, dpi = 300)
cat("Salvo: output/FIGURA2_PAINELB_forest_REVISED.tif e .png\n")
