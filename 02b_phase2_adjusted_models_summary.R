# 02b_phase2_adjusted_models_summary.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

out_dir <- "output"
obj <- readRDS(file.path(out_dir, "phase2_models.rds"))

sink(file.path(out_dir, "adjusted_models_summary.txt"), split = TRUE)

for (bm in c("nfl_fits", "s100b_fits")) {
  label <- ifelse(bm == "nfl_fits", "NfL", "S100B")
  fits <- obj[[bm]]
  for (mname in c("M1_PRISM", "M2_sedation", "M3_PRISM_sedation")) {
    cat(sprintf("\n\n================ %s : %s ================\n", label, mname))
    print(summary(fits[[mname]]))
    cat(sprintf("Singular: %s\n", isSingular(fits[[mname]], tol = 1e-4)))
  }
}

sink()
cat("Saved to output/adjusted_models_summary.txt\n")
