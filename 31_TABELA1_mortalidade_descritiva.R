# 31_TABELA1_mortalidade_descritiva.R
# See README.md and script_output_map.csv for inputs and outputs.
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))

n <- nrow(raw)
n_death <- sum(raw$obitooverall == 1, na.rm = TRUE)
cat(sprintf("Overall in-hospital mortality, n/N: %d/%d\n", n_death, n))

n_death30 <- sum(raw[["OBITO<30"]] == 1, na.rm = TRUE)
cat(sprintf("In-hospital mortality <=30 days, n/N: %d/%d\n", n_death30, n))

n_alive30 <- n - n_death30
nonsurv <- raw[raw$obitooverall == 1, ]
n_death_gt30 <- sum(nonsurv[["OBITO<30"]] == 0, na.rm = TRUE)
cat(sprintf("In-hospital mortality >30 days, n/N (among patients alive at 30 days): %d/%d\n",
            n_death_gt30, n_alive30))

med_ttd <- median(nonsurv$time90)
q1_ttd  <- quantile(nonsurv$time90, 0.25)
q3_ttd  <- quantile(nonsurv$time90, 0.75)
cat(sprintf("Time to death among non-survivors, days, median (IQR): %.0f (%.0f-%.0f)\n",
            med_ttd, q1_ttd, q3_ttd))

sorted_death_times <- sort(nonsurv$time90)
cat("\nAll death times (time90), sorted, illustrating the bimodal pattern:\n")
print(sorted_death_times)

result <- data.frame(
  Row = c("Overall in-hospital mortality, n/N",
          "Time to death among non-survivors, days, median (IQR)",
          "In-hospital mortality <=30 days, n/N",
          "In-hospital mortality >30 days, n/N"),
  Value = c(sprintf("%d/%d", n_death, n),
            sprintf("%.0f (%.0f-%.0f)", med_ttd, q1_ttd, q3_ttd),
            sprintf("%d/%d", n_death30, n),
            sprintf("%d/%d", n_death_gt30, n_alive30))
)
write.csv(result, file.path(out_dir, "TABELA1_mortalidade_descritiva.csv"), row.names = FALSE)
cat("\n=== Table 1 mortality rows - reproduced ===\n")
print(result, row.names = FALSE)
