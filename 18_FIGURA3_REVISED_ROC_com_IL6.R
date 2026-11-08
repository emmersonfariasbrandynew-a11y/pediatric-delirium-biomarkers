# 18_FIGURA3_REVISED_ROC_com_IL6.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages(library(pROC))
options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
set.seed(20260818)
R_BOOT <- 5000

build_df <- function(day) {
  if (day == "D1") {
    feno <- raw$fenotipo_D1
    vals <- list(NfL = raw[["NFLD1,2"]], S100B = raw[["SB100D1,2"]], `IL-6` = raw[["IL-6 (pg/mL)D1"]])
  } else {
    feno <- raw$fenotipo_D5
    vals <- list(NfL = raw[["NFLD5,2"]], S100B = raw[["SB100D5,2"]], `IL-6` = raw[["il6 deliriumd5"]])
  }
  df <- data.frame(ID = raw$ID, phenotype = feno)
  for (nm in names(vals)) df[[nm]] <- vals[[nm]]
  df <- df[!is.na(df$phenotype), ]
  df$delirium <- ifelse(df$phenotype == "G2_delirium", 1, 0)
  df
}
df_d1 <- build_df("D1"); df_d5 <- build_df("D5")

analyze <- function(df, yvar) {
  df2 <- df[!is.na(df[[yvar]]), ]
  roc_obj <- roc(df2$delirium, df2[[yvar]], direction = "<", quiet = TRUE)
  apparent <- as.numeric(auc(roc_obj))
  ci_delong <- withCallingHandlers(ci.auc(roc_obj, method = "delong"),
                                    warning = function(w) invokeRestart("muffleWarning"))
  best <- coords(roc_obj, "best", best.method = "youden",
                  ret = c("threshold","sensitivity","specificity","ppv","npv"), transpose = FALSE)
  sens <- best$sensitivity[1]; spec <- best$specificity[1]
  lr_pos <- sens/(1-spec); lr_neg <- (1-sens)/spec
  n <- nrow(df2)
  set.seed(20260818)
  boot_aucs <- numeric(R_BOOT); boot_youden <- numeric(R_BOOT)
  for (i in seq_len(R_BOOT)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    bs <- df2[idx, ]
    if (length(unique(bs$delirium)) < 2) { boot_aucs[i] <- NA; boot_youden[i] <- NA; next }
    r <- tryCatch(roc(bs$delirium, bs[[yvar]], direction = "<", quiet = TRUE), error = function(e) NULL)
    if (is.null(r)) { boot_aucs[i] <- NA; boot_youden[i] <- NA; next }
    boot_aucs[i] <- as.numeric(auc(r))
    bb <- tryCatch(coords(r, x = best$threshold[1], input = "threshold",
                           ret = c("sensitivity","specificity"), transpose = FALSE), error = function(e) NULL)
    boot_youden[i] <- if (is.null(bb)) NA else (bb$sensitivity[1] + bb$specificity[1] - 1)
  }
  boot_aucs_ok <- boot_aucs[!is.na(boot_aucs)]
  boot_ci <- quantile(boot_aucs_ok, c(0.025, 0.975))
  boot_auc_mean <- mean(boot_aucs_ok)
  youden_apparent <- sens + spec - 1
  boot_y_ok <- boot_youden[!is.na(boot_youden)]
  youden_ci <- quantile(boot_y_ok, c(0.025, 0.975))
  list(roc = roc_obj,
       AUC = apparent, ci_lo = ci_delong[1], ci_hi = ci_delong[3],
       AUC_boot = boot_auc_mean, boot_lo = boot_ci[1], boot_hi = boot_ci[2],
       Cutoff = round(best$threshold[1], 2), Youden = round(youden_apparent, 2),
       youden_lo = youden_ci[1], youden_hi = youden_ci[2],
       SensPct = sens*100, SpecPct = spec*100, PPV = best$ppv[1]*100, NPV = best$npv[1]*100,
       LRpos = lr_pos, LRneg = lr_neg)
}

bms <- c("NfL","S100B","IL-6")
cores  <- c(NfL="black", S100B="darkgoldenrod", `IL-6`="purple")
pchs   <- c(NfL=17, S100B=16, `IL-6`=18)

res_d1 <- lapply(bms, function(b) analyze(df_d1, b)); names(res_d1) <- bms
res_d5 <- lapply(bms, function(b) analyze(df_d5, b)); names(res_d5) <- bms

draw_table <- function(res, bms) {
  plot.new()
  rows <- c("Variable","AUC (DeLong 95% CI)","Bootstrap AUC (95% CI)","Cutoff","Youden index (95% CI)",
            "LR+","LR-")
  fmt <- function(r) c(
    sprintf("%.3f (%.3f-%.3f)", r$AUC, r$ci_lo, r$ci_hi),
    sprintf("%.3f (%.3f-%.3f)", r$AUC_boot, r$boot_lo, r$boot_hi),
    sprintf("%.2f", r$Cutoff),
    sprintf("%.2f (%.2f-%.2f)", r$Youden, r$youden_lo, r$youden_hi),
    sprintf("%.2f", r$LRpos), sprintf("%.2f", r$LRneg))
  ncol_ <- length(bms) + 1
  xpos <- c(0.00, 0.24, 0.53, 0.80)
  ypos <- seq(0.95, 0.05, length.out = length(rows))
  for (j in seq_along(rows)) text(xpos[1], ypos[j], rows[j], adj = c(0,0.5), font = ifelse(j==1,2,1), cex = 0.85)
  for (i in seq_along(bms)) {
    r <- res[[bms[i]]]
    vals <- c(bms[i], fmt(r))
    for (j in seq_along(vals)) text(xpos[i+1], ypos[j], vals[j], adj = c(0,0.5), font = ifelse(j==1,2,1), cex = 0.85)
  }
}

png(file.path(out_dir, "FIGURA3_REVISED_com_IL6.png"), width = 3200, height = 3200, res = 300)
layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE), widths = c(0.8,1.5))
par(mar = c(4.5,4.5,3,1))
plot(res_d1$NfL$roc, col = cores["NfL"], lwd = 2, main = "Day 1", legacy.axes = TRUE)
mtext("A", side = 3, line = 1, adj = 0, font = 2, cex = 1.4)
for (b in bms[-1]) plot(res_d1[[b]]$roc, col = cores[b], lwd = 2, add = TRUE)
legend("bottomright", legend = bms, col = cores[bms], lwd = 2, bty = "n")
par(mar = c(1,0.5,3,0.5)); draw_table(res_d1, bms)
par(mar = c(4.5,4.5,3,1))
plot(res_d5$NfL$roc, col = cores["NfL"], lwd = 2, main = "Day 5", legacy.axes = TRUE)
mtext("B", side = 3, line = 1, adj = 0, font = 2, cex = 1.4)
for (b in bms[-1]) plot(res_d5[[b]]$roc, col = cores[b], lwd = 2, add = TRUE)
legend("bottomright", legend = bms, col = cores[bms], lwd = 2, bty = "n")
par(mar = c(1,0.5,3,0.5)); draw_table(res_d5, bms)
dev.off()

tiff(file.path(out_dir, "FIGURA3_REVISED_com_IL6.tif"), width = 13200, height = 13200, res = 1200, compression = "lzw")
layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE), widths = c(0.8,1.5))
par(mar = c(4.5,4.5,3,1))
plot(res_d1$NfL$roc, col = cores["NfL"], lwd = 2, main = "Day 1", legacy.axes = TRUE)
mtext("A", side = 3, line = 1, adj = 0, font = 2, cex = 1.4)
for (b in bms[-1]) plot(res_d1[[b]]$roc, col = cores[b], lwd = 2, add = TRUE)
legend("bottomright", legend = bms, col = cores[bms], lwd = 2, bty = "n")
par(mar = c(1,0.5,3,0.5)); draw_table(res_d1, bms)
par(mar = c(4.5,4.5,3,1))
plot(res_d5$NfL$roc, col = cores["NfL"], lwd = 2, main = "Day 5", legacy.axes = TRUE)
mtext("B", side = 3, line = 1, adj = 0, font = 2, cex = 1.4)
for (b in bms[-1]) plot(res_d5[[b]]$roc, col = cores[b], lwd = 2, add = TRUE)
legend("bottomright", legend = bms, col = cores[bms], lwd = 2, bty = "n")
par(mar = c(1,0.5,3,0.5)); draw_table(res_d5, bms)
dev.off()

cat("Salvo: output/FIGURA3_REVISED_com_IL6.tif e .png\n")
