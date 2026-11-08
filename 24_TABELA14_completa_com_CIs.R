# 24_TABELA14_completa_com_CIs.R
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
    vals <- list(NfL = raw[["NFLD1,2"]], S100B = raw[["SB100D1,2"]],
                 `TNF-alfa` = raw[["TNF (pg/mL)D1"]], `IL-6` = raw[["IL-6 (pg/mL)D1"]])
  } else {
    feno <- raw$fenotipo_D5
    vals <- list(NfL = raw[["NFLD5,2"]], S100B = raw[["SB100D5,2"]],
                 `TNF-alfa` = raw[["tnfdeliriumd5"]], `IL-6` = raw[["il6 deliriumd5"]])
  }
  df <- data.frame(ID = raw$ID, phenotype = feno)
  for (nm in names(vals)) df[[nm]] <- vals[[nm]]
  df <- df[!is.na(df$phenotype), ]
  df$delirium <- ifelse(df$phenotype == "G2_delirium", 1, 0)
  df
}
df_d1 <- build_df("D1"); df_d5 <- build_df("D5")

analyze_full <- function(df, yvar) {
  df2 <- df[!is.na(df[[yvar]]), ]
  roc_obj <- roc(df2$delirium, df2[[yvar]], direction = "<", quiet = TRUE)
  apparent <- as.numeric(auc(roc_obj))
  ci_delong <- withCallingHandlers(ci.auc(roc_obj, method = "delong"),
                                    warning = function(w) invokeRestart("muffleWarning"))
  best <- coords(roc_obj, "best", best.method = "youden",
                  ret = c("threshold","sensitivity","specificity","ppv","npv"), transpose = FALSE)
  sens <- best$sensitivity[1]; spec <- best$specificity[1]
  ppv <- best$ppv[1]; npv <- best$npv[1]
  lr_pos <- sens/(1-spec); lr_neg <- (1-sens)/spec
  youden <- sens + spec - 1
  n <- nrow(df2)
  set.seed(20260818)

  boot_aucs <- numeric(R_BOOT); boot_sens <- numeric(R_BOOT); boot_spec <- numeric(R_BOOT)
  boot_ppv <- numeric(R_BOOT); boot_npv <- numeric(R_BOOT); boot_lrpos <- numeric(R_BOOT)
  boot_lrneg <- numeric(R_BOOT); boot_youden <- numeric(R_BOOT)
  for (i in seq_len(R_BOOT)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    bs <- df2[idx, ]
    if (length(unique(bs$delirium)) < 2) { boot_aucs[i] <- NA; next }
    r <- tryCatch(roc(bs$delirium, bs[[yvar]], direction = "<", quiet = TRUE), error = function(e) NULL)
    if (is.null(r)) { boot_aucs[i] <- NA; next }
    boot_aucs[i] <- as.numeric(auc(r))
    bb <- tryCatch(coords(r, x = best$threshold[1], input = "threshold",
                           ret = c("sensitivity","specificity","ppv","npv"), transpose = FALSE),
                    error = function(e) NULL)
    if (is.null(bb)) { boot_sens[i] <- NA; boot_spec[i] <- NA; boot_ppv[i] <- NA; boot_npv[i] <- NA
      boot_lrpos[i] <- NA; boot_lrneg[i] <- NA; boot_youden[i] <- NA; next }
    bsens <- bb$sensitivity[1]; bspec <- bb$specificity[1]
    boot_sens[i] <- bsens; boot_spec[i] <- bspec
    boot_ppv[i] <- bb$ppv[1]; boot_npv[i] <- bb$npv[1]
    boot_lrpos[i] <- if (bspec < 1) bsens/(1-bspec) else NA
    boot_lrneg[i] <- if (bspec > 0) (1-bsens)/bspec else NA
    boot_youden[i] <- bsens + bspec - 1
  }
  ci_of <- function(v) { v <- v[!is.na(v) & is.finite(v)]; if (length(v) < 100) return(c(NA,NA)); quantile(v, c(0.025,0.975)) }

  auc_ci_b <- ci_of(boot_aucs)
  sens_ci <- ci_of(boot_sens); spec_ci <- ci_of(boot_spec)
  ppv_ci <- ci_of(boot_ppv); npv_ci <- ci_of(boot_npv)
  lrpos_ci <- ci_of(boot_lrpos); lrneg_ci <- ci_of(boot_lrneg)
  youden_ci <- ci_of(boot_youden)

  data.frame(
    Biomarker = yvar, n = n,
    AUC = round(apparent,3), AUC_lo = round(ci_delong[1],3), AUC_hi = round(ci_delong[3],3),
    AUC_boot = round(mean(boot_aucs, na.rm=TRUE),3), AUC_boot_lo = round(auc_ci_b[1],3), AUC_boot_hi = round(auc_ci_b[2],3),
    Cutoff = round(best$threshold[1],2),
    Youden = round(youden,2), Youden_lo = round(youden_ci[1],2), Youden_hi = round(youden_ci[2],2),
    Sens = round(sens*100), Sens_lo = round(sens_ci[1]*100), Sens_hi = round(sens_ci[2]*100),
    Spec = round(spec*100), Spec_lo = round(spec_ci[1]*100), Spec_hi = round(spec_ci[2]*100),
    PPV = round(ppv*100), PPV_lo = round(ppv_ci[1]*100), PPV_hi = round(ppv_ci[2]*100),
    NPV = round(npv*100), NPV_lo = round(npv_ci[1]*100), NPV_hi = round(npv_ci[2]*100),
    LRpos = round(lr_pos,2), LRpos_lo = round(lrpos_ci[1],2), LRpos_hi = round(lrpos_ci[2],2),
    LRneg = round(lr_neg,2), LRneg_lo = round(lrneg_ci[1],2), LRneg_hi = round(lrneg_ci[2],2)
  )
}

bms <- c("NfL","S100B","TNF-alfa","IL-6")
res_d1 <- do.call(rbind, lapply(bms, function(b) analyze_full(df_d1, b)))
res_d5 <- do.call(rbind, lapply(bms, function(b) analyze_full(df_d5, b)))
res_d1$Day <- "D1"; res_d5$Day <- "D5"

tab <- rbind(res_d1, res_d5)
write.csv(tab, file.path(out_dir, "TABELA14_COMPLETA_com_todos_ICs.csv"), row.names = FALSE)
cat("=== TABELA 14 completa (todos os ICs) ===\n\n")
print(tab, row.names = FALSE)
