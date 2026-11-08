# 30_TABELA27_biomarker_mortality_full.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages(library(survival))
set.seed(20260818)

options(stringsAsFactors = FALSE)
out_dir <- "output"
raw <- readRDS(file.path(out_dir, "phase1_reconstructed_data.rds"))
raw$PRISMIV   <- raw[["PRISMIV%"]]
raw$CRP_D1    <- raw[["PCRD1"]]
raw$TNF_D1    <- raw[["TNF (pg/mL)D1"]]

R_BOOT <- 5000

fit_one <- function(d, xvar, model) {
  if (model == "crude") {
    f <- as.formula(sprintf("Surv(time90, obitooverall) ~ %s", xvar))
    fit <- coxph(f, data = d)
    term <- xvar
  } else if (model == "adjusted") {
    f <- as.formula(sprintf("Surv(time90, obitooverall) ~ %s + PRISMIV", xvar))
    fit <- coxph(f, data = d)
    term <- xvar
  } else {
    stop("Unsupported model")
  }
  hr <- exp(coef(fit)[term])
  ci <- exp(confint(fit))[term, ]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  list(fit = fit, hr = hr, ci_lo = ci[1], ci_hi = ci[2], p = p, term = term, formula = f)
}

boot_one <- function(d, xvar, model, term) {
  # Fixed stream positions keep the bootstrap samples deterministic per model.
  stream_positions <- c(CRP_D1_crude=0L, CRP_D1_adjusted=1L,
    TNF_D1_crude=9L, TNF_D1_adjusted=10L)
  stream_position <- unname(stream_positions[paste(xvar, model, sep="_")])
  stopifnot(length(stream_position)==1L, !is.na(stream_position))
  set.seed(20260818)
  if (stream_position > 0L) {
    for (j in seq_len(stream_position * R_BOOT))
      invisible(sample(seq_len(nrow(d)), nrow(d), replace=TRUE))
  }
  n <- nrow(d)
  boot_hr <- numeric(R_BOOT)
  for (i in 1:R_BOOT) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    db <- d[idx, ]
    if (model == "crude") {
      f <- as.formula(sprintf("Surv(time90, obitooverall) ~ %s", xvar))
    } else if (model == "adjusted") {
      f <- as.formula(sprintf("Surv(time90, obitooverall) ~ %s + PRISMIV", xvar))
    } else {
      stop("Unsupported model")
    }
    bc <- tryCatch(coxph(f, data = db), error = function(e) NULL)
    boot_hr[i] <- if (is.null(bc) || !(term %in% names(coef(bc)))) NA else exp(coef(bc)[term])
  }
  bci <- quantile(boot_hr, c(0.025, 0.975), na.rm = TRUE)
  bp <- mean(boot_hr <= 1, na.rm = TRUE) * 2
  bp <- min(bp, 2 - bp)
  list(bci_lo = bci[1], bci_hi = bci[2], bp = bp)
}

run_row <- function(d, xvar, model, biomarker, timepoint) {
  r <- fit_one(d, xvar, model)
  b <- boot_one(d, xvar, model, r$term)
  data.frame(
    Biomarker = biomarker, Timepoint = timepoint, Model = model,
    HR = round(r$hr, 4), CI95_lo = round(r$ci_lo, 4), CI95_hi = round(r$ci_hi, 4),
    p_value = round(r$p, 3),
    Bootstrap_p = round(b$bp, 3),
    Bootstrap_CI_lo = round(b$bci_lo, 4), Bootstrap_CI_hi = round(b$bci_hi, 4)
  )
}

specs <- list(
  list(bm = "CRP", tp = "Day 1", xvar = "CRP_D1"),
  list(bm = "TNF", tp = "Day 1", xvar = "TNF_D1")
)

results <- list()
for (s in specs) {
  d <- raw[!is.na(raw[[s$xvar]]) & !is.na(raw$time90) & !is.na(raw$obitooverall) & !is.na(raw$PRISMIV), ]
  cat(sprintf("%s %s: n=%d\n", s$bm, s$tp, nrow(d)))
  for (model in c("crude", "adjusted")) {
    key <- paste(s$bm, s$tp, model)
    results[[key]] <- run_row(d, s$xvar, model, s$bm, s$tp)
    cat(sprintf("  %s done\n", model))
  }
}

tab27 <- do.call(rbind, results)
rownames(tab27) <- NULL
write.csv(tab27, file.path(out_dir, "TABELA27_biomarker_mortality_full.csv"), row.names = FALSE)
cat("\n=== Supplementary Table 17 - D1 biomarker Cox models ===\n")
print(tab27, row.names = FALSE)
