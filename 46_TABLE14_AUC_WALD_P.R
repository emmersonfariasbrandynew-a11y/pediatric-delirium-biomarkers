# Two-sided Wald tests for H0: AUC = 0.5, using DeLong variance.
# Apparent AUCs are exploratory; these tests do not correct selection optimism.
if (!requireNamespace('pROC', quietly = TRUE)) stop('Install pROC.')
d <- readRDS(file.path('output','phase1_reconstructed_data.rds'))
rows <- list()
for(day in c('D1','D5')) {
 y <- as.integer(d[[paste0('fenotipo_',day)]] == 'G2_delirium')
 cols <- if(day == 'D1') c(NfL='NFLD1,2',S100B='SB100D1,2',`TNF-alfa`='TNF (pg/mL)D1',`IL-6`='IL-6 (pg/mL)D1') else c(NfL='NFLD5,2',S100B='SB100D5,2',`TNF-alfa`='tnfdeliriumd5',`IL-6`='il6 deliriumd5')
 for(b in names(cols)) {
  x <- d[[cols[[b]]]]; keep <- complete.cases(y,x)
  r <- pROC::roc(y[keep],x[keep],direction='<',quiet=TRUE)
  a <- as.numeric(pROC::auc(r)); v <- as.numeric(pROC::var(r,method='delong'))
  if (!is.finite(v) || v <= 0) stop('Nonpositive or nonfinite DeLong variance.')
  z <- (a-0.5)/sqrt(v); p <- 2*pnorm(-abs(z))
  ci <- as.numeric(pROC::ci.auc(r,method='delong'))
  rows[[length(rows)+1]] <- data.frame(Biomarker=b,Day=day,n_cases=sum(y[keep]==1),n_controls=sum(y[keep]==0),AUC=a,AUC_lo=ci[1],AUC_hi=ci[3],DeLong_variance=v,z=z,p_value=p,p_display=if(p<0.001)'p < 0.001' else sprintf('p = %.3f',p),method='Two-sided Wald test; DeLong variance; H0 AUC=0.5')
 }
}
result <- do.call(rbind,rows)
stopifnot(nrow(result)==8, all(result$p_value>=0 & result$p_value<=1))
write.csv(result,file.path('output','TABLE14_AUC_WALD_P.csv'),row.names=FALSE)
print(result[,c('Biomarker','Day','AUC','p_value','p_display')],row.names=FALSE)
