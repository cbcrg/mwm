anova.censReg <- function(mod1, mod0, ...)
 {
  lC <- logLik(mod1)
  dfC <- attr(lC, "df")
  lR <- logLik(mod0)
  dfR <- attr(lR, "df")
  df <- dfC - dfR
  L <- -2*(lR - lC)
  pval <- pchisq(L, df, lower=FALSE)
  attr(pval, "df") <- df
  unclass(pval)
}
