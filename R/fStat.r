#' Perform row f-tests
#' Notes:
#' Output:
#' A list: (1.) fstat (2.) df1 (3.) df2
#' @param alpha matrix of median adjusted quantiles (genes by samples)
#' @param groupFactor groups factor / character
#' @param levels (optional) specify order of levels
#' @export
fStat = function (alpha, groupFactor, levels=NULL) {
  
  fit0 = fitCoeffs(alpha, groupFactor=1, lambda=0, levels=NULL)
  pred0 = t(predicted(fit0$design, fit0$betas))
  
  fit1 = fitCoeffs(alpha, groupFactor=groupFactor, lambda=0, levels=levels)
  pred1 = t(predicted(fit1$design, fit1$betas))
  
  SSR = rowSums((pred1 - pred0)^2)
  SSE = rowSums((alpha - pred1)^2)
  
  n = ncol(alpha)
  k0 = ncol(fit0$design)
  k1 = ncol(fit1$design)
  
  df1 = k1 - k0 # between
  df2 = n - k1 # within
  
  MSR = SSR / df1
  MSE = SSE / df2
  
  list(fstat=MSR/MSE, df1=df1, df2=df2)
}
