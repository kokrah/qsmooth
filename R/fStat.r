# Perform row f-tests

# Notes:
# 1.

# Input:
# y = matrix of median adjusted quantiles (genes by samples)
# groups = groups factor / character
# levels = (optional) specify order of levels

# Output:
# A list: (1.) fstat (2.) df1 (3.) df2

fStat = function (y, groups, levels=NULL) {
  
  fit0 = fitCoeffs(y, groups=1, lambda=0, levels=NULL)
  pred0 = predicted(fit0$design, fit0$betas)
  
  fit1 = fitCoeffs(y, groups=groups, lambda=0, levels=levels)
  pred1 = predicted(fit1$design, fit1$betas)
  
  SSR = rowSums((pred1 - pred0)^2)
  SSE = rowSums((y - pred1)^2)
  
  n = ncol(y)
  k0 = ncol(fit0$design)
  k1 = ncol(fit1$design)
  
  df1 = k1 - k0
  df2 = n - k1
  
  MSR = SSR / df1
  MSE = SSE / df2
  
  list(fstat=MSR/MSE, df1=df1, df2=df2)
}
