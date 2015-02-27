# Compute residuals

# Notes:
# 1. 

# Input:
# y = matrix of median adjusted quantiles (genes by samples)
# groups = groups factor / character
# lambda = shrinkage parameter
# levels = (optional) specify order of levels

# Output:
# returns a matrix of residuals

residuals = function (y, groups, lambda, levels=NULL) {
    
  fitModel = fitCoeffs(y, groups, lambda, levels)
  
  y - predicted(fitModel$design, fitModel$betas)
    
}