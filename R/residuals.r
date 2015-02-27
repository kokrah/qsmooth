#' Compute residuals
#' Notes:
#' 1. 
#' Output:
#' returns a matrix of residuals
#' @param y matrix of median adjusted quantiles (genes by samples)
#' @param groups groups factor / character
#' @param lambda shrinkage parameter
#' @param levels (optional) specify order of levels
#' @export
residuals = function (y, groups, lambda, levels=NULL) {
    
  fitModel = fitCoeffs(y, groups, lambda, levels)
  
  y - predicted(fitModel$design, fitModel$betas)
    
}