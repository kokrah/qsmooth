#' Compute quantile statistics
#' 
#' @param exprs for counts use log2(raw counts + 1)), for MA use log2(raw intensities)
#' @param groups groups to which samples belong (character vector)
#' @param refType type of quantile reference (defualt="mean")
#' @param groupLoc type of group location estimate (default="mean")
#' @param window window size for running median
#' @export 
qstats = function (exprs, groups, refType="mean", groupLoc="mean", window=99) {
  
  # 1. Compute sample quantiles
  Q = apply(exprs, 2, sort) 
  
  # 2. Compute quantile reference
  if (refType == "median") {
    Qref = matrixStats::rowMedians(Q)  
  }
    
  if (refType == "mean") {
    Qref = rowMeans(Q)  
  }
  
  # 3. Compute group location estimates
  uGroups = unique(groups)
  
  if (groupLoc=="mean") {
    QBETAS = c()
    
    for (g in uGroups) {
      QBETAS = cbind(QBETAS, rowMeans(Q[, g==groups]))
    }
    
    colnames(QBETAS) = uGroups
  }
  
  if (groupLoc=="median") {
    QBETAS = c()
    
    for (g in uGroups) {
       QBETAS = cbind(QBETAS, matrixStats::rowMedians(Q[, g==groups]))
    }

    colnames(QBETAS) = uGroups
  }
  
  # 4. Compute group variance estimates
  if (groupLoc=="mean") {
    TAU = matrixStats::rowVars(QBETAS)  
  }
  
  if (groupLoc=="median") {
    TAU = (matrixStats::rowMads(QBETAS))^2
  }
  
  # 5. Compute within group variance estimates
  if (groupLoc=="mean") {
    SIGMA = c()
    
    for (g in uGroups) {
      SIGMA = cbind(SIGMA, matrixStats::rowVars(Q[, g==groups]))
    }
    
    SIGMA = rowMeans(SIGMA)
  }

  if (groupLoc=="median") {
    SIGMA = c()
    
    for (g in uGroups) {
      SIGMA = cbind(SIGMA, (matrixStats::rowMads(Q[, g==groups]))^2)
    }
    
    SIGMA = matrixStats::rowMedians(SIGMA)
  }
  
  
  # 6. Compute weights 
  roughWeights = SIGMA / (SIGMA + TAU)
  
  # if both SIGMA and TAU are 0, put 100% weigth on Ref.
  roughWeights[SIGMA < 10^(-6) & TAU < 10^(-6)] = 1 

#   u = (1:nrow(Q) - 0.5) / nrow(Q)
#   l = lowess(u, roughWeights, f=span)
#   f = approxfun(l, rule = 2)
# #   smoothWeights = f(u)
  smoothWeights = runmed(roughWeights, k = window, endrule="constant")
#   smoothWeights[smoothWeights > 1] = 1
  
  # 7. List results
  list(Q=Q, Qref=Qref, QBETAS=QBETAS, TAU=TAU, SIGMA=SIGMA,
       roughWeights=roughWeights, smoothWeights=smoothWeights,
       X = model.matrix(~0+factor(groups, levels=uGroups)))
  
}