#' Quantile shrinkage normalization
#' 
#' @param exprs for counts use log2(raw counts + 1)), for MA use log2(raw intensities)
#' @param groups groups to which samples belong (character vector)
#' @param refType type of quantile reference (defualt="mean")
#' @param groupLoc type of group location estimate (default="mean")
#' @param window window size for running median
#' @param verbose show info
#' @param groupCol character vector indicating the group color for each sample  
#' @param plot plot weights? (default=FALSE)
#' @param ... plot options
#' @export 
qshrink = function (exprs, groups, refType="mean", groupLoc="mean", window=99,
                    verbose=FALSE, groupCol=NULL, plot=FALSE, ...) {
  # 1. Compute quantile stats
  res = qstats(exprs, groups, refType=refType, groupLoc=groupLoc, window=window)
  
  QBETAS = res$QBETAS
  Qref = res$Qref
  X = res$X
  w = res$smoothWeights
  
  # 2. Compute weighted quantiles
  wQBETAS = QBETAS * (1 - w)
  
  wQBETAS = X %*% t(wQBETAS)
  
  wQref = Qref * w
  
  wQref = matrix(rep(1, nrow(X)), ncol=1) %*% t(wQref)
  
  normExprs = t(wQBETAS + wQref)
  
  # 3. Re-order
  RANKS = t(matrixStats::colRanks(exprs, ties.method = "average"))
  
  for (k in 1:ncol(normExprs)) {
    x = normExprs[, k]
    normExprs[, k] = x[RANKS[, k]]
  }
  
  # 4. Average ties  
  normExprs = aveTies(RANKS, normExprs)
  
  rownames(normExprs) = rownames(exprs)
  colnames(normExprs) = colnames(exprs)
  
  # 5. Plot weights
  if (plot) {
    
    oldpar = par(mar=c(4, 4, 1.5, 0.5))
    
    lq = length(Qref)
    u = (1:lq - 0.5) / lq
    
    if (length(u) > 10000) { # do not plot more than 10000 points
      
      sel = sample(1:lq, 10000)
        
      plot(u[sel], w[sel], pch=".", main="Quantile reference weights",
           xlab="u (normalized gene ranks)", ylab="Weight", ylim=c(0, 1), ...)
      
    }else{
      
      plot(u, w, pch=".", main="Quantile reference weights",
           xlab="u (normalized gene ranks)", ylab="Weight", ylim=c(0, 1), ...)
      
    }

    abline(h=0.5, v=0.5, col="red", lty=2)
    
    par(oldpar)
  
  }
  
  normExprs
}
