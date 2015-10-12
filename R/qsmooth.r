#' Quantile shrinkage normalization
#' 
#' @param exprs for counts use log2(raw counts + 1)), for MA use log2(raw intensities)
#' @param groups groups to which samples belong (character vector)
#' @param norm.factors scaling normalization factors
#' @param plot plot weights? (default=FALSE)
#' @param window window size for running median as a fraction on the number of rows of exprs
#' @export 
qsmooth = function (exprs, groups, norm.factors=NULL, plot=FALSE, window=0.05) {
  
  # Scale normalization step
  if (is.null(norm.factors)) {
    dat = exprs
  }else{
    dat = t(t(exprs) - norm.factors)
  }
  
  # Compute quantile stats
  qs = qstats(dat, groups, window=window)
  Qref = qs$Qref 
  Qhat = qs$Qhat  
  w = qs$smoothWeights
  
  # Weighted quantiles
  normExprs = w * Qref + (1 - w) * Qhat
  
  # Re-order
  RANKS = apply(exprs, 2, rank, ties.method="average")
  
  for (k in 1:ncol(normExprs)) {
    x = normExprs[, k]
    normExprs[, k] = x[RANKS[, k]]
  }
  
  # Average ties  
  normExprs = aveTies(RANKS, normExprs)
  
  # Plot weights
  if (plot) {
    
    oldpar = par(mar=c(4, 4, 1.5, 0.5))
    
    lq = length(Qref)
    u = (1:lq - 0.5) / lq
    
    if (length(u) > 10000) { # do not plot more than 10000 points
      
      sel = sample(1:lq, 10000)
      
      plot(u[sel], w[sel], pch=".", main="Quantile reference weights",
           xlab="u (normalized gene ranks)", ylab="Weight", ylim=c(0, 1))
      
    }else{
      
      plot(u, w, pch=".", main="Quantile reference weights",
           xlab="u (normalized gene ranks)", ylab="Weight", ylim=c(0, 1))
      
    }
    
    abline(h=0.5, v=0.5, col="red", lty=2)
    
    par(oldpar)
    
  }
  
  rownames(normExprs) = rownames(exprs)
  colnames(normExprs) = colnames(exprs)
  normExprs
}
