#' Internal function: average ties
#' 
#' @param ranks ranks
#' @param y values 
ave.ties = function (ranks, y) {
  tab = table(ranks)
  sel = tab > 1
  
  if (sum(sel) != 0) {
    
    ties = as.numeric(names(tab[sel]))
    
    for (k in ties) {
      sel = ranks==k
      y[sel] = mean(y[sel])
    }
    
  }

  y    
}

#' Internal function: average ties
#' 
#' @param RANKS matrix of ranks 
#' @param normExprs normalized values 
aveTies = function (RANKS, normExprs) {
  normExprs
  
  for (k in ncol(RANKS)) {
    normExprs[,k] = ave.ties(RANKS[,k], normExprs[,k])
  }
  
  normExprs
}
