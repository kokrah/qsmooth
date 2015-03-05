#' Fit alpha group coefficients
#' Notes:
#' 1. Currently using ridge regression
#' Output:
#' returns a list: (1.) a matrix of fitted group means for each gene (2.) design matrix
#'
#' @param alpha matrix of reference adjusted quantiles (genes by samples)
#' @param groupFactor groups factor / character
#' @param lambda shrinkage parameter
#' @param levels (optional) specify order of levels
#' @export
fitCoeffs = function (alpha, groupFactor, lambda, levels=NULL) {
  
  Y = as.matrix(alpha)
  
  if (length(groupFactor)==1) {
    
    X = matrix(1, ncol(Y))
    colnames(X) = "Average"
    
  }else{
   
    
    if (is.null(levels)) {
      
      groupFactor = factor(groupFactor)  
      
    }else{
      
      groupFactor = factor(groupFactor, levels=levels)
      
    }
    
    X = model.matrix(~ 0 + groupFactor)
    colnames(X) = levels(groupFactor)
    
  }
  
  if (length(lambda) > 1) {
    
    n = ncol(Y)
    betas = apply(cbind(Y, lambda), 1, 
                  function(z) solve(t(X) %*% X + z[n+1] * diag(1, ncol(X))) %*% t(X) %*% z[1:n])
    
  }else{
    
    if (ncol(Y) == 1) { # if y is a vector
      betas = solve(t(X) %*% X + lambda * diag(1, ncol(X))) %*% t(X) %*% Y  
    }else{
      betas = solve(t(X) %*% X + lambda * diag(1, ncol(X))) %*% t(X) %*% t(Y)  
    }
    
  }

  betas = t(betas)
  
  colnames(betas) = colnames(X)
  
  list(betas=betas, design=X)
}
