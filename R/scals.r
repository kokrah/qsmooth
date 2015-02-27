#' Compute sample quantiles, the quantile reference, and alpha matrix.
#' Notes:
#' 1. 
#' Output:
#' A list: (1.) log2(sample quantiles) (2.) log2(quantile ref.) (3.) alphas
#' @param data (positive values: raw intensities or pseudo-counts)
#' @export 
scals <- function (data) {
  
  # Compute sample quantiles
  Q <- apply(data, 2, sort) 
  log2.Q <- log2(Q)
  
  # Compute quantile reference
  log2.Qref <- rowMeans(log2.Q)
  
  # Compute alpha matrix
  alpha <- log2.Q - log2.Qref 
  
  list(log2.Q = log2.Q,
       log2.Qref = log2.Qref,
       alpha = alpha)
}