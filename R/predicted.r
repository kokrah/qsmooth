#' Compute predicted alpha values
#' Notes:
#' 1. Takes in results of fitCoeffs()
#' # Output:
#' matrix of predicted values
#' @param design = design matrix
#' @param betas = matrix of betas
#' @export 
predicted = function (design, betas) {
  design %*% t(betas)
}
