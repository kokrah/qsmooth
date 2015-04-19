#' Scale normalization
#' 
#' @param exprs counts + 1 matrix or PM values
#' @param type type of scaling (default = "trim.mean")
#' @export 
scaleData = function (exprs, type="trim.mean") {
  
  if (type=="trim.mean") {
    res = t(t(exprs) / apply(exprs, 2, mean, trim=0.25))
  }
  
  if (type=="AH") {
    res = t(t(exprs) / matrixStats::colMedians(as.matrix(exprs / exp(rowMeans(log(exprs))))))
  }
  
  if (type=="median") {
    res = t(t(exprs) / matrixStats::colMedians(as.matrix(exprs)))
  }
  
  if (type=="mean") {
    res = t(t(exprs) / colMeans(exprs))
  }
  
  res = res / min(res)
  res
}