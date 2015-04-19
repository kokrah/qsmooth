#' Filter out low counts
#' 
#' @param counts counts matrix (do not add 1)
#' @param thresh minimum counts per million to determine expression (default=1)
#' @param minSamples minimum number of samples where gene is required to be expressed. 
#' This should be set to the numer of samples in the smallest group of interest. (default=2)
#' @export
filterCounts = function (counts, thresh = 1, minSamples = 2) {
  cpms = t(t(counts + 0.5) / (colSums(counts) + 1)) * 1e06
  keep = rowSums(cpms > thresh) >= minSamples
  counts = counts[keep, ]
  counts
}