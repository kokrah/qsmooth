#' Median scale data
#' Notes:
#' 1. 
#' Output:
#' Median scaled data
#' @param data (positive values: raw intensities or pseudo-counts)
#' @export
medScal = function (data) {
  meds = apply(data, 2, median, na.rm=TRUE)
  sweep(data, 2, meds, "/")
}