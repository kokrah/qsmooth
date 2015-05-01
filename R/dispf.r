#' Generate genewise dispersion 
#' 
#' @param mu base line counts
#' @param a asymptotic dispersion
#' @param chiDF degree of freedom gene specific disp. factor
dispf = function (mu, a, chiDF) {
  delta = chiDF / rchisq(length(mu), df=chiDF)
  delta * (a + (1 / mu))^2
}