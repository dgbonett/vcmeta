# ================= Miscellaneous Functions ============

# ci.fisher
#' Fisher confidence interval for any type of correlation.
#' 
#'
#' @description 
#' This function computes a confidence interval for any type
#' correlation using an estimated correlation and its standard 
#' error. This function should be used with the meta.ave.gen
#' function when the effect size is a correlation. Use the 
#' estimated average correlation and it standard error from 
#' meta.ave.gen (when the effect size is a correlation) in the 
#' ci.fisher function to obtain a more accurate confidence interval
#' for the population average correlation. 
#'
#'
#' @param alpha alpha value for 1-alpha confidence
#' @param cor   estimate of correlation 
#' @param se    standard error of estimated correlation
#' 
#' 
#' @return A 2-element vector with lower and upper bounds of the confidence
#' interval
#' 
#' @examples 
#' ci.fisher(0.05, 0.50, .10)
#'
#' # Should return:
#' # [1] 0.2802723 0.6699402
#' 
#' 
#' @importFrom Rdpack reprompt
#' @importFrom stats qnorm
#' @export
ci.fisher <- function(alpha, cor, se) {
  z <- qnorm(1 - alpha/2)
  zr <- log((1 + cor)/(1 - cor))/2
  ll0 <- zr - z*se/(1 - cor^2)
  ul0 <- zr + z*se/(1 - cor^2)
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- c(ll, ul)
  return(out)
}


#' Compute Pearson correlation between paired measurements from t statistic
#' 
#' @description 
#' This function computes the Pearson correlation between paired 
#' measurements using a reported paired-samples t statistic and 
#' other sample information. This correlation estimate is needed 
#' in several functions that analyze mean differences and 
#' standardized mean differences in paired-samples studies.
#' 
#' 
#' @param m1 sample mean for measurement 1 
#' @param m2 sample mean for measurement 2 
#' @param sd1 sample standard deviation for measurement 1 
#' @param sd2 sample standard deviation for measurement 2 
#' @param t value for paired-samples t-test
#' @param n sample size
#'  
#' @return
#' Returns the sample Pearson correlation between the two paired measurements
#' 
#' @examples
#' cor.from.t(9.4, 9.8, 1.26, 1.40, 2.27, 30)
#'
#' # Should return:
#' # [1] 0.7415209
#' 
#' @export
cor.from.t <- function(m1, m2, sd1, sd2, t, n) {
  cor <- ((sd1^2 + sd2^2) - n*(m1 - m2)^2/t^2)/(2*sd1*sd2)
  return (cor)
}
