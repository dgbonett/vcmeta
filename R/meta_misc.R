# ================= Miscellaneous Functions ============

# meta.ave.fisher =============================================================
#' Fisher confidence interval for an average correlation.
#' 
#'
#' @description 
#' This function should be used with the \link[vcmeta]{meta.ave.gen}
#' function when the effect size is a correlation. Use the estimated average
#' correlation and it standard error from meta.ave.gen in this function to 
#' obtain a more accurate confidence interval for the population average 
#' correlation. 
#'
#'
#' @param alpha alpha value for 1-alpha confidence
#' @param cor   estimate of average correlation 
#' @param se    standard error of average correlation
#' 
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of average correlation (from input) 
#' * LL - lower limit of the confidence interval
#' * UL - lower limit of the confidence interval
#' 
#' @examples 
#' meta.ave.fisher(0.05, 0.376, .054)
#'
#' # Should return:
#' # Estimate        LL        UL
#' #    0.376 0.2656039 0.4766632
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.fisher <- function(alpha, cor, se) {
  z <- qnorm(1 - alpha/2)
  zr <- log((1 + cor)/(1 - cor))/2
  ll0 <- zr - z*se/(1 - cor^2)
  ul0 <- zr + z*se/(1 - cor^2)
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- t(c(cor, ll, ul))
  colnames(out) <- c("Estimate", "LL", "UL")
  rownames(out) <- ""
  return(out)
}


#  cor.from.t
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
#' @param m1 estimated mean for measurement 1 
#' @param m2 estimated mean for measurement 2 
#' @param sd1 estimated standard deviation for measurement 1 
#' @param sd2 estimated standard deviation for measurement 2 
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


# meta.chitest 
#' Computes a chi-square test of effect-size homogeneity
#' 
#'
#' @description
#' Computes a chi-square test of effect size homogeneity and p-value using 
#' effect-size estimates and their standard errors from two or more studies.
#' This test should not be used to justify the use of a constant coeffient
#' (fixed-effect) meta-analysis. This test can be used to justify the 
#' estimation of an average effect size in a varying coefficient model.
#' 
#' 
#' @param    est  	vector of effect-size estimates
#' @param    se		vector of effect-size standard errors
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Q - chi-square test statitic
#' * df - degrees of freedom
#' * p - p-value
#' 
#'  
#' @examples
#' est <- c(.297, .324, .281, .149) 
#' se <- c(.082, .051, .047, .094)
#' meta.chitest(est, se)
#'
#' # Should return:
#' #         Q df         p
#' #  2.706526  3 0.4391195
#' 
#' 
#' @references
#'  \insertRef{Borenstein2009}{vcmeta}
#'
#'
#' @importFrom stats pchisq
#' @export
meta.chitest <- function(est, se) {
 df <- length(est) - 1
 w <- 1/se^2
 ave <- sum(w*est)/sum(w)
 Q <- sum(w*(est - ave)*(est - ave))
 p <- 1 - pchisq(Q, df)
 out <- t(c(Q, df, p))
 colnames(out) <- c("Q", "df", "p")
 rownames(out) <- c("")
 return(out)
} 


