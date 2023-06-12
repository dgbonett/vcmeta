# ================= Miscellaneous Functions ================
# ci.fisher
#' Fisher confidence interval for any type of correlation.
#' 
#'
#' @description 
#' This function computes a confidence interval for any type
#' correlation using an estimated correlation and its standard 
#' error. This function should be used with the \link[vcmeta]{meta.ave.gen}
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


#  cor.from.t =============================================================
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


# meta.chitest ========================================================
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


#  stdmean2.from.t ============================================================
#' Computes Cohen's d from pooled-variance t statistic
#' 
#' @description 
#' This function computes Cohen's d for a 2-group design (which is a 
#' standardized mean difference with a weighted variance standardizer) using 
#' a pooled-variance independent-samples t statistic and the two sample sizes. 
#' This function also computes an equal-variance standard error for Cohen's d. 
#' 
#' 
#' @param t  	pooled-variance t statistic  
#' @param n1 	sample size for group 1 
#' @param n2 	sample size for group 2 
#'  
#' @return
#' Returns Cohen's d and its equal-variance standard error
#' 
#' @examples
#' stdmean2.from.t(3.27, 25, 25)
#'
#' # Should return:
#' #       Estimate       SE
#' # [1,] 0.9439677 0.298801
#' 
#' @export
stdmean2.from.t <- function(t, n1, n2) {
  d <- t*sqrt((n1 + n2)^2/(n1*n2*(n1 + n2 - 2)))
  se <- sqrt(d^2*(1/(n1 - 1) + 1/(n2 - 1))/8 + 1/n1 + 1/n2)
  out <- t(c(d, se))
  colnames(out) <- c("Estimate", "SE")
  return (out)
}

