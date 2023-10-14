# meta.ave.fisher ============================================================
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


#  cor.from.t =============================================================
#' Computes Pearson correlation between paired measurements from t statistic
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
#' #                Estimate
#' # Correlation:  0.7415209
#' 
#'
#' @export
cor.from.t <- function(m1, m2, sd1, sd2, t, n) {
  out <- t(((sd1^2 + sd2^2) - n*(m1 - m2)^2/t^2)/(2*sd1*sd2))
  colnames(out) <- c("Estimate")
  rownames(out) <- c("Correlation: ")
  return (out)
}


# meta.chitest ========================================================
#' Computes a chi-square test of effect-size homogeneity
#' 
#'
#' @description
#' Computes a chi-square test of effect size homogeneity and p-value using 
#' effect-size estimates and their standard errors from two or more studies.
#' This test should not be used to justify the use of a constant coeffient
#' (fixed-effect) meta-analysis. 
#' 
#' 
#' @param    est  	vector of effect-size estimates
#' @param    se		  vector of effect-size standard errors
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
 rownames(out) <- ""
 return(out)
} 


#  stdmean2.from.t ============================================================
#' Computes Cohen's d from pooled-variance t statistic
#' 
#'
#' @description 
#' This function computes Cohen's d for a 2-group design (which is a 
#' standardized mean difference with a weighted variance standardizer) using 
#' a pooled-variance independent-samples t statistic and the two sample sizes. 
#' This function also computes the standard error for Cohen's d. The Cohen's d
#' estimate and standard error assume equality of population variances. 
#' 
#' 
#' @param t  	   pooled-variance t statistic  
#' @param n1 	   sample size for group 1 
#' @param n2 	   sample size for group 2 
#'  
#' @return
#' Returns Cohen's d and its equal-variance standard error
#' 
#' @examples
#' stdmean2.from.t(3.27, 25, 25)
#'
#' # Should return:
#' #             Estimate       SE
#' # Cohen's d  0.9439677 0.298801
#' 
#' @export
stdmean2.from.t <- function(t, n1, n2) {
  d <- t*sqrt(1/n1 + 1/n2)
  se <- sqrt(d^2*(1/(n1 - 1) + 1/(n2 - 1))/8 + 1/n1 + 1/n2)
  out <- t(c(d, se))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Cohen's d: ")
  return (out)
}


#' table.from.odds ============================================================
#' Computes the cell frequencies in a 2x2 table using the marginal proportions
#' and odds ratio 
#'                   
#'                         
#' @description 
#' This function computes the cell proportions and frequencies in a 2x2 
#' contingency table using the reported marginal proportions, estimated odds 
#' ratio, and total sample size. The cell frequncies could then be used to
#' compute other measures of effect size. In the output, "cell ij" refers to
#' row i and column j.
#' 
#' 
#' @param    p1row	    marginal proportion for row 1
#' @param    p1col	    marginal proportion for column 1 
#' @param    or         estimated odds ratio
#' @param    n          total sample size
#' 
#' 
#' @return A 2-row matrix. The rows are:
#' * Row 1 gives the four computed cell proportions 
#' * Row 2 gives the four computed cell frequencies
#'
#'
#' The columns are:
#' * Cell 11 - proportion and frequency for cell 11
#' * Cell 12 - proportion and frequency for cell 12
#' * Cell 21 - proportion and frequency for cell 21
#' * Cell 22 - proportion and frequency for cell 22
#'    
#' 
#' @examples
#' table.from.odds(.17, .5, 3.18, 100)
#'
#' # Should return:
#' #                cell 11    cell 12    cell 21    cell 22
#' # Proportion:  0.1233262 0.04667383  0.3766738  0.4533262
#' # Frequency:  12.0000000 5.00000000 38.0000000 45.0000000
#' 
#' 
#' @references
#' \insertRef{Bonett2007}{vcmeta}           
#' 
#' 
#' @export
table.from.odds <- function(p1row, p1col, or, n){
 if (or <= 0) {stop("the odds ratio must be greater than 0")}
 p2row <- 1 - p1row
 if (or != 1){
  a <- or*(p1row + p1col) + p2row - p1col
  b <- sqrt(a^2 - 4*p1row*p1col*or*(or - 1))
  p11 <- (a - b)/(2*(or - 1))}
 else {
  p11 <- p1row*p1col
 }
 p12 <- p1row - p11
 p21 <- p1col - p11
 p22 <- 1 - (p11 + p12 + p21)
 f11 <- round(n*p11)
 f12 <- round(n*p12)
 f21 <- round(n*p21)
 f22 <- n - (f11 + f12 + f21)
 out1 <- t(c(p11, p12, p21, p22))
 out2 <- t(c(f11, f12, f21, f22))
 out <- rbind(out1, out2)
 colnames(out) <- c("cell 11", "cell 12", " cell 21", "cell 22")
 rownames(out) <- c("Proportion:", "Frequency:")
 return(out)
}


#' table.from.phi ============================================================
#' Computes the cell frequencies in a 2x2 table using the marginal proportions
#' and phi correlation 
#'                   
#'                         
#' @description 
#' This function computes the cell proportions and frequencies in a 2x2 
#' contingency table using the reported marginal proportions, estimated phi  
#' correlation, and total sample size. The cell frequncies could then be used  
#' to compute other measures of effect size. In the output, "cell ij" refers 
#' to row i and column j. 
#' 
#' 
#' @param    p1row	    marginal proportion for row 1
#' @param    p1col	    marginal proportion for column 1 
#' @param    phi        estimated phi correlation
#' @param    n          total sample size
#' 
#' 
#' @return A 2-row matrix. The rows are:
#' * Row 1 gives the four computed cell proportions 
#' * Row 2 gives the four computed cell frequencies
#'
#'
#' The columns are:
#' * Cell 11 - proportion and frequency for cell 11
#' * Cell 12 - proportion and frequency for cell 12
#' * Cell 21 - proportion and frequency for cell 21
#' * Cell 22 - proportion and frequency for cell 22
#'    
#' 
#' @examples
#' table.from.phi(.28, .64, .38, 200)
#'
#' # Should return:
#' #                cell 11   cell 12    cell 21    cell 22
#' # Proportion:  0.2610974 0.0189026  0.3789026  0.3410974
#' # Frequency   52.0000000 4.0000000 76.0000000 68.0000000
#' 
#' 
#' @export                                  
table.from.phi <- function(p1row, p1col, phi, n){
 if (abs(phi) > 1) {stop("phi must be between -1 and 1")}
 p2row <- 1 - p1row
 p2col <- 1 - p1col
 phimax <- sqrt(p1col*p2row/(p1row*p2col))
 if (phimax > 1) {phimax = 1/phimax}
 phimin <- sqrt(p2col*p2row/(p1row*p1col))
 if (phimin > 1) {phimin = 1/phimin}
 if (phi > phimax) {stop("phi is too large for given marginal proportions")}
 if (phi < -phimin) {stop("phi is too small for given marginal proportions")}
 a <- sqrt(p1row*p2row*p1col*p2col)
 p11 <- a*phi + p1row*p1col
 p12 <- p1row - p11
 p21 <- p1col - p11
 p22 <- 1 - (p11 + p12 + p21)
 f11 <- round(n*p11)
 f12 <- round(n*p12)
 f21 <- round(n*p21)
 f22 <- n - (f11 + f12 + f21)
 out1 <- t(c(p11, p12, p21, p22))
 out2 <- t(c(f11, f12, f21, f22))
 out <- rbind(out1, out2)
 colnames(out) <- c("cell 11", "cell 12", " cell 21", "cell 22")
 rownames(out) <- c("Proportion:", "Frequency")
 return(out)
}

