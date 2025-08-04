#  meta.lm.mean2 ===========================================================
#' Meta-regression analysis for 2-group mean differences
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 2-group
#' mean difference. The estimates are OLS estimates with robust standard
#' errors that accommodate residual heteroscedasticity.  
#' 
#'  
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     m1    	vector of estimated means for group 1 
#' @param     m2    	vector of estimated means for group 2
#' @param     sd1   	vector of estimated SDs for group 1
#' @param     sd2   	vector of estimated SDs for group 2
#' @param     n1    	vector of group 1 sample sizes
#' @param     n2    	vector of group 2 sample sizes
#' @param     X     	matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix. The first row is for the intercept with one additional 
#' row per predictor. The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * t - t-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom 
#' 
#' 
#' @examples
#' n1 <- c(65, 30, 29, 45, 50)
#' n2 <- c(67, 32, 31, 20, 52)
#' m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
#' m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
#' sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
#' sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
#' x1 <- c(4, 6, 7, 7, 8)
#' x2 <- c(1, 0, 0, 0, 1)
#' X <- matrix(cbind(x1, x2), 5, 2)
#' meta.lm.mean2(.05, m1, m2, sd1, sd2, n1, n2, X)
#' 
#' # Should return:
#' #    Estimate        SE         t     p         LL        UL  df
#' # b0   -15.20 3.4097610 -4.457791 0.000 -21.902415 -8.497585 418
#' # b1     2.35 0.4821523  4.873979 0.000   1.402255  3.297745 418
#' # b2     2.85 1.5358109  1.855697 0.064  -0.168875  5.868875 418
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
meta.lm.mean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, X) {
  m <- length(m1)
  nt <- sum(n1 + n2)
  var <- sd1^2/n1 + sd2^2/n2
  d <- m1 - m2
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%d
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  df <- nt - q
  crit <- qt(1 - alpha/2, df)
  ll <- b - crit*se
  ul <- b + crit*se
  t <- b/se
  p <- round(2*(1 - pt(abs(t), df)), digits = 3)
  out <- cbind(b, se, t, p, ll, ul, df)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  rownames(out) <- row
  return(out)
}


#  meta.lm.stdmean2 ==========================================================
#' Meta-regression analysis for 2-group standardized mean differences
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 2-group
#' standardized mean difference. The estimates are OLS estimates with
#' robust standard errors that accommodate residual heteroscedasticity. 
#' Use the unweighted variance standardizer for 2-group experimental
#' designs, and use the weighted variance standardizer for 2-group
#' nonexperimental designs. A single-group standardizer can be used
#' in either experimental or nonexperimental designs.
#'  
#
#' @param     alpha 	alpha level for 1-alpha confidence
#' @param     m1    	vector of estimated means for group 1 
#' @param     m2    	vector of estimated means for group 2
#' @param     sd1   	vector of estimated SDs for group 1
#' @param     sd2   	vector of estimated SDs for group 2
#' @param     n1    	vector of group 1 sample sizes
#' @param     n2    	vector of group 2 sample sizes
#' @param     X     	matrix of predictor values
#' @param stdzr
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#' * set to 3 for square root weighted average variance standardizer
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n1 <- c(65, 30, 29, 45, 50)
#' n2 <- c(67, 32, 31, 20, 52)
#' m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
#' m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
#' sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
#' sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
#' x1 <- c(4, 6, 7, 7, 8)
#' X <- matrix(x1, 5, 1)
#' meta.lm.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, X, 0)
#' 
#' # Should return:
#' #      Estimate        SE         z p         LL         UL
#' # b0 -1.6988257 0.4108035 -4.135373 0 -2.5039857 -0.8936657
#' # b1  0.2871641 0.0649815  4.419167 0  0.1598027  0.4145255
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
meta.lm.stdmean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, X, stdzr) {
  df1 <- n1 - 1
  df2 <- n2 - 1
  m <- length(m1)
  nt <- sum(n1 + n2)
  z <- qnorm(1 - alpha/2)
  v1 <- sd1^2
  v2 <- sd2^2
  if (stdzr == 0) {
    s1 <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s1
    du <- (1 - 3/(4*n - 9))*d
    var <- d^2*(v1^2/df1 + v2^2/df2)/(8*s1^4) + (v1/df1 + v2/df2)/s1^2
  } else if (stdzr == 1) { 
    d <- (m1 - m2)/sd1
    du <- (1 - 3/(4*n1 - 5))*d
    var <- d^2/(2*df1) + 1/df1 + v2/(df2*v1)
  } else if (stdzr == 2) {
    d <- (m1 - m2)/sd2
    du <- (1 - 3/(4*n2 - 5))*d
    var <- d^2/(2*df2) + 1/df2 + v1/(df1*v2)
  } else {
    s2 <- sqrt((df1*v1 + df2*v2)/(df1 + df2))
    d <- (m1 - m2)/sd2
    du <- (1 - 3/(4*n - 9))*d
    var <- d^2*(1/df1 + 1/df2)/8 + (v1/n1 + v2/n2)/s2^2
  }
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%d
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return(out)
}


#  meta.lm.mean.ps ==========================================================
#' Meta-regression analysis for paired-samples mean differences
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a paired-samples
#' mean difference. The estimates are OLS estimates with robust standard
#' errors that accommodate residual heteroscedasticity. 
#' 
#'  
#' @param     alpha 	alpha level for 1-alpha confidence
#' @param     m1    	vector of estimated means for measurement 1 
#' @param     m2    	vector of estimated means for measurement 2
#' @param     sd1   	vector of estimated SDs for measurement 1
#' @param     sd2   	vector of estimated SDs for measurement 2
#' @param     cor   	vector of estimated correlations for paired measurements
#' @param     n     	vector of sample sizes
#' @param     X     	matrix of predictor values
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * t - t-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @examples
#' n <- c(65, 30, 29, 45, 50)
#' cor <- c(.87, .92, .85, .90, .88)
#' m1 <- c(20.1, 20.5, 19.3, 21.5, 19.4)
#' m2 <- c(10.4, 10.2, 8.5, 10.3, 7.8)
#' sd1 <- c(9.3, 9.9, 10.1, 10.5, 9.8)
#' sd2 <- c(7.8, 8.0, 8.4, 8.1, 8.7)
#' x1 <- c(2, 3, 3, 4, 4)
#' X <- matrix(x1, 5, 1)
#' meta.lm.mean.ps(.05, m1, m2, sd1, sd2, cor, n, X)
#' 
#' # Should return:
#' #    Estimate        SE        t     p        LL        UL  df
#' # b0     8.00 1.2491990 6.404104 0.000 5.5378833 10.462117 217
#' # b1     0.85 0.3796019 2.239188 0.026 0.1018213  1.598179 217
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
meta.lm.mean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, X) {
  m <- length(m1)
  nt <- sum(n)
  var <- (sd1^2 + sd2^2 - 2*cor*sd1*sd2)/n
  d <- m1 - m2
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%d
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  df <- nt - q
  t <- qt(1 - alpha/2, df)
  ll <- b - t*se
  ul <- b + t*se
  t <- b/se
  p <- round(2*(1 - pt(abs(t), df)), digits = 3)
  out <- cbind(b, se, t, p, ll, ul, df)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  rownames(out) <- row
  return(out)
}


#  meta.lm.stdmean.ps =======================================================
#' Meta-regression analysis for paired-samples standardized mean differences
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a paired-samples
#' standardized mean difference. The estimates are OLS estimates with  
#' robust standard errors that accommodate residual heteroscedasticity. 
#' 
#'  
#' @param     alpha 	alpha level for 1-alpha confidence
#' @param     m1    	vector of estimated means for measurement 1 
#' @param     m2    	vector of estimated means for measurement 2
#' @param     sd1   	vector of estimated SDs for measurement 1
#' @param     sd2   	vector of estimated SDs for measurement 2
#' @param     cor   	vector of estimated correlations for paired measurements
#' @param     n     	vector of sample sizes
#' @param     X     	matrix of predictor values
#' @param stdzr
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for measurement 1 SD standardizer 
#' * set to 2 for measurement 2 SD standardizer 
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' 
#' n <- c(65, 30, 29, 45, 50)
#' cor <- c(.87, .92, .85, .90, .88)
#' m1 <- c(20.1, 20.5, 19.3, 21.5, 19.4)
#' m2 <- c(10.4, 10.2, 8.5, 10.3, 7.8)
#' sd1 <- c(9.3, 9.9, 10.1, 10.5, 9.8)
#' sd2 <- c(7.8, 8.0, 8.4, 8.1, 8.7)
#' x1 <- c(2, 3, 3, 4, 4)
#' X <- matrix(x1, 5, 1)
#' meta.lm.stdmean.ps(.05, m1, m2, sd1, sd2, cor, n, X, 0)
#' 
#' # Should return:
#' #      Estimate         SE         z     p         LL        UL
#' # b0 1.01740253 0.25361725 4.0115667 0.000  0.5203218 1.5144832
#' # b1 0.04977943 0.07755455 0.6418635 0.521 -0.1022247 0.2017836
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.stdmean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, X, stdzr) {
  m <- length(m1)
  df <- n - 1
  z <- qnorm(1 - alpha/2)
  v1 <- sd1^2
  v2 <- sd2^2
  vd <- v1 + v2 - 2*cor*sd1*sd2
  if (stdzr == 0) {
    s <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s
    du <- sqrt((n - 2)/df)*d
    var <- d^2*(v1^2 + v2^2 + 2*cor^2*v1*v2)/(8*df*s^4) + vd/(df*s^2)
  } else if (stdzr == 1) { 
    d <- (m1 - m2)/sd1
    du <- (1 - 3/(4*df - 1))*d
    var <- d^2/(2*df) + vd/(df*v1)
  } else {
    d <- (m1 - m2)/sd2
    du <- (1 - 3/(4*df - 1))*d
    var <- d^2/(2*df) + vd/(df*v2)
  }
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%d
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return(out)
}


#  meta.lm.meanratio2 =======================================================
#' Meta-regression analysis for 2-group log mean ratios
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 2-group
#' log mean ratio. The estimates are OLS estimates with robust standard
#' errors that accommodate residual heteroscedasticity. The exponentiated 
#' slope estimate for a predictor variable describes a multiplicative
#' change in the mean ratio associated with a 1-unit increase in that 
#' predictor variable, controlling for all other predictor variables
#' in the model.
#' 
#'  
#' @param     alpha 	alpha level for 1-alpha confidence
#' @param     m1    	vector of estimated means for group 1 
#' @param     m2    	vector of estimated means for group 2
#' @param     sd1   	vector of estimated SDs for group 1
#' @param     sd2   	vector of estimated SDs for group 2
#' @param     n1    	vector of group 1 sample sizes
#' @param     n2    	vector of group 2 sample sizes
#' @param     X     	matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * exp(Estimate) - the exponentiated estimate
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' 
#' 
#' @references
#' \insertRef{Bonett2020}{vcmeta}
#'
#'
#' @examples
#' n1 <- c(65, 30, 29, 45, 50)
#' n2 <- c(67, 32, 31, 20, 52)
#' m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
#' m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
#' sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
#' sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
#' x1 <- c(4, 6, 7, 7, 8)
#' X <- matrix(x1, 5, 1)
#' meta.lm.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, X)
#' 
#' # Should return:
#' #       Estimate         SE          LL          UL         z p
#' # b0 -0.40208954 0.09321976 -0.58479692 -0.21938216 -4.313351 0
#' # b1  0.06831545 0.01484125  0.03922712  0.09740377  4.603078 0
#' #    exp(Estimate)  exp(LL)   exp(UL)
#' # b0     0.6689208 0.557219 0.8030148
#' # b1     1.0707030 1.040007 1.1023054
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.meanratio2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, X) {
  m <- length(m1)
  var1 <- sd1^2/(n1*m1^2) 
  var2 <- sd2^2/(n2*m2^2)
  var <- var1 + var2
  y <- log(m1/m2)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%y
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  crit <- qnorm(1 - alpha/2)
  ll <- b - crit*se
  ul <- b + crit*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul, exp(b), exp(ll), exp(ul))
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL", 
    "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- row
  return(out)
}


#  meta.lm.meanratio.ps =====================================================
#' Meta-regression analysis for paired-samples log mean ratios
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a paired-samples
#' log mean ratio. The estimates are OLS estimates with robust standard
#' errors that accommodate residual heteroscedasticity. The exponentiated 
#' slope estimate for a predictor variable describes a multiplicative
#' change in the mean ratio associated with a 1-unit increase in that 
#' predictor variable, controlling for all other predictor variables
#' in the model.
#' 
#' 
#' @param     alpha 	alpha level for 1-alpha confidence
#' @param     m1    	vector of estimated means for measurement 1 
#' @param     m2    	vector of estimated means for measurement 2
#' @param     sd1   	vector of estimated SDs for measurement 1
#' @param     sd2   	vector of estimated SDs for measurement 2
#' @param     cor   	vector of estimated correlations for paired measurements
#' @param     n     	vector of sample sizes
#' @param     X     	matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'  * exp(Estimate) - the exponentiated estimate 
#'  * exp(LL) - lower limit of the exponentiated confidence interval
#'  * exp(UL) - upper limit of the exponentiated confidence interval
#' 
#' 
#' @references
#' \insertRef{Bonett2020}{vcmeta}
#' 
#' 
#' @examples
#' n <- c(65, 30, 29, 45, 50)
#' cor <- c(.87, .92, .85, .90, .88)
#' m1 <- c(20.1, 20.5, 19.3, 21.5, 19.4)
#' m2 <- c(10.4, 10.2, 8.5, 10.3, 7.8)
#' sd1 <- c(9.3, 9.9, 10.1, 10.5, 9.8)
#' sd2 <- c(7.8, 8.0, 8.4, 8.1, 8.7)
#' x1 <- c(2, 3, 3, 4, 4)
#' X <- matrix(x1, 5, 1)
#' meta.lm.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, X)
#' 
#' # Should return: 
#' #      Estimate         SE           LL        UL        z     p
#' # b0 0.50957008 0.13000068  0.254773424 0.7643667 3.919749 0.000
#' # b1 0.07976238 0.04133414 -0.001251047 0.1607758 1.929697 0.054
#' #     exp(Estimate)   exp(LL)  exp(UL)
#' # b0       1.664575 1.2901693 2.147634
#' # b1       1.083030 0.9987497 1.174422
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.meanratio.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, X) {
  m <- length(m1)
  var <- (sd1^2/m1^2 + sd2^2/m2^2 - 2*cor*sd1*sd2/(m1*m2))/n
  y <- log(m1/m2)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%y
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  crit <- qnorm(1 - alpha/2)
  ll <- b - crit*se
  ul <- b + crit*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul, exp(b), exp(ll), exp(ul))
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL", 
    "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- row
  return(out)
}


#  meta.lm.cor.gen ==========================================================
#' Meta-regression analysis for correlations 
#' 
#'  
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 
#' Fisher-transformed correlation. The correlations can be of different types 
#' (e.g., Pearson, partial, Spearman). The estimates are OLS estimates
#' with robust standard errors that accommodate residual heteroscedasticity. 
#' This function uses estimated correlations and their standard errors as 
#' input. The correlations are Fisher-transformed and hence the parameter
#' estimates do not have a simple interpretation. However, the hypothesis 
#' test results can be used to decide if a population slope is either 
#' positive or negative.
#' 
#' 
#' @param     alpha	 alpha level for 1-alpha confidence
#' @param     cor		 vector of estimated correlations 
#' @param     se		 number of control variables
#' @param     X		   matrix of predictor values
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' 
#' cor <- c(.40, .65, .60, .45)
#' se <- c(.182, .114, .098, .132)
#' x1 <- c(18, 25, 23, 19)
#' X <- matrix(x1, 4, 1)
#' meta.lm.cor.gen(.05, cor, se, X)
#' 
#' # Should return: 
#' #       Estimate         SE          z     p
#' # b0 -0.47832153 0.63427931 -0.7541181 0.451
#' # b1  0.05047154 0.02879859  1.7525699 0.080
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.cor.gen <- function(alpha, cor, se, X) {
  m <- length(cor)
  z <- qnorm(1 - alpha/2)
  zcor <- log((1 + cor)/(1 - cor))/2
  zvar <- se^2/(1 - cor^2)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%zcor
  V <- diag(zvar)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p")
  rownames(out) <- row
  return (out)
}

#  meta.lm.cor ==============================================================
#' Meta-regression analysis for Pearson or partial correlations
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 
#' Fisher-transformed Pearson or partial correlation. The estimates are OLS
#' estimates with robust standard errors that accommodate residual heteroscedasticity.
#' The correlations are Fisher-transformed and hence the parameter estimates
#' do not have a simple interpretation. However, the hypothesis test results
#' can be used to decide if a population slope is either positive or negative.
#' 
#' 
#' @param     alpha	alpha level for 1-alpha confidence
#' @param     n     vector of sample sizes
#' @param     cor		vector of estimated Pearson or partial correlations 
#' @param     s     number of control variables
#' @param     X		  matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - Standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' 
#' n <- c(55, 190, 65, 35)
#' cor <- c(.40, .65, .60, .45)
#' q <- 0
#' x1 <- c(18, 25, 23, 19)
#' X <- matrix(x1, 4, 1)
#' meta.lm.cor(.05, n, cor, q, X)
#' 
#' # Should return: 
#' #       Estimate         SE         z     p           LL         UL
#' # b0 -0.47832153 0.48631509 -0.983563 0.325 -1.431481595 0.47483852
#' # b1  0.05047154 0.02128496  2.371231 0.018  0.008753794 0.09218929
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.cor <- function(alpha, n, cor, s, X) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  zcor <- log((1 + cor)/(1 - cor))/2
  zvar <- 1/(n - 3 - s)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%zcor
  V <- diag(zvar)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.lm.spear ============================================================
#' Meta-regression analysis for Spearman correlations 
#' 
#'  
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 
#' Fisher-transformed Spearman correlation. The estimates are OLS estimates 
#' with robust standard errors that accommodate residual heteroscedasticity. 
#' The correlations are Fisher-transformed and hence the parameter 
#' estimates do not have a simple interpretation. However, the hypothesis
#' test results can be used to decide if a population slope is either 
#' positive or negative.
#' 
#' 
#' @param     alpha	 alpha level for 1-alpha confidence
#' @param     n      vector of sample sizes
#' @param     cor		 vector of estimated Spearman correlations 
#' @param     X		   matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' 
#' n <- c(150, 200, 300, 200, 350)
#' cor <- c(.14, .29, .16, .21, .23)
#' x1 <- c(18, 25, 23, 19, 24)
#' X <- matrix(x1, 5, 1)
#' meta.lm.spear(.05, n, cor, X)
#' 
#' # Should return: 
#' #       Estimate         SE          z     p           LL         UL
#' # b0 -0.08920088 0.26686388 -0.3342561 0.738 -0.612244475 0.43384271
#' # b1  0.01370866 0.01190212  1.1517825 0.249 -0.009619077 0.03703639
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.spear <- function(alpha, n, cor, X) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  zcor <- log((1 + cor)/(1 - cor))/2
  zvar <- (1 + cor^2/2)/(n - 3)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%zcor
  V <- diag(zvar)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return (out)
}

#  meta.lm.semipart =========================================================
#' Meta-regression analysis for semipartial correlations
#' 
#'  
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a Fisher-transformed
#' semipartial correlation. The estimates are OLS estimates with robust
#' standard errors that accommodate residual heteroscedasticity.  The
#' correlations are Fisher-transformed and hence the parameter estimates
#' do not have a simple interpretation. However, the hypothesis test results
#' can be used to decide if a population slope is either positive or negative.
#'  
#' 
#' @param     alpha	 alpha level for 1-alpha confidence
#' @param     n      vector of sample sizes
#' @param     cor		 vector of estimated semipartial correlations
#' @param     r2   	 vector of squared multiple correlations for a model that
#' includes the IV and all control variables
#' @param     X		   matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' 
#' n <- c(128, 97, 210, 217)
#' cor <- c(.35, .41, .44, .39)
#' r2 <- c(.29, .33, .36, .39)
#' x1 <- c(18, 25, 23, 19)
#' X <- matrix(x1, 4, 1)
#' meta.lm.semipart(.05, n, cor, r2, X)
#' 
#' # Should return: 
#' #      Estimate        SE         z     p          LL         UL
#' # b0 0.19695988 0.3061757 0.6432905 0.520 -0.40313339 0.79705315
#' # b1 0.01055584 0.0145696 0.7245114 0.469 -0.01800004 0.03911172
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.semipart <- function(alpha, n, cor, r2, X) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  r0 <- r2 - cor^2
  zcor <- log((1 + cor)/(1 - cor))/2
  zvar = (r2^2 - 2*r2 + r0 - r0^2 + 1)/((1 - cor^2)^2*(n - 3))
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%zcor
  V <- diag(zvar)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return (out)
}

#  meta.lm.cronbach =========================================================
#' Meta-regression analysis for Cronbach reliabilities 
#' 
#'
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a log-complement
#' Cronbach reliablity. The estimates are OLS estimates with robust standard 
#' errors that accommodate residual heteroscedasticity. The exponentiated slope 
#' estimate for a predictor variable describes a multiplicative change in 
#' non-reliability associated with a 1-unit increase in that predictor 
#' variable, controlling for all other predictor variables in the model.
#' 
#' 
#' @param     alpha	 alpha level for 1-alpha confidence
#' @param     n		   vector of sample sizes 
#' @param     rel		 vector of estimated reliabilities 
#' @param     r		   number of measurements (e.g., items) 
#' @param     X		   matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - exponentiated OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the exponentiated confidence interval
#' * UL - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' n <- c(583, 470, 546, 680)
#' rel <- c(.91, .89, .90, .89)
#' x1 <- c(1, 0, 0, 0)
#' X <- matrix(x1, 4, 1)
#' meta.lm.cronbach(.05, n, rel, 10, X)
#' 
#' # Should return:
#' #      Estimate         SE          z     p         LL          UL
#' # b0 -2.2408328 0.03675883 -60.960391 0.000 -2.3128788 -2.16878684
#' # b1 -0.1689006 0.07204625  -2.344336 0.019 -0.3101087 -0.02769259
#' 
#' 
#' @references
#' \insertRef{Bonett2010}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.cronbach <- function(alpha, n, rel, r, X) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  hn <- m/sum(1/n)
  a <- ((r - 2)*(m - 1))^.25
  log.rel <- log(1 - rel) - log(hn/(hn - 1))
  var.rel <- 2*r*(1 - rel)^2/((r - 1)*(n - 2 - a))
  var.log <- var.rel/(1 - rel)^2
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%log.rel
  V <- diag(var.log)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.lm.oddsratio ========================================================
#' Meta-regression analysis for odds ratios 
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a log odds
#' ratio. The estimates are OLS estimates with robust standard errors
#' that accommodate residual heteroscedasticity. The exponentiated 
#' slope estimate for a predictor variable describes a multiplicative
#' change in the odds ratio associated with a 1-unit increase in that 
#' predictor variable, controlling for all other predictor variables
#' in the model.
#' 
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f1     	vector of group 1 frequency counts
#' @param     f2     	vector of group 2 frequency counts
#' @param     n1     	vector of group 1 sample sizes 
#' @param     n2     	vector of group 2 sample sizes
#' @param     X     	matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * exp(Estimate) - the exponentiated estimate 
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' n1 <- c(204, 201, 932, 130, 77)
#' n2 <- c(106, 103, 415, 132, 83)
#' f1 <- c(24, 40, 93, 14, 5)
#' f2 <- c(12, 9, 28, 3, 1)
#' x1 <- c(4, 4, 5, 3, 26)
#' x2 <- c(1, 1, 1, 0, 0)
#' X <- matrix(cbind(x1, x2), 5, 2)
#' meta.lm.oddsratio(.05, f1, f2, n1, n2, X)
#' 
#' # Should return:
#' #        Estimate         SE           z     p         LL         UL
#' # b0  1.541895013 0.69815801  2.20851868 0.027  0.1735305 2.91025958
#' # b1 -0.004417932 0.04840623 -0.09126784 0.927 -0.0992924 0.09045653
#' # b2 -1.071122269 0.60582695 -1.76803337 0.077 -2.2585213 0.11627674
#' #    exp(Estimate)   exp(LL)   exp(UL)
#' # b0     4.6734381 1.1894969 18.361564
#' # b1     0.9955918 0.9054779  1.094674
#' # b2     0.3426238 0.1045049  1.123307
#' 
#' 
#' @references
#' \insertRef{Bonett2015}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.oddsratio <- function(alpha, f1, f2, n1, n2, X) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  lor <- log((f1 + .5)*(n2 - f2 + .5)/((f2 + .5)*(n1 - f1 + .5)))
  var <- 1/(f1 + .5) + 1/(f2 + .5) + 1/(n1 - f1 + .5) + 1/(n2 - f2 + .5)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%lor
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  exp.b <- exp(b)
  exp.ll <- exp(ll)
  exp.ul <- exp(ul)
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul, exp.b, exp.ll, exp.ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL", 
    "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- row
  return (out)
}


#  meta.lm.propratio2 =======================================================
#' Meta-regression analysis for 2-group proportion ratios 
#' 
#'
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 2-group log
#' proportion ratio. The estimates are OLS estimates with robust standard 
#' errors that accommodate residual heteroscedasticity. The exponentiated 
#' slope estimate for a predictor variable describes a multiplicative 
#' change in the proportion ratio associated with a 1-unit increase in 
#' that predictor variable, controlling for all other predictor variables
#' in the model.
#' 
#' 
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f1     	vector of group 1 frequency counts
#' @param     f2     	vector of group 2 frequency counts
#' @param     n1     	vector of group 1 sample sizes 
#' @param     n2     	vector of group 2 sample sizes
#' @param     X      	matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * exp(Estimate) - the exponentiated estimate 
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' n1 <- c(204, 201, 932, 130, 77)
#' n2 <- c(106, 103, 415, 132, 83)
#' f1 <- c(24, 40, 93, 14, 5)
#' f2 <- c(12, 9, 28, 3, 1)
#' x1 <- c(4, 4, 5, 3, 26)
#' x2 <- c(1, 1, 1, 0, 0)
#' X <- matrix(cbind(x1, x2), 5, 2)
#' meta.lm.propratio2(.05, f1, f2, n1, n2, X)
#' 
#' # Should return:
#' #         Estimate         SE           z     p          LL         UL
#' # b0  1.4924887636 0.69172794  2.15762393 0.031  0.13672691 2.84825062
#' # b1  0.0005759509 0.04999884  0.01151928 0.991 -0.09741998 0.09857188
#' # b2 -1.0837844594 0.59448206 -1.82307345 0.068 -2.24894789 0.08137897
#' #     exp(Estimate)   exp(LL)   exp(UL)
#' # b0      4.4481522 1.1465150 17.257565
#' # b1      1.0005761 0.9071749  1.103594
#' # b2      0.3383128 0.1055102  1.084782
#' 
#' 
#' @references
#' \insertRef{Price2008}{vcmeta}
#'
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.propratio2 <- function(alpha, f1, f2, n1, n2, X) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  p1 <- (f1 + 1/4)/(n1 + 7/4) 
  p2 <- (f2 + 1/4)/(n2 + 7/4)
  lrr <- log(p1/p2)
  v1 <- 1/(f1 + 1/4 + (f1 + 1/4)^2/(n1 - f1 + 3/2))
  v2 <- 1/(f2 + 1/4 + (f2 + 1/4)^2/(n2 - f2 + 3/2))
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%lrr
  V <- diag(v1 + v2)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  exp.b <- exp(b)
  exp.ll <- exp(ll)
  exp.ul <- exp(ul)
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul, exp.b, exp.ll, exp.ul)
  row <- t(t(paste0(rep("b", q), seq(1:q)-1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL", 
    "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- row
  return (out)
}


#  meta.lm.prop2 ============================================================
#' Meta-regression analysis for 2-group proportion differences 
#' 
#'
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a 2-group
#' proportion difference. The estimates are OLS estimates with
#' robust standard errors that accommodate residual heteroscedasticity. 
#' 
#' 
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f1     	vector of group 1 frequency counts
#' @param     f2     	vector of group 2 frequency counts
#' @param     n1     	vector of group 1 sample sizes 
#' @param     n2     	vector of group 2 sample sizes
#' @param     X      	matrix of predictor values
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' f1 <- c(24, 40, 93, 14, 5)
#' f2 <- c(12, 9, 28, 3, 1)
#' n1 <- c(204, 201, 932, 130, 77)
#' n2 <- c(106, 103, 415, 132, 83)
#' x1 <- c(4, 4, 5, 3, 26)
#' x2 <- c(1, 1, 1, 0, 0)
#' X <- matrix(cbind(x1, x2), 5, 2)
#' meta.lm.prop2(.05, f1, f2, n1, n2, X)
#' 
#' # Should return:
#' #        Estimate          SE          z     p          LL          UL
#' # b0  0.089756283 0.034538077  2.5987632 0.009  0.02206290 0.157449671
#' # b1 -0.001447968 0.001893097 -0.7648672 0.444 -0.00515837 0.002262434
#' # b2 -0.034670988 0.034125708 -1.0159786 0.310 -0.10155615 0.032214170
#' 
#' 
#' @references
#' \insertRef{Bonett2014}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.prop2 <- function(alpha, f1, f2, n1, n2, X) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  p1 <- (f1 + 1/m)/(n1 + 2/m) 
  p2 <- (f2 + 1/m)/(n2 + 2/m)
  rd <- p1 - p2
  v1 <- p1*(1 - p1)/(n1 + 2/m)
  v2 <- p2*(1 - p2)/(n2 + 2/m) 
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%rd
  V <- diag(v1 + v2)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.lm.prop.ps ==========================================================
#' Meta-regression analysis for paired-samples proportion differences 
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients 
#' in a meta-regression model where the dependent variable is a 
#' paired-samples proportion difference. The estimates are OLS 
#' estimates with robust standard errors that accommodate residual 
#' heteroscedasticity.  		
#' 				
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f11    	vector of frequency counts in cell 1,1
#' @param     f12    	vector of frequency counts in cell 1,2
#' @param     f21    	vector of frequency counts in cell 2,1
#' @param     f22    	vector of frequency counts in cell 2,2
#' @param     X      	matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' f11 <- c(40, 20, 25, 30)
#' f12 <- c(3, 2, 2, 1)
#' f21 <- c(7, 6, 8, 6)
#' f22 <- c(26, 25, 13, 25)
#' x1 <- c(1, 1, 4, 6)
#' x2 <- c(1, 1, 0, 0)
#' X <- matrix(cbind(x1, x2), 4, 2)
#' meta.lm.prop.ps(.05, f11, f12, f21, f22, X)
#' 
#' # Should return: 
#' #       Estimate         SE          z     p          LL         UL
#' # b0 -0.21113402 0.21119823 -0.9996960 0.317 -0.62507494 0.20280690
#' # b1  0.02185567 0.03861947  0.5659236 0.571 -0.05383711 0.09754845
#' # b2  0.12575138 0.17655623  0.7122455 0.476 -0.22029248 0.47179524
#' 
#' 
#' @references
#' \insertRef{Bonett2012}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.prop.ps <- function(alpha, f11, f12, f21, f22, X) {
  m <- length(f11)
  z <- qnorm(1 - alpha/2)
  n <- f11 + f12 + f21 + f22
  p12 <- (f12 + 1/m)/(n + 2/m) 
  p21 <- (f21 + 1/m)/(n + 2/m)
  rd <- p12 - p21
  var <- (p12 + p21 - rd^2)/(n + 2/m)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%rd
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.lm.agree ============================================================
#' Meta-regression analysis for G agreement indices
#' 
#'
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a G-index of
#' agreement. The estimates are OLS estimates with robust standard errors 
#' that accomodate residual heteroscedasticity. 
#' 
#' 
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f11    	vector of frequency counts in cell 1,1
#' @param     f12    	vector of frequency counts in cell 1,2
#' @param     f21    	vector of frequency counts in cell 2,1
#' @param     f22    	vector of frequency counts in cell 2,2
#' @param     X      	matrix of predictor values
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval 
#' 
#' 
#' @examples
#' f11 <- c(40, 20, 25, 30)
#' f12 <- c(3, 2, 2, 1)
#' f21 <- c(7, 6, 8, 6)
#' f22 <- c(26, 25, 13, 25)
#' x1 <- c(1, 1, 4, 6)
#' x2 <- c(1, 1, 0, 0)
#' X <- matrix(cbind(x1, x2), 4, 2)
#' meta.lm.agree(.05, f11, f12, f21, f22, X)
#' 
#' # Should return:
#' #     Estimate         SE         z     p          LL        UL
#' # b0 0.1904762 0.38772858 0.4912617 0.623 -0.56945786 0.9504102
#' # b1 0.0952381 0.07141957 1.3335013 0.182 -0.04474169 0.2352179
#' # b2 0.4205147 0.32383556 1.2985438 0.194 -0.21419136 1.0552207
#' 
#' 
#' @references 
#' \insertRef{Bonett2022}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.agree <- function(alpha, f11, f12, f21, f22, X) {
  m <- length(f11)
  z <- qnorm(1 - alpha/2)
  n <- f11 + f12 + f21 + f22
  p0 <- (f11 + f22 + 2/m)/(n + 4/m)
  g <- 2*p0 - 1 
  var <- 4*p0*(1 - p0)/(n + 4/m)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%g
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.lm.mean1 ============================================================
#' Meta-regression analysis for 1-group means
#' 
#'
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a mean
#' from one group. The estimates are OLS estimates with robust
#' standard errors that accomodate residual heteroscedasticity. 
#' 
#' 
#' @param     alpha 	alpha level for 1-alpha confidence
#' @param     m     	vector of estimated means
#' @param     sd    	vector of estimated standard deviations
#' @param     n     	vector of sample sizes
#' @param     X     	matrix of predictor values
#' 
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * t - t-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @examples
#' n <- c(25, 15, 30, 25, 40)
#' m <- c(20.1, 20.5, 19.3, 21.5, 19.4)
#' sd <- c(10.4, 10.2, 8.5, 10.3, 7.8)
#' x1 <- c(1, 1, 0, 0, 0)
#' x2 <- c( 12, 13, 11, 13, 15)
#' X <- matrix(cbind(x1, x2), 5, 2)
#' meta.lm.mean1(.05, m, sd, n, X)
#' 
#' # Should return: 
#' #       Estimate        SE          t     p         LL        UL  df
#' # b0 19.45490196 6.7873381 2.86635227 0.005  6.0288763 32.880928 132
#' # b1  0.25686275 1.9834765 0.12950128 0.897 -3.6666499  4.180375 132
#' # b2  0.04705882 0.5064693 0.09291544 0.926 -0.9547876  1.048905 132
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
meta.lm.mean1 <- function(alpha, m, sd, n, X) {
  k <- length(m)
  var <- sd^2/n
  x0 <- matrix(c(1), k, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%m
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  df <- sum(n) - q
  t <- qt(1 - alpha/2, df)
  ll <- b - t*se
  ul <- b + t*se
  t <- b/se
  p <- round(2*(1 - pt(abs(t), df)), digits = 3)
  out <- cbind(b, se, t, p, ll, ul, df)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  rownames(out) <- row
  return(out)
}


#  meta.lm.prop1 ============================================================
#' Meta-regression analysis for 1-group proportions
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is a proportion
#' from one group. The estimates are OLS estimates with robust
#' standard errors that accomodate residual heteroscedasticity. 
#' 
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f      	vector of frequency counts
#' @param     n      	vector of sample sizes
#' @param     X     	matrix of predictor values
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' f <- c(38, 26, 24, 15, 45, 38)
#' n <- c(80, 60, 70, 50, 180, 200)
#' x1 <- c(10, 15, 18, 22, 24, 30)
#' X <- matrix(x1, 6, 1)
#' meta.lm.prop1(.05, f, n, X)
#' 
#' # Should return: 
#' #       Estimate         SE         z p          LL           UL
#' # b0  0.63262816 0.06845707  9.241239 0  0.49845477  0.766801546
#' # b1 -0.01510565 0.00290210 -5.205076 0 -0.02079367 -0.009417641
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.prop1 <- function(alpha, f, n, X) {
  z <- qnorm(1 - alpha/2)
  k <- length(f)
  p <- (f + 2/k)/(n + 4/k)
  var <- p*(1 - p)/(n + 4/k)
  x0 <- matrix(c(1), k, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%p
  V <- diag(var)
  se <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*se
  ul <- b + z*se
  z <- b/se
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, se, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return(out)
}


#  meta.lm.gen ==============================================================
#' Meta-regression analysis for any type of effect size
#' 
#' 
#' @description
#' This function estimates the intercept and slope coefficients in a
#' meta-regression model where the dependent variable is any type of
#' effect size. The estimates are OLS estimates with robust standard 
#' errors that accomodate residual heteroscedasticity. 
#' 
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     est    	vector of parameter estimates 
#' @param     se     	vector of standard errors
#' @param     X      	matrix of predictor values
#' 
#' @return
#' Returns a matrix.  The first row is for the intercept with one additional 
#' row per predictor.  The matrix has the following columns:
#' * Estimate - OLS estimate
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#'   
#' @examples
#' est <- c(4.1, 4.7, 4.9, 5.7, 6.6, 7.3)
#' se <- c(1.2, 1.5, 1.3, 1.8, 2.0, 2.6)
#' x1 <- c(10, 20, 30, 40, 50, 60)
#' x2 <- c(1, 1, 1, 0, 0, 0)
#' X <- matrix(cbind(x1, x2), 6, 2)
#' meta.lm.gen(.05, est, se, X)
#' 
#' # Should return:
#' #      Estimate         SE           z     p         LL         UL
#' # b0  3.5333333 4.37468253  0.80767766 0.419 -5.0408869 12.1075535
#' # b1  0.0600000 0.09058835  0.66233679 0.508 -0.1175499  0.2375499
#' # b2 -0.1666667 2.81139793 -0.05928249 0.953 -5.6769054  5.3435720
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
meta.lm.gen <- function(alpha, est, se, X) {
  m <- length(est)
  z <- qnorm(1 - alpha/2)
  x0 <- matrix(c(1), m, 1)
  X <- cbind(x0, X)
  q <- ncol(X)
  M <- solve(t(X)%*%X)
  b <- M%*%t(X)%*%est
  V <- diag(se^2)
  seb <- sqrt(diag(M%*%t(X)%*%V%*%X%*%M))
  ll <- b - z*seb
  ul <- b + z*seb
  z <- b/seb
  p <- round(2*(1 - pnorm(abs(z))), digits = 3)
  out <- cbind(b, seb, z, p, ll, ul)
  row <- t(t(paste0(rep("b", q), seq(1:q) - 1)))
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- row
  return(out)
}

