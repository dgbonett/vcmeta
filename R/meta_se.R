# ================= Effect Size Standard Errors ============

# se.mean2 
#' Computes the standard error for a 2-group mean difference
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' 2-group mean difference using the two sample means, sample standard
#' deviations, and sample sizes. The effect size estimate and standard 
#' error output from this function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' applications where compatible mean differences from a combination of 2-group
#' and paired-samples experiments are used in the meta-analysis. 
#' 
#' 
#' @param    m1		sample mean for group 1 
#' @param    m2		sample mean for group 2 
#' @param    sd1	sample standard deviation for group 1
#' @param    sd2	sample standard deviation for group 2
#' @param    n1		group 1 sample size
#' @param    n2		group 2 sample size
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of mean difference
#' * SE - standard error
#' 
#'  
#' @examples
#' se.mean2(21.9, 16.1, 3.82, 3.21, 40, 40)
#'
#  # Should return:
#' #                  Estimate        SE
#' # Mean difference:      5.8 0.7889312
#' 
#' 
#' @references
#' \insertRef{Snedecor1980}{vcmeta}
#'
#'
#' @export
se.mean2 <- function(m1, m2, sd1, sd2, n1, n2) {
  d <- m1 - m2
  se <- sqrt(sd1^2/n1 + sd2^2/n2)
  out <- t(c(d, se))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Mean difference:")
  return(out)
}


# se.mean.ps
#' Computes the standard error for a paired-samples mean difference
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' paired-samples mean difference using the two sample means, sample 
#  standard deviations, Pearson correlation, and sample size. The 
#' effect size estimate and standard error output from this function
#' can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' applications where compatible mean differences from a combination of 2-group
#' and paired-samples experiments are used in the meta-analysis. 
#' 
#' 
#' @param    m1		sample mean for measurement 1 
#' @param    m2		sample mean for measurement 2 
#' @param    sd1	sample standard deviation for measurement 1
#' @param    sd2	sample standard deviation for measurement 2
#' @param    cor	sample correlation for measurements 1 and 2
#' @param    n		sample size
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of mean difference
#' * SE - standard error
#' 
#' 
#' @examples
#' se.mean.ps(23.9, 25.1, 1.76, 2.01, .78, 25)
#'
#' # Should return:
#' #                  Estimate        SE
#' # Mean difference:     -1.2 0.2544833
#' 
#' @references
#' \insertRef{Snedecor1980}{vcmeta}
#' 
#' 
#' @export
se.mean.ps <- function(m1, m2, sd1, sd2, cor, n) {
  d <- m1 - m2
  se <- sqrt((sd1^2 + sd2^2 - 2*cor*sd1*sd2)/n)
  out <- t(c(d, se))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Mean difference:")
  return(out)
}


# se.stdmean2 		
#' Computes the standard error for a 2-group standardized mean difference
#' 
#' 
#' @description
#' Use the square root average variance standardizer (stdzr = 0) for 2-group
#' experimental designs. Use the square root weighted variance standardizer
#' (stdzr = 3) for 2-group nonexperimental designs with simple random sampling.
#' The single-group standardizers (stdzr = 1 and stdzr = 2) can be used with
#' either 2-group experimental or nonexperimental designs. The effect size 
#' estimate and standard error output from this function can be used as input
#' in the \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} 
#' functions in applications where compatible standardized mean differences from a 
#' combination of 2-group and paired-samples experiments are used in the meta-analysis. 
#'
#' 
#' @param    m1		  sample mean for group 1 
#' @param    m2		  sample mean for group 2 
#' @param    sd1		sample standard deviation for group 1
#' @param    sd2		sample standard deviation for group 2
#' @param    n1		  group 1 sample size
#' @param    n2		  group 2 sample size
#' @param stdzr
#' * set to 0 for square root average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#' * set to 3 for square root weighted variance standardizer
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of standardized mean difference
#' * SE - standard error
#' 
#' 
#' @examples
#' se.stdmean2(21.9, 16.1, 3.82, 3.21, 40, 40, 0)
#'
#' # Should return: 
#' #                               Estimate        SE
#' # Standardized mean difference: 1.643894 0.2629049
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#'
#'
#' @export
se.stdmean2 <- function(m1, m2, sd1, sd2, n1, n2, stdzr) {
  df1 <- n1 - 1
  df2 <- n2 - 1
  if (stdzr == 0) {
    s <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s
    se <- sqrt(d^2*(sd1^4/df1 + sd2^4/df2)/(8*s^4) + (sd1^2/df1 + sd2^2/df2)/s^2) 
  }
  else if (stdzr == 1) {
    s <- sd1
    d <- (m1 - m2)/s
    se <- sqrt(d^2/(2*df1) + 1/df1 + sd2^2/(df2*sd1^2))
  }
  else if (stdzr == 2) {
    s <- sd2
    d <- (m1 - m2)/s
    se <- sqrt(d^2/(2*df2) + 1/df2 + sd1^2/(df1*sd2^2))
  }
  else {
    s <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s
    se <- sqrt(d^2*(1/df1 + 1/df2)/8 + 1/n1 + 1/n2)
  }
  out <- t(c(d, se))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Standardized mean difference:")
  return(out)
}


# se.stdmean.ps 		
#' Computes the standard error for a paired-samples standardized mean difference
#' 
#'
#' @description
#' The effect size estimate and standard error output from this function can be used 
#' as input in the \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, and 
#' \link[vcmeta]{meta.lm.gen} functions in applications where compatible standardized 
#' mean differences from a combination of 2-group and paired-samples experiments are 
#' used in the meta-analysis. 
#' 
#' 
#' @param    m1		sample mean for measurement 1 
#' @param    m2		sample mean for measurement 2 
#' @param    sd1	sample standard deviation for measurement 1
#' @param    sd2	sample standard deviation for measurement 2
#' @param    cor	sample correlation for measurements 1 and 2
#' @param    n		sample size
#' @param stdzr
#' * set to 0 for square root average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#'  
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of standardized mean difference
#' * SE - standard error
#' 
#' 
#' @examples
#' se.stdmean.ps(23.9, 25.1, 1.76, 2.01, .78, 25, 0)
#'
#' # Should return: 
#' #                                  Estimate        SE
#' # Standardizedd mean difference: -0.6352097 0.1602852
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @export
se.stdmean.ps <- function(m1, m2, sd1, sd2, cor, n, stdzr) {
  df <- n - 1
  v1 <- sd1^2
  v2 <- sd2^2
  vd <- v1 + v2 - 2*cor*sd1*sd2
  if (stdzr == 0) {
    s <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s
    se <- sqrt(d^2*(v1^2 + v2^2 + 2*cor^2*v1*v2)/(8*df*s^4) + vd/(df*s^2))
  }
  else if (stdzr == 1) {
    s <- sd1
    d <- (m1 - m2)/s
    se <- sqrt(d^2/(2*df) + vd/(df*v1))
  }
  else {
    s <- sd2
    d <- (m1 - m2)/s
    se <- sqrt(d^2/(2*df) + vd/(df*v2))
  }
  out <- t(c(d, se))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Standardized mean difference:")
  return(out)
}


# se.cor
#' Computes the standard error for a Pearson or partial correlation 
#' 
#'
#' @description
#' This function can be used to compute the standard error of a 
#' Pearson or partial correlation using the sample correlation, sample 
#' size, and number of control variables. The effect size estimate and 
#' standard error output from this function can be used as input in the
#' \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, and
#' \link[vcmeta]{meta.lm.gen} functions in applications where a combination
#' of different types of correlations are used in the meta-analysis. 
#' 
#' 
#' @param    cor	sample Pearson or partial correlation  
#' @param    s		number of control variables (0 for Pearson)  
#' @param    n		sample size
#'   
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of Pearson or partial correlation
#' * SE - standard error
#' 
#' 
#' @examples
#' se.cor(.40, 0, 55)
#'
#' # Should return: 
#' #              Estimate       SE
#' # Correlation:      0.4 0.116487
#' 
#' 
#' @references
#' \insertRef{Bonett2008a}{vcmeta}
#' 
#' 
#' @export
se.cor <- function(cor, q, n) {
  se.cor <- sqrt((1 - cor^2)^2/(n - 3 - s))
  out <- t(c(cor, se.cor))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Correlation:")
  return(out)
}


# se.spear
#' Computes the standard error for a Spearman correlation 
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' Spearman correlation using the sample correlation and sample 
#' size. The standard error from this function can be used as input 
#' in the \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, and
#' \link[vcmeta]{meta.lm.gen} functions in applications where a combination
#' of different types of correlations are used in the meta-analysis. 
#' 
#' 
#' @param    cor		sample Spearman correlation  
#' @param    n		  sample size
#'   
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of Spearman correlation
#' * SE - standard error
#' 
#' 
#' @examples
#' se.spear(.40, 55)
#'
#' # Should return: 
#' #                       Estimate        SE
#' # Spearman correlation:      0.4 0.1210569
#' 
#' 
#' @references
#' \insertRef{Bonett2000}{vcmeta}
#'
#'
#' @export
se.spear <- function(cor, n) {
  se.cor <- sqrt((1 - cor^2)^2*(1 + cor^2/2)/(n - 3))
  out <- t(c(cor, se.cor))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Correlation:")
  return(out)
}


# se.semipart
#' Computes the standard error for a semipartial correlation 
#' 
#'
#' @description
#' This function can be used to compute the standard error of a 
#' semipartial correlation using the sample correlation, sample 
#' size, and squared multiple correlation for the full model.
#' The effect size estimate and standard error output from this 
#' function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions
#' in applications where a combination of different types of correlations
#' are used in the meta-analysis. 
#' 
#' 
#' @param    cor	sample semipartial correlation  
#' @param    r2		squared multiple correlation for full model  
#' @param    n		sample size
#'   
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of semipartial correlation
#' * SE - standard error
#' 
#' 
#' @examples
#' se.semipartial(.40, .25, 60)
#'
#' # Should return: 
#' #                          Estimate        SE
#' # Semipartial correlation:      0.4 0.1063262
#' 
#' 
#' @export
se.semipartial <- function(cor, r2, n) {
 r0 <- r2 - cor^2
 a <- r2^2 - 2*r2 + r0 - r0^2 + 1
 se.cor <- sqrt(a/(n - 3))
 out <- t(c(cor, se.cor))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- c("Semipartial correlation:")
 return(out)
}

# se.pbcor
#' Computes the standard error for a point-biserial correlation 
#' 
#'
#' @description
#' The function computes a point-biserial correlation and its standard 
#' error for two types of point-biserial correlations in 2-group designs
#' using the sample means, sample standard deviations, and samples sizes.
#' One type of point-biserial correlation uses an unweighted average of
#' variances and is appropriate for 2-group experimental designs. The 
#' other type of point-biserial correlation uses an weighted average of
#' variances and is appropriate for 2-group nonexperimental designs with
#' simple random sampling. This function is useful in a meta-analysis of
#' compatible point-biserial correlations where some studies used a 2-group 
#' experimental design and other studies used a 2-group nonexperimental 
#' design. The effect size estimate and standard error output from this 
#' function can  be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions.
#'
#' 
#' @param    m1		sample mean for group 1 
#' @param    m2		sample mean for group 2 
#' @param    sd1	sample standard deviation for group 1
#' @param    sd2	sample standard deviation for group 2
#' @param    n1		group 1 sample size
#' @param    n2		group 2 sample size
#' @param    type		
#' * set to 1 for weighted variance average
#' * set to 2 for unweighted variance average
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of point-biserial correlation
#' * SE - standard error
#' 
#' 
#' @examples
#' se.pbcor(21.9, 16.1, 3.82, 3.21, 40, 40, 1)
#'
#' #  Should return: 
#' #                               Estimate         SE
#' #  Point-biserial correlation: 0.6349786 0.05981325
#' 
#' 
#' @references
#' \insertRef{Bonett2020b}{vcmeta}
#' 
#' 
#' @export
se.pbcor <- function(m1, m2, sd1, sd2, n1, n2, type) {
  df1 <- n1 - 1
  df2 <- n2 - 1
  if (type == 1) {
    u <- n1/(n1 + n2)
    s <- sqrt((df1*sd1^2 + df2*sd2^2)/(df1 + df2))
    d <- (m1 - m2)/s
    c <- 1/(u*(1 - u))
    cor <- d/sqrt(d^2 + c)
    se.d <- sqrt(d^2*(1/df1 + 1/df2)/8 + 1/n1 + 1/n2)
    se.cor <- (c/(d^2 + c)^(3/2))*se.d                                                
  } else {
    s <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s
    cor <- d/sqrt(d^2 + 4)
    se.d <- sqrt(d^2*(sd1^4/df1 + sd2^4/df2)/(8*s^4) + (sd1^2/df1 + sd2^2/df2)/s^2) 
    se.cor <- (4/(d^2 + 4)^(3/2))*se.d                                                
  }
  out <- t(c(cor, se.cor))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Point-biserial correlation:")
  return(out)
}


# se.odds
#' Computes the standard error for a log odds ratio 
#' 
#'
#' @description 
#' This function computes a log odds ratio and its standard error using
#' the frequency counts and sample sizes in a 2-group design. These
#' frequency counts and sample sizes can be obtained from a 2x2 
#' contingency table. This function is useful in a meta-analysis of
#' odds ratios where some studies report the sample odds ratio and its
#' standard error and other studies only report the frequency counts
#' (or a 2x2 contingency table. The log odds ratio and standard error 
#' output from this function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions.
#' 
#' 
#' @param    f1		number of participants who have the outcome of interest in group 1 
#' @param    f2		number of participants who have the outcome of interest in group 2   
#' @param    n1		group 1 sample size
#' @param    n2		group 2 sample size
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of log odds ratio
#' * SE - standard error
#' 
#' 
#' @examples
#' se.odds(36, 50, 21, 50)
#'
#' # Should return: 
#' #                 Estimate        SE
#' # Log odds ratio: 1.239501 0.4204435
#' 
#' 
#' @references
#' \insertRef{Bonett2015}{vcmeta}
#' 
#' 
#' @export
se.odds <- function(f1, n1, f2, n2) {
  log.OR <- log((f1 + .5)*(n2 - f2 + .5)/((f2 + .5)*(n1 - f1 + .5)))
  se.log.OR <- sqrt(1/(f1 + .5) + 1/(f2 + .5) + 1/(n1 - f1 + .5) + 1/(n2 - f2 + .5))
  out <- t(c(log.OR, se.log.OR))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- "Log odds ratio:"
  return(out)
}


# se.meanratio2
#' Computes the standard error for a 2-group log mean ratio
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' 2-group mean ratio using the two sample means, sample standard
#' deviations, and sample sizes. The effect size estimate and standard 
#' error output from this function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' application where compatible mean ratios from a combination of 2-group
#' and paired-samples experiments are used in the meta-analysis.
#' 
#' 
#' @param    m1		sample mean for group 1 
#' @param    m2		sample mean for group 2 
#' @param    sd1	sample standard deviation for group 1
#' @param    sd2	sample standard deviation for group 2
#' @param    n1		group 1 sample size
#' @param    n2		group 2 sample size
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of log mean ratio
#' * SE - standard error
#' @return
#'  
#' @examples
#' se.meanratio2(21.9, 16.1, 3.82, 3.21, 40, 40)
#'
#' # Should return:
#' #                  Estimate       SE
#' # Log mean ratio: 0.3076674 0.041886
#' 
#' 
#' @references
#' \insertRef{Bonett2020}{vcmeta}
#' 
#' 
#' @export
se.meanratio2 <- function(m1, m2, sd1, sd2, n1, n2) {
  logratio <- log(m1/m2)
  var1 <- sd1^2/(n1*m1^2) 
  var2 <- sd2^2/(n2*m2^2)
  se <- sqrt(var1 + var2)
  out <- t(c(logratio, se))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- "Log mean ratio:"
  return(out)
}


# se.meanratio.ps
#' Computes the standard error for a paired-samples log mean ratio
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' paired-samples mean ratio using the two sample means, sample 
#  standard deviations, Pearson correlation, and sample size. The 
#' effect size estimate and standard error output from this function
#' can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' application where compatible mean ratios from a combination of 2-group
#' and paired-samples experiments are used in the meta-analysis. 
#' 
#' 
#' @param    m1		sample mean for measurement 1 
#' @param    m2		sample mean for measurement 2 
#' @param    sd1	sample standard deviation for measurement 1
#' @param    sd2	sample standard deviation for measurement 2
#' @param    cor	smple correlation for measurements 1 and 2 
#' @param    n		sample size
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of log mean ratio
#' * SE - standard error
#' @return
#'  
#' @examples
#' se.meanratio.ps(21.9, 16.1, 3.82, 3.21, .748, 40)
#'
#' # Should return:
#' #                  Estimate         SE
#' # Log mean ratio: 0.3076674 0.02130161
#' 
#' 
#' @references
#' \insertRef{Bonett2020}{vcmeta}
#' 
#' 
#' @export
se.meanratio.ps <- function(m1, m2, sd1, sd2, cor, n) {
 logratio <- log(m1/m2)
 var1 <- sd1^2/(n*m1^2) 
 var2 <- sd2^2/(n*m2^2)
 cov <- cor*sd1*sd2/(n*m1*m2)
 se <- sqrt(var1 + var2 - 2*cov)
 out <- t(c(logratio, se))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- "Log mean ratio:"
 return(out)
}


# se.slope
#' Computes a slope and standard error
#' 
#'
#' @description 
#' This function can be used to compute a slope and its standard error
#' for a simple linear regression model using the sample Pearson
#' correlation and the standard deviations of response and predictor 
#' variables. This function is useful in a meta-analysis of slopes of 
#' a simple linear regression model where some studies report the Pearson 
#' correlation but not the slope.
#' 
#'
#' @param    cor		sample Pearson correlation  
#' @param    sdy		sample standard deviation of the response variable
#' @param    sdx		sample standard deviation of the predictor variable
#' @param    n		  sample size
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimate of slope
#' * SE - standard error
#' 
#' 
#' @examples
#' se.slope(.392, 4.54, 2.89, 60)
#'
#' # Should return: 
#' #         Estimate        SE
#' # Slope: 0.6158062 0.1897647
#' 
#' 
#' @references
#' \insertRef{Snedecor1980}{vcmeta}
#'
#'
#' @export
se.slope <- function(cor, sdy, sdx, n) {
  slope <- cor*sdy/sdx
  se.slope <- sqrt((sdy^2*(1 - cor^2)*(n - 1))/(sdx^2*(n - 1)*(n - 2)))
  out <- t(c(slope, se.slope))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- c("Slope:")
  return(out)
}



