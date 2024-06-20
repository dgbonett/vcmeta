# se.mean2 =================================================================
#' Computes the standard error for a 2-group mean difference
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' 2-group mean difference using the estimated means, estimated 
#' standard deviations, and sample sizes. The effect size estimate 
#' and standard error output from this function can be used as input
#' in the \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, 
#' and \link[vcmeta]{meta.lm.gen} functions in applications where 
#' compatible mean differences from a combination of 2-group
#' and paired-samples experiments are used in the meta-analysis. 
#' Equality of variances is not asumed.
#' 
#' 
#' @param    m1		estimated mean for group 1 
#' @param    m2		estimated mean for group 2 
#' @param    sd1	estimated standard deviation for group 1
#' @param    sd2	estimated standard deviation for group 2
#' @param    n1		sample size for group 1
#' @param    n2		sample size for group 2
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated mean difference
#' * SE - standard error
#' 
#'  
#' @examples
#' se.mean2(21.9, 16.1, 3.82, 3.21, 40, 40)
#'
#  # Should return:
#' #                   Estimate        SE
#' # Mean difference:       5.8 0.7889312
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
  rownames(out) <- "Mean difference: "
  return(out)
}


# se.mean.ps =================================================================
#' Computes the standard error for a paired-samples mean difference
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' paired-samples mean difference using the estimated means, 
#' estimated standard deviations, estimated Pearson correlation, 
#' and sample size. The effect size estimate and standard error 
#' output from this function can be used as input in the
#' \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, 
#' and \link[vcmeta]{meta.lm.gen} functions in applications where 
#' compatible mean differences from a combination of 2-group
#' and paired-samples experiments are used in the meta-analysis. 
#' Equality of variances is not assumed.
#' 
#' 
#' @param    m1		estimated mean for measurement 1 
#' @param    m2		estimated mean for measurement 2 
#' @param    sd1	estimated standard deviation for measurement 1
#' @param    sd2	estimated standard deviation for measurement 2
#' @param    cor	estimated correlation for measurements 1 and 2
#' @param    n		sample size
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated mean difference
#' * SE - standard error
#' 
#' 
#' @examples
#' se.mean.ps(23.9, 25.1, 1.76, 2.01, .78, 25)
#'
#' # Should return:
#' #                   Estimate        SE
#' # Mean difference:      -1.2 0.2544833
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
  rownames(out) <- "Mean difference: "
  return(out)
}


# se.stdmean2 ================================================================		
#' Computes the standard error for a 2-group standardized mean difference
#' 
#' 
#' @description
#' This function computes the standard error of a 2-group standardized
#' mean difference using the sample sizes and the estimated means
#  and standard deviations. Use the square root average variance
#' standardizer (stdzr = 0) for 2-group experimental designs. Use the 
#' square root weighted variance standardizer (stdzr = 3) for 2-group 
#' nonexperimental designs with simple random sampling. The single-group 
#' standardizers (stdzr = 1 and stdzr = 2) can be used with either 
#' 2-group experimental or nonexperimental designs. The effect size 
#' estimate and standard error output from this function can be used as 
#' input in the \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen},
#' and \link[vcmeta]{meta.lm.gen} functions in applications where compatible
#' standardized mean differences from a combination of 2-group and 
#' paired-samples experiments are used in the meta-analysis. Equality 
#' of variances is not assumed.
#'
#' 
#' @param    m1		  estimated mean for group 1 
#' @param    m2		  estimated mean for group 2 
#' @param    sd1		estimated standard deviation for group 1
#' @param    sd2		estimated standard deviation for group 2
#' @param    n1		  sample size for group 1
#' @param    n2		  sample size for group 2
#' @param stdzr
#' * set to 0 for square root average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#' * set to 3 for square root weighted variance standardizer
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated standardized mean difference
#' * SE - standard error
#' 
#' 
#' @examples
#' se.stdmean2(21.9, 16.1, 3.82, 3.21, 40, 40, 0)
#'
#' # Should return: 
#' #                                Estimate        SE
#' # Standardized mean difference:  1.643894 0.2629049
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#'
#'
#' @seealso \link[vcmeta]{se.cohen}
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
  rownames(out) <- "Standardized mean difference: "
  return(out)
}


# se.stdmean.ps ==============================================================	
#' Computes the standard error for a paired-samples standardized mean 
#' difference
#' 
#'
#' @description
#' This function computes the standard error of a paired-samples standardized
#' mean difference using the sample size and estimated means, standard 
#' deviations, and estimated correlation. The effect size estimate and standard error
#' output from this function can be used as input in the \link[vcmeta]{meta.ave.gen},
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' applications where compatible standardized mean differences from a combination
#' of 2-group and paired-samples experiments are used in the meta-analysis. 
#' Equality of variances is not assumed.
#' 
#' 
#' @param    m1		estimated mean for measurement 1 
#' @param    m2		estimated mean for measurement 2 
#' @param    sd1	estimated standard deviation for measurement 1
#' @param    sd2	estimated standard deviation for measurement 2
#' @param    cor	estimated correlation for measurements 1 and 2
#' @param    n		sample size
#' @param stdzr
#' * set to 0 for square root average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#'  
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated standardized mean difference
#' * SE - standard error
#' 
#' 
#' @examples
#' se.stdmean.ps(23.9, 25.1, 1.76, 2.01, .78, 25, 0)
#'
#' # Should return: 
#' #                                   Estimate        SE
#' # Standardizedd mean difference:  -0.6352097 0.1602852
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
  rownames(out) <- "Standardized mean difference: "
  return(out)
}


# se.cor ==========================================================
#' Computes the standard error for a Pearson or partial correlation 
#' 
#'
#' @description
#' This function can be used to compute the standard error of a 
#' Pearson or partial correlation using the estimated correlation, 
#' sample size, and number of control variables. The correlation, 
#' along with the standard error output from this function, can be used 
#' as input in the \link[vcmeta]{meta.ave.gen} function in applications
#' where a combination of different types of correlations are used in
#' the meta-analysis. 
#' 
#' 
#' @param    cor	estimated Pearson or partial correlation  
#' @param    s		number of control variables (set to 0 for Pearson)  
#' @param    n		sample size
#'   
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - Pearson or partial correlation (from input)
#' * SE - standard error
#' 
#' 
#' @examples
#' se.cor(.40, 0, 55)
#'
#' # Should return: 
#' #               Estimate       SE
#' # Correlation:       0.4 0.116487
#' 
#' 
#' @references
#' \insertRef{Bonett2008a}{vcmeta}
#' 
#' 
#' @export
se.cor <- function(cor, s, n) {
  se.cor <- sqrt((1 - cor^2)^2/(n - 3 - s))
  out <- t(c(cor, se.cor))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- "Correlation: "
  return(out)
}


# se.spear ===================================================================
#' Computes the standard error for a Spearman correlation 
#' 
#' 
#' @description
#' This function can be used to compute the Bonett-Wright standard
#' error of a Spearman correlation using the estimated correlation
#' and sample size. The standard error from this function can be used
#' as input in the \link[vcmeta]{meta.ave.gen} function in applications 
#' where a combination of different types of correlations are used in 
#' the meta-analysis. 
#' 
#' 
#' @param    cor		estimated Spearman correlation  
#' @param    n		  sample size
#'   
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - Spearman correlation (from input)
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
  rownames(out) <- "Spearman correlation: "
  return(out)
}


# se.semipart ================================================================
#' Computes the standard error for a semipartial correlation 
#' 
#'
#' @description
#' This function can be used to compute the standard error of a 
#' semipartial correlation using the estimated correlation, sample 
#' size, and squared multiple correlation for the full model.
#' The effect size estimate and standard error output from this 
#' function can be used as input in the \link[vcmeta]{meta.ave.gen} 
#' function in applications where a combination of different types
#' of correlations are used in the meta-analysis. 
#' 
#' 
#' @param    cor	estimated semipartial correlation  
#' @param    r2		estimated squared multiple correlation for full model  
#' @param    n		sample size
#'   
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - semipartial correlation (from input)
#' * SE - standard error
#' 
#' 
#' @examples
#' se.semipartial(.40, .25, 60)
#'
#' # Should return: 
#' #                           Estimate        SE
#' # Semipartial correlation:       0.4 0.1063262
#' 
#' 
#' @export
se.semipartial <- function(cor, r2, n) {
 r0 <- r2 - cor^2
 a <- r2^2 - 2*r2 + r0 - r0^2 + 1
 se.cor <- sqrt(a/(n - 3))
 out <- t(c(cor, se.cor))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- "Semipartial correlation: "
 return(out)
}


# se.pbcor ==============================================================
#' Computes the standard error for a point-biserial correlation 
#' 
#'
#' @description
#' This function computes a point-biserial correlation and its standard 
#' error for two types of point-biserial correlations in 2-group designs
#' using the estimated means, estimated standard deviations, and samples
#' sizes. Equality of variances is not assumed. One type of point-biserial
#' correlation uses an unweighted average of variances and is recommended
#' for 2-group experimental designs. The other type of point-biserial 
#' correlation uses a weighted average of variances and is recommended for
#' 2-group nonexperimental designs with simple random sampling (but not 
#' stratified random sampling). This function is useful in a meta-analysis
#' of compatible point-biserial correlations where some studies used a 
#' 2-group experimental design and other studies used a 2-group 
#' nonexperimental design. The effect size estimate and standard error 
#' output from this function can  be used as input in the
#' \link[vcmeta]{meta.ave.gen} function.
#'
#' 
#' @param    m1		estimated mean for group 1 
#' @param    m2		estimated mean for group 2 
#' @param    sd1	estimated standard deviation for group 1
#' @param    sd2	estimated standard deviation for group 2
#' @param    n1		sample size for group 1
#' @param    n2		sample size for group 2
#' @param    type		
#' * set to 1 for weighted variance average
#' * set to 2 for unweighted variance average
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated point-biserial correlation
#' * SE - standard error
#' 
#' 
#' @examples
#' se.pbcor(21.9, 16.1, 3.82, 3.21, 40, 40, 1)
#'
#' #  Should return: 
#' #                                Estimate         SE
#' #  Point-biserial correlation:  0.6349786 0.05981325
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
  rownames(out) <- "Point-biserial correlation: "
  return(out)
}


# se.odds ====================================================================
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
#' or a 2x2 contingency table. The log odds ratio and standard error 
#' output from this function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions.
#' 
#' 
#' @param    f1		number of participants who have the outcome in group 1 
#' @param    f2		number of participants who have the outcome in group 2   
#' @param    n1		sample size for group 1
#' @param    n2		sample size for group 2 
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated log odds ratio
#' * SE - standard error
#' 
#' 
#' @examples
#' se.odds(36, 50, 21, 50)
#'
#' # Should return: 
#' #                  Estimate        SE
#' # Log odds ratio:  1.239501 0.4204435
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
  rownames(out) <- "Log odds ratio: "
  return(out)
}


# se.meanratio2 =========================================================
#' Computes the standard error for a 2-group log mean ratio
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' 2-group log mean ratio using the estimated means, estimated standard
#' deviations, and sample sizes. The log mean estimate and standard 
#' error output from this function can be used as input in the
#' \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, and
#' \link[vcmeta]{meta.lm.gen} functions in application where compatible
#' mean ratios from a combination of 2-group and paired-samples experiments
#' are used in the meta-analysis. Equality of variances is not assumed.
#' 
#' 
#' @param    m1		estimated mean for group 1 
#' @param    m2		estimated mean for group 2 
#' @param    sd1	estimated standard deviation for group 1
#' @param    sd2	estimated standard deviation for group 2
#' @param    n1		sample size for group 1
#' @param    n2		sample size for group 2
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated log mean ratio
#' * SE - standard error
#'  
#' @examples
#' se.meanratio2(21.9, 16.1, 3.82, 3.21, 40, 40)
#'
#' # Should return:
#' #                   Estimate       SE
#' # Log mean ratio:  0.3076674 0.041886
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
  rownames(out) <- "Log mean ratio: "
  return(out)
}


# se.meanratio.ps =============================================================
#' Computes the standard error for a paired-samples log mean ratio
#' 
#' 
#' @description
#' This function can be used to compute the standard error of a 
#' paired-samples log mean ratio using the estimated means, estimated
#'  standard deviations, estimated Pearson correlation, and sample 
#' size. The log-mean estimate and standard error output from
#' this function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' applications where compatible mean ratios from a combination of 2-group
#' and paired-samples experiments are used in the meta-analysis. 
#' Equality of variances is not assumed.
#' 
#' 
#' @param    m1		estimated mean for measurement 1 
#' @param    m2		estimated mean for measurement 2 
#' @param    sd1	estimated standard deviation for measurement 1
#' @param    sd2	estimated standard deviation for measurement 2
#' @param    cor	estimated correlation for measurements 1 and 2 
#' @param    n		sample size
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated log mean ratio
#' * SE - standard error
#'  
#' @examples
#' se.meanratio.ps(21.9, 16.1, 3.82, 3.21, .748, 40)
#'
#' # Should return:
#' #                   Estimate         SE
#' # Log mean ratio:  0.3076674 0.02130161
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
 rownames(out) <- "Log mean ratio: "
 return(out)
}


# se.slope =================================================================
#' Computes a slope and standard error
#' 
#'
#' @description 
#' This function can be used to compute a slope and its standard error
#' for a simple linear regression model (random-x model) using the estimated 
#' Pearson correlation and the estimated standard deviations of the 
#' response and predictor variables. This function is useful in a meta-analysis
#' of slopes of a simple linear regression model where some studies report
#' the Pearson correlation but not the slope.
#' 
#'
#' @param    cor		estimated Pearson correlation  
#' @param    sdy		estimated standard deviation of the response variable
#' @param    sdx		estimated standard deviation of the predictor variable
#' @param    n		  sample size
#'   
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated slope
#' * SE - standard error
#' 
#' 
#' @examples
#' se.slope(.392, 4.54, 2.89, 60)
#'
#' # Should return: 
#' #          Estimate        SE
#' # Slope:  0.6158062 0.1897647
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
  rownames(out) <- "Slope: "
  return(out)
}


# se.prop2 =================================================================== 
#' Computes the estimate and standard error for a 2-group proportion 
#' difference
#' 
#' 
#' @description
#' This function can be used to compute the Agresti-Caffo standard 
#' error of a 2-group proportion difference using the frequency 
#' counts and sample sizes. The effect size estimate and standard 
#' error output from this function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' applications where compatible proportion differences from a combination of 
#' 2-group and paired-samples studies are used in the meta-analysis. 
#' 
#' 
#' @param    f1   number of participants in group 1 who have the outcome
#' @param    f2	  number of participants in group 2 who have the outcome
#' @param    n1	  sample size for group 1
#' @param    n2	  sample size for group 2
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated proportion difference
#' * SE - standard error
#' 
#'  
#' @examples
#' se.prop2(31, 16, 40, 40)
#'
#' # Should return:
#' #                          Estimate        SE
#' # Proportion difference:  0.3571429 0.1002777
#' 
#' 
#' @references
#' \insertRef{Agresti2000}{vcmeta}
#'
#'
#' @export
se.prop2 <- function(f1, f2, n1, n2) {
 p1 <- (f1 + 1)/(n1 + 2)
 p2 <- (f2 + 1)/(n2 + 2)
 est <- p1 - p2
 se <- sqrt(p1*(1 - p1)/(n1 + 2) + p2*(1 - p2)/(n2 + 2))
 out <- t(c(est, se))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- "Proportion difference: "
 return(out)
}


# se.prop.ps ==============================================================
#' Computes the estimate and standard error for a paired-samples
#' proportion difference
#' 
#' 
#' @description
#' This function can be used to compute the Bonett-Price standard error 
#' of a paired-samples proportion difference using the frequency counts 
#' from a 2 x 2 contingency table. The effect size estimate and standard 
#' error output from this function can be used as input in the \link[vcmeta]{meta.ave.gen}, 
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' applications where compatible proportion differences from a combination of
#' 2-group and paired-samples studies are used in the meta-analysis. 
#' 
#' 
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated proportion difference
#' * SE - standard error
#' 
#'  
#' @examples
#' se.prop.ps(16, 64, 5, 15)
#'
#' # Should return:
#' #                          Estimate         SE
#' # Proportion difference:  0.5784314 0.05953213
#' 
#' 
#' @references
#' \insertRef{Bonett2012}{vcmeta}
#'
#'
#' @export
se.prop.ps <- function(f00, f01, f10, f11) {
 n <- f00 + f01 + f10 + f11
 p01 <- (f01 + 1)/(n + 2)
 p10 <- (f10 + 1)/(n + 2)
 est <- p01 - p10
 se <- sqrt(((p01 + p10) - (p01 - p10)^2)/(n + 2))
 out <- t(c(est, se))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- "Proportion difference: "
 return(out)
}


# se.ave.mean2.dep ============================================================
#' Computes the standard error for the average of 2-group mean differences from 
#' two parallel measurement response variables in the same sample 
#'                       
#' 
#' @description
#' In a study that reports a 2-group mean difference for two response
#' variables that satisfy the conditions of parallel measurments, this function
#' can be used to compute the standard error of the average of the two mean 
#' differences using the two estimated means, estimated standard deviations, 
#' estimated within-group correlation between the two response variables, and 
#' the two sample sizes. The average mean difference and standard error output 
#' from this function can then be used as input in the
#' \link[vcmeta]{meta.ave.gen}, \link[vcmeta]{meta.lc.gen}, and 
#' \link[vcmeta]{meta.lm.gen} functions in a meta-analysis where some studies
#' have used one of the two parallel response variables and other studies have
#' used the other parallel response variable. Equality of variances is not
#' assumed.
#' 
#' 
#' @param    m1A	   estimated mean for variable A in group 1 
#' @param    m2A	   estimated mean for variable A in group 2 
#' @param    sd1A	   estimated standard deviation for variable A in group 1
#' @param    sd2A	   estimated standard deviation for variable A in group 2
#' @param    m1B	   estimated mean for variable B in group 1 
#' @param    m2B	   estimated mean for variable B in group 2 
#' @param    sd1B	   estimated standard deviation for variable B in group 1
#' @param    sd2B	   estimated standard deviation for variable B in group 2
#' @param    rAB	   estimated within-group correlation between variables A and B
#' @param    n1		   sample size for group 1
#' @param    n2		   sample size for group 2
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated average mean difference
#' * SE - standard error 
#' * VAR(A) - variance of mean difference for variable A
#' * VAR(B) - variance of mean difference for variable B
#' * COV(A,B) - covariance of mean differences for variables A and B
#' 
#'  
#' @examples
#' se.ave.mean2.dep(21.9, 16.1, 3.82, 3.21, 24.8, 17.1, 3.57, 3.64, .785, 40, 40)
#'
#' # Should return:
#' #                          Estimate        SE    VAR(A)    VAR(B)  COV(A,B)
#' # Average mean difference:     6.75 0.7526878 0.6224125 0.6498625 0.4969403
#' 
#' 
#' @export
se.ave.mean2.dep <- function(m1A, m2A, sd1A, sd2A, m1B, m2B, sd1B, sd2B, rAB, n1, n2) {
  m1 <- (m1A + m1B)/2
  m2 <- (m2A + m2B)/2
  est <- m1 - m2
  v1 <- sd1A^2/n1 + sd2A^2/n2 
  v2 <- sd1B^2/n1 + sd2B^2/n2
  cov <- rAB*sd1A*sd1B/n1 + rAB*sd2A*sd2B/n2
  se <- sqrt((v1 + v2 + 2*cov)/4)
  out <- t(c(est, se, v1, v2, cov))
  colnames(out) <- c("Estimate", "SE", "VAR(A)", "VAR(B)", "COV(A,B)")
  rownames(out) <- "Average mean difference: "
  return(out)
}


# se.ave.cor.over =============================================================
#' Computes the standard error for the average of two Pearson correlations with 
#' one variable in common that have been estimated from the same sample 
#'     
#' 
#' @description
#' In a study that reports the sample size and three correlations (cor12, cor13, 
#' and cor23 where variable 1 is called the "overlapping" variable), and 
#' variables 2 and 3 are different measurements of the same attribute, this 
#' function can be used to compute the average of cor12 and cor13 and its 
#' standard error. The average correlation and the standard error from this 
#' function can be used as input in the \link[vcmeta]{meta.ave.gen} function
#' in a meta-analysis where some studies have reported cor12 and other studies
#' have reported cor13. 
#' 
#' 
#' @param    cor12	estimated correlation between variables 1 and 2 
#' @param    cor13	estimated correlation between variables 1 and 3 
#' @param    cor23	estimated correlation between variables 2 and 3
#' @param    n		  sample size
#' 
#' 
#' @return
#' Returns a two-row matrix. The first row gives results for the average 
#' correlation and the second row gives the results with a Fisher
#' transformation. The columns are:
#' * Estimate - estimated average of cor12 and cor13
#' * SE - standard error 
#' * VAR(cor12) - variance of cor12 
#' * VAR(cor13) - variance of cor13
#' * COV(cor12,cor13) - covariance of cor12 and cor13
#' 
#'  
#' @examples
#' se.ave.cor.over(.462, .518, .755, 100)
#'
#' # Should return:
#' #                Estimate         SE  VAR(cor12) VAR(cor13) COV(cor12,cor13)
#' # Correlation:  0.4900000 0.07087351 0.006378045 0.00551907      0.004097553
#' # Fisher:       0.5360603 0.09326690 0.010309278 0.01030928      0.007119936
#'
#'
#' @export
se.ave.cor.over <- function(cor12, cor13, cor23, n) {
  est1 <- (cor12 + cor13)/2
  cov1 <- ((cor23 - cor12*cor13/2)*(1 - cor12^2 - cor13^2 - cor23^2) + cor23^3)/(n - 3)
  v1 <- (1 - cor12^2)^2/(n - 3)
  v2 <- (1 - cor13^2)^2/(n - 3)
  se1 <- sqrt((v1 + v2 + 2*cov1)/4)
  est2 <- log((1 + est1)/(1 - est1))/2 
  se2 <- se1/(1 - est1^2)
  cov2 <- cov1/((1 - cor12^2)*(1 - cor13^2))
  v1.z <- 1/(n - 3)
  v2.z <- 1/(n - 3)
  out1 <- t(c(est1, se1, v1, v2, cov1))
  out2 <- t(c(est2, se2, v1.z, v2.z, cov2))
  out <- rbind(out1, out2)
  colnames(out) <- c("Estimate", "SE", "VAR(cor12)", "VAR(cor13)", "COV(cor12,cor13)")
  rownames(out) <- c("Correlation: ", "Fisher: ")
  return(out)
}


# se.ave.cor.nonover ==========================================================
#' Computes the standard error for the average of two Pearson correlations with 
#' no variables in common that have been estimated from the same sample 
#'
#' 
#' @description
#' In a study that reports the sample size and six correlations (cor12, cor34,
#' cor13, cor14, cor23, and cor24) where variables 1 and 3 are different 
#' measurements of one attribute and variables 2 and 4 are different 
#' measurements of a second attribute, this function can be used to compute the 
#' average of cor12 and cor34 and its standard error. Note that cor12 and cor34
#' have no variable in common (i.e., no "overlapping" variable). The average 
#' correlation and the standard error from this function can be used as 
#' input in the \link[vcmeta]{meta.ave.gen} function in a meta-analysis where
#' some studies have reported cor12 and other studies have reported cor34. 
#' 
#' 
#' @param    cor12	  estimated correlation between variables 1 and 2 
#' @param    cor34	  estimated correlation between variables 3 and 4 
#' @param    cor13	  estimated correlation between variables 1 and 3
#' @param    cor14	  estimated correlation between variables 1 and 4
#' @param    cor23	  estimated correlation between variables 2 and 3
#' @param    cor24	  estimated correlation between variables 2 and 4
#' @param    n		    sample size
#' 
#' 
#' @return
#' Returns a two-row matrix. The first row gives results for the average 
#' correlation and the second row gives the results with a Fisher
#' transformation. The columns are:
#' * Estimate - estimated average of cor12 and cor34
#' * SE - standard error 
#' * VAR(cor12) - variance of cor12 
#' * VAR(cor34) - variance of cor34
#' * COV(cor12,cor34) - covariance of cor12 and cor34
#' 
#'  
#' @examples
#' se.ave.cor.nonover(.357, .398, .755, .331, .347, .821, 100)
#'
#' # Should return:
#' #               Estimate         SE VAR(cor12)  VAR(cor34) COV(cor12,cor34)
#' # Correlation:  0.377500 0.07768887 0.00784892 0.007301895      0.004495714
#' # Fisher:       0.397141 0.09059993 0.01030928 0.010309278      0.006122153
#'
#'
#' @export
se.ave.cor.nonover <- function(cor12, cor34, cor13, cor14, cor23, cor24, n) {
  est1 <- (cor12 + cor34)/2
  c1 <- (cor12*cor34)*(cor13^2 + cor14^2 + cor23^2 + cor24^2)/2 + cor13*cor24 + cor14*cor23
  c2 <- (cor12*cor13*cor14 + cor12*cor23*cor24 + cor13*cor23*cor34 + cor14*cor24*cor34)
  cov1 <- (c1 - c2)/(n - 3)
  v1 <- (1 - cor12^2)^2/(n - 3)
  v2 <- (1 - cor34^2)^2/(n - 3)
  se1 <- sqrt((v1 + v2 + 2*cov1)/4)
  est2 <- log((1 + est1)/(1 - est1))/2 
  se2 <- se1/(1 - est1^2)
  cov2 <- cov1/((1 - cor12^2)*(1 - cor34^2))
  v1.z <- 1/(n - 3)
  v2.z <- 1/(n - 3)
  out1 <- t(c(est1, se1, v1, v2, cov1))
  out2 <- t(c(est2, se2, v1.z, v2.z, cov2))
  out <- rbind(out1, out2)
  colnames(out) <- c("Estimate", "SE", "VAR(cor12)", "VAR(cor34)", "COV(cor12,cor34)")
  rownames(out) <- c("Correlation: ", "Fisher: ")
  return(out)
}


#  se.tetra ==================================================================
#' Computes the standard error for a tetrachoric correlation approximation  
#'
#'
#' @description
#' This function can be used to compute an estimate of a tetrachoric 
#' correlation approximation and its standard error using the frequency counts
#' from a 2 x 2 contingency table for two artifically dichotomous variables.
#' A tetrachoric approximation could be compatible with a Pearson correlation
#' in a meta-analysis. The tetrachoric approximation and the standard error 
#' from this function can be used as input in the \link[vcmeta]{meta.ave.gen} 
#' function in a meta-analysis where some studies have reported Pearson 
#' correlations between quantitative variables x and y and other studies have
#' reported a 2 x 2 contingency table for dichotomous measurements of variables
#' x and y. 
#'
#'
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#'
#'
#' @references
#' \insertRef{Bonett2005}{vcmeta}
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated tetrachoric approximation
#' * SE - standard error
#'
#'
#' @examples
#' se.tetra(46, 15, 54, 85)
#'
#' # Should return:
#' #                Estimate         SE 
#' # Tetrachoric:  0.5135167 0.09358336
#'
#'
#' @export
se.tetra <- function(f00, f01, f10, f11) {
 n <- f00 + f01 + f10 + f11
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 r1 <- (f00 + f01 + 1)/(n + 2)
 r2 <- (f10 + f11 + 1)/(n + 2)
 c1 <- (f00 + f10 + 1)/(n + 2)
 c2 <- (f01 + f11 + 1)/(n + 2)
 pmin <- min(c1, c2, r1, r2)
 c <- (1 - abs(r1 - c1)/5 - (.5 - pmin)^2)/2
 lor <- log(or)
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 tetra <- cos(3.14159/(1 + or^c))
 k <- (3.14159*c*or^c)*sin(3.14159/(1 + or^c))/(1 + or^c)^2
 se <- k*se.lor
 out <- t(c(tetra, se))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- "Tetrachoric: "
 return(out)
}


#  se.biphi ==================================================================
#' Computes the standard error for a biserial-phi correlation  
#'
#'
#' @description
#' This function can be used to compute an estimate of a biserial-phi  
#' correlation and its standard error using the frequency counts from a 2 x 2
#' contingency table where one variable is naturally dichotomous and the other
#' variable is artifically dichotomous. A biserial-phi correlation could be 
#' compatible with a point-biserial correlation in a meta-analysis. The 
#' biserial-phi estimate and the standard error from this function can be used 
#' as input in the \link[vcmeta]{meta.ave.gen} function in a meta-analysis 
#' where a point-biserial correlation has been obtained in some studies and
#' a biserial-phi correlation has been obtained in other studies.  
#'
#'
#' @param   f1     number of participants in group 1 who have the attribute
#' @param   f2     number of participants in group 2 who have the attribute
#' @param   n1     sample size for group 1
#' @param   n2     sample size for group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated biserial-phi correlation 
#' * SE - standard error
#'
#'
#' @examples
#' se.biphi(34, 22, 50, 50)
#'
#' # Should return:
#' #               Estimate        SE 
#' # Biserial-phi:  0.27539 0.1074594
#'
#'
#' @export
se.biphi <- function(f1, f2, n1, n2) {
 if (f1 > n1) {stop("f cannot be greater than n")}
 if (f2 > n2) {stop("f cannot be greater than n")}
 f00 <- f1
 f10 <- n1 - f1
 f01 <- f2
 f11 <- n2 - f2
 p1 <- n1/(n1 + n2)
 p2 <- n2/(n1 + n2)
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 lor <- log(or)
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 c <- 2.89/(p1*p2)
 biphi <- lor/sqrt(lor^2 + c)
 se.biphi <- sqrt(c^2/(lor^2 + c)^3)*se.lor
 out <- t(c(biphi, se.biphi))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- "Biserial-phi: "
 return(out)
}


# se.cohen ====================================================================		
#' Computes the standard error for Cohen's d
#' 
#' 
#' @description
#' This function computes the standard error of Cohen's d using only the two 
#' sample sizes and an estimate of Cohen's d. Cohen's d and its standard error 
#' assume equal variances. The estimate of Cohen's d, with the standard error
#' output from this function, can be used as input in the \link[vcmeta]{meta.ave.gen},
#' \link[vcmeta]{meta.lc.gen}, and \link[vcmeta]{meta.lm.gen} functions in 
#' applications where different types of compatible standardized mean 
#' differences are used in the meta-analysis. 
#'
#' 
#' @param    d		  estimated Cohen's d
#' @param    n1		  sample size for group 1
#' @param    n2		  sample size for group 2
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - Cohen's d (from input)
#' * SE - standard error
#' 
#' 
#' @examples
#' se.cohen(.78, 35, 50)
#'
#' # Should return: 
#' #            Estimate        SE
#' # Cohen's d:     0.78 0.2288236
#'
#'
#' @seealso \link[vcmeta]{se.stdmean2}
#'
#'
#' @export
se.cohen <- function(d, n1, n2) {
  df1 <- n1 - 1
  df2 <- n2 - 1
  se <- sqrt(d^2*(1/df1 + 1/df2)/8 + 1/n1 + 1/n2)
  out <- t(c(d, se))
  colnames(out) <- c("Estimate", "SE")
  rownames(out) <- "Cohen's d: "
  return(out)
}


#  se.bscor ===================================================================
#' Computes the standard error for a biserial correlation 
#' 
#'
#' @description
#' This function computes a biserial correlation and its standard error. A
#' biserial correlation can be used when one variable is quantitative and the 
#' other variable has been artifically dichotmized. The biserial correlation
#' estimates the correlation between an observable quantitative variable and
#' an unobserved quantitative variable that is measured on a dichotomous
#' scale. This function requires the estimated mean, estimated standard 
#' deviation, and samples size from each level of the dichotomized variable. 
#' This function is useful in a meta-analysis of Pearson correlations where
#' some studies report a Pearson correlation and other studies report the
#' information needed to compute a biserial correlation. The biserial 
#' correlation and standard error output from this function can be used as 
#' input in the \link[vcmeta]{meta.ave.gen} function.
#'
#'
#' @details
#' This function computes a point-biserial correlation and its standard error
#' as a function of a standardized mean difference with a weighted variance
#' standardizer. Then the point-biserial estimate is transformed into a 
#' biserial correlation using the traditional adjustment. The adjustment is 
#' also applied to the point-biserial standard error to obtain the standard 
#' error for the biserial correlation. 
#' 
#' The biserial correlation assumes that the observed quantitative variable 
#' and the unobserved quantitative variable have a bivariate normal 
#' distribution. Bivariate normality is a crucial assumption underlying the
#' transformation of a point-biserial correlation to a biserial correlation.
#' Bivariate normality also implies equal variances of the observed 
#' quantitative variable at each level of the dichotomized variable, and this
#' assumption is made in the computation of the standard error.
#'
#' 
#' @param    m1		estimated mean for level 1 
#' @param    m2		estimated mean for level 2
#' @param    sd1	estimated standard deviation for level 1
#' @param    sd2	estimated standard deviation for level 2
#' @param    n1		sample size for level 1
#' @param    n2		sample size for level 2
#' 
#' 
#' @return
#' Returns a one-row matrix:
#' * Estimate - estimated biserial correlation
#' * SE - standard error
#' 
#' 
#' @examples
#' se.bscor(21.9, 16.1, 3.82, 3.21, 40, 40)
#'
#' #  Should return: 
#' #                          Estimate         SE
#' #  Biserial correlation:  0.8018318 0.07451665
#' 
#' 
#' @references
#' \insertRef{Bonett2020b}{vcmeta}
#' 
#' 
se.bscor <- function(m1, m2, sd1, sd2, n1, n2) {
 df1 <- n1 - 1
 df2 <- n2 - 1
 u <- n1/(n1 + n2)
 a <- sqrt(u*(1 - u))/dnorm(qnorm(u))
 s <- sqrt((df1*sd1^2 + df2*sd2^2)/(df1 + df2))
 d <- (m1 - m2)/s
 c <- (df1 + df2)/((n1 + n2)*(u*(1 - u)))
 pbcor <- d/sqrt(d^2 + c)
 bscor <- pbcor*a
 if (bscor > 1) {bscor = .99999}
 if (bscor < -1) {bscor = -.99999}
 se.d <- sqrt(d^2*(1/n1 + 1/n2)/8 + 1/n1 + 1/n2)
 se.pbcor <- (c/(d^2 + c)^(3/2))*se.d  
 se.bscor <- se.pbcor*a  
 out <- t(c(bscor, se.bscor))
 colnames(out) <- c("Estimate", "SE")
 rownames(out) <- "Biserial correlation: "
 return(out)
}


