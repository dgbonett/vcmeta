#  meta.ave.mean2 ==========================================================
#' Confidence interval for an average mean difference from 2-group studies 
#'
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average mean difference from two or more 2-group studies. A Satterthwaite 
#' adjustment to the degrees of freedom is used to improve the accuracy of the
#' confidence intervals. Equality of variances within or across studies is not
#' assumed.
#'
#'  
#' @param alpha  	 alpha level for 1-alpha confidence
#' @param m1     	 vector of estimated means for group 1 
#' @param m2     	 vector of estimated means for group 2 
#' @param sd1    	 vector of estimated SDs for group 1
#' @param sd2    	 vector of estimated SDs for group 2
#' @param n1     	 vector of group 1 sample sizes
#' @param n2     	 vector of group 2 sample sizes
#' @param bystudy  logical to also return each study estimate (TRUE) or not
#'
#'
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy 
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @references 
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @examples
#' m1 <- c(7.4, 6.9)
#' m2 <- c(6.3, 5.7)
#' sd1 <- c(1.72, 1.53)
#' sd2 <- c(2.35, 2.04)
#' n1 <- c(40, 60)
#' n2 <- c(40, 60)
#' meta.ave.mean2(.05, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
#'
#' # Should return:
#' #         Estimate        SE        LL       UL     df
#' # Average     1.15 0.2830183 0.5904369 1.709563 139.41
#' # Study 1     1.10 0.4604590 0.1819748 2.018025  71.47
#' # Study 2     1.20 0.3292036 0.5475574 1.852443 109.42
#' 
#' 
#' @importFrom stats qt
#' @export
meta.ave.mean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE) {
  m <- length(m1)
  v1 <- sd1^2
  v2 <- sd2^2
  var <- v1/n1 + v2/n2
  d <- m1 - m2
  ave <- sum(d)/m
  se <- sqrt(sum(var)/m^2)
  u1 <- sum(var)^2
  u2 <- sum(v1^2/(n1^3 - n1^2) + v2^2/(n2^3 - n2^2))
  df <- round(u1/u2, 2)
  t <- qt(1 - alpha/2, df)
  ll <- ave - t*se
  ul <- ave + t*se
  out <- cbind(ave, se, ll, ul, df)
  row <- "Average"
  if (bystudy) {
    se <- sqrt(var)
    u1 <- var^2
    u2 <- v1^2/(n1^3 - n1^2) + v2^2/(n2^3 - n2^2)
    df <- u1/u2
    t <- qt(1 - alpha/2, df)
    ll <- d - t*se
    ul <- d + t*se
    row2 <- t(t(paste(rep("Study", m), seq(1, m))))
    row <- rbind(row, row2)
    out2 <- cbind(d, se, ll, ul, df)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) <- row 
  return(out)
}


#  meta.ave.stdmean2 ==========================================================
#' Confidence interval for an average standardized mean difference
#' from 2-group studies
#'
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average standardized mean difference from two or more 2-group studies.
#' Square root unweighted variances, square root weighted variances, and 
#' single group standard deviation are options for the standardizer. 
#' Equality of variances within or across studies is not assumed.
#'
#'
#' @param alpha	  alpha level for 1-alpha confidence
#' @param m1		  vector of estimated means for group 1
#' @param m2		  vector of estimated means for group 2
#' @param sd1		  vector of estimated SDs for group 1
#' @param sd2		  vector of estimated SDs for group 2
#' @param n1		  vector of group 1 sample sizes
#' @param n2		  vector of group 2 sample sizes
#' @param stdzr
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#' * set to 3 for square root weighted average variance standardizer
#' @param bystudy  logical to also return each study estimate (TRUE) or not
#'
#'
#' @return
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references 
#' \insertRef{Bonett2009a}{vcmeta}
#'
#'
#' @examples
#' m1 <- c(21.9, 23.1, 19.8)
#' m2 <- c(16.1, 17.4, 15.0)
#' sd1 <- c(3.82, 3.95, 3.67)
#' sd2 <- c(3.21, 3.30, 3.02)
#' n1 <- c(40, 30, 24)
#' n2 <- c(40, 28, 25)
#' meta.ave.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, 0, bystudy = TRUE)
#'
#' # Should return: 
#' #         Estimate      SE     LL     UL
#' # Average   1.5261 0.17343 1.1862 1.8661
#' # Study 1   1.6439 0.26290 1.1286 2.1592
#' # Study 2   1.5661 0.30563 0.9671 2.1652
#' # Study 3   1.4283 0.32892 0.7836 2.0729
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.stdmean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, stdzr, bystudy = TRUE) {
  df1 <- n1 - 1
  df2 <- n2 - 1
  m <- length(m1)
  z <- qnorm(1 - alpha/2)
  v1 <- sd1^2
  v2 <- sd2^2
  if (stdzr == 0) {
    s1 <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s1
    du <- (1 - 3/(4*(n1 + n2) - 9))*d
    ave <- sum(du)/m
    var <- d^2*(v1^2/df1 + v2^2/df2)/(8*s1^4) + (v1/df1 + v2/df2)/s1^2
    se <- sqrt(sum(var)/m^2)
  } else if (stdzr == 1) { 
    d <- (m1 - m2)/sd1
    du <- (1 - 3/(4*n1 - 5))*d
    ave <- sum(du)/m
    var <- d^2/(2*df1) + 1/df1 + v2/(df2*v1)
    se <- sqrt(sum(var)/m^2)
  } else if (stdzr == 2) {
    d <- (m1 - m2)/sd2
    du <- (1 - 3/(4*n2 - 5))*d
    ave <- sum(du)/m
    var <- d^2/(2*df2) + 1/df2 + v1/(df1*v2)
    se <- sqrt(sum(var)/m^2)
  } else {
    s2 <- sqrt((df1*v1 + df2*v2)/(df1 + df2))
    d <- (m1 - m2)/s2
    du <- (1 - 3/(4*(n1 + n2) - 9))*d
    ave <- sum(du)/m
    var <- d^2*(1/df1 + 1/df2)/8 + (v1/n1 + v2/n2)/s2^2
    se <- sqrt(sum(var)/m^2)
  }
  ll <- ave - z*se
  ul <- ave + z*se
  out <- cbind(round(ave, 4), round(se, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    if (stdzr == 0) {
      se <- sqrt(d^2*(v1^2/df1 + v2^2/df2)/(8*s1^4) + (v1/df1 + v2/df2)/s1^2)
    } else if (stdzr == 1) { 
      se <- sqrt(d^2/(2*df1) + 1/df1 + v2/(df2*v1))
    } else if (stdzr == 2) {
      se <- sqrt(d^2/(2*df2) + 1/df2 + v1/(df1*v2))
    } else {
      se <- sqrt(d^2*(1/df1 + 1/df2)/8 + (v1/n1 + v2/n2)/s2^2)
    }
    ll <- d - z*se
    ul <- d + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(d, 4), round(se, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return(out)
}


#  meta.ave.mean.ps ==========================================================
#' Confidence interval for an average mean difference from paired-samples studies
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average mean difference from two or more paired-samples studies.
#' A Satterthwaite adjustment to the degrees of freedom is used to improve 
#' the accuracy of the confidence interval for the average effect size. 
#' Equality of variances within or across studies is not assumed.
#'
#'
#' @param   alpha		 alpha level for 1-alpha confidence
#' @param   m1		   vector of estimated means for measurement 1 
#' @param   m2		   vector of estimated means for measurement 2 
#' @param   sd1		   vector of estimated SDs for measurement 1
#' @param   sd2		   vector of estimated SDs for measurement 2
#' @param   cor		   vector of estimated correlations for paired measurements
#' @param   n		     vector of sample sizes
#' @param   bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @examples
#' m1 <- c(53, 60, 53, 57)
#' m2 <- c(55, 62, 58, 61)
#' sd1 <- c(4.1, 4.2, 4.5, 4.0)
#' sd2 <- c(4.2, 4.7, 4.9, 4.8)
#' cor <- c(.72, .78, .81, .85)
#' n <- c(30, 50, 30, 70)
#' meta.ave.mean.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
#' 
#' # Should return:
#' #         Estimate        SE        LL         UL     df
#' # Average    -3.25 0.2340603 -3.713965 -2.7860352 107.66
#' # Study 1    -2.00 0.5672507 -3.160158 -0.8398421  29.00
#' # Study 2    -2.00 0.4227434 -2.849535 -1.1504653  49.00
#' # Study 3    -5.00 0.5335104 -6.091151 -3.9088487  29.00
#' # Study 4    -4.00 0.3023716 -4.603215 -3.3967852  69.00
#' 
#' 
#' @importFrom stats qt
#' @export
meta.ave.mean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, bystudy = TRUE) {
  m <- length(m1)
  v1 <- sd1^2
  v2 <- sd2^2
  v <- (v1 + v2 - 2*cor*sd1*sd2)/n
  d <- m1 - m2
  ave <- sum(d)/m
  se <- sqrt(sum(v)/m^2)
  u1 <- sum(v)^2
  u2 <- sum(v^2/(n - 1))
  df <- round(u1/u2, 2)
  t <- qt(1 - alpha/2, df)
  ll <- ave - t*se
  ul <- ave + t*se
  out <- cbind(ave, se, ll, ul, df)
  row <- "Average"
  if (bystudy) {
    se <- sqrt(v)
    df <- n - 1
    t <- qt(1 - alpha/2, df)
    ll <- d - t*se
    ul <- d + t*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(d, se, ll, ul, df)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) <- row 
  return(out)
}


#  meta.ave.stdmean.ps ==========================================================
#' Confidence interval for an average standardized mean difference from 
#' paired-samples studies  
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average standardized mean difference from two or more paired-samples
#' studies. Squrare root Unweighted variances and a single condition standard
#' deviation are options for the standardizer. Equality of variances within
#' or across studies is not assumed.
#'
#'
#' @param   alpha		alpha level for 1-alpha confidence
#' @param   m1		  vector of estimated means for measurement 1 
#' @param   m2		  vector of estimated means for measurement 2 
#' @param   sd1		  vector of estimated SDs for measurement 1
#' @param   sd2		  vector of estimated SDs for measurement 2
#' @param   cor		  vector of estimated correlations for paired measurements
#' @param   n		    vector of sample sizes
#' @param   stdzr		
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for measurement 1 SD standardizer 
#' * set to 2 for measurement 2 SD standardizer 
#' @param   bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @examples
#' m1 <- c(23.9, 24.1)
#' m2 <- c(25.1, 26.9)
#' sd1 <- c(1.76, 1.58)
#' sd2 <- c(2.01, 1.76)
#' cor <- c(.78, .84)
#' n <- c(25, 30)
#' meta.ave.stdmean.ps(.05, m1, m2, sd1, sd2, cor, n, 1, bystudy = TRUE)
#' 
#' # Should return: 
#' #         Estimate      SE      LL      UL
#' # Average  -1.1931 0.15680 -1.5004 -0.8858
#' # Study 1  -0.6818 0.17738 -1.0295 -0.3342
#' # Study 2  -1.7722 0.25862 -2.2790 -1.2653
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.stdmean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, stdzr, bystudy = TRUE) {
  df <- n - 1
  m <- length(m1)
  z <- qnorm(1 - alpha/2)
  v1 <- sd1^2
  v2 <- sd2^2
  vd <- v1 + v2 - 2*cor*sd1*sd2
  if (stdzr == 0) {
    s <- sqrt((sd1^2 + sd2^2)/2) 
    d <- (m1 - m2)/s
    du <- sqrt((n - 2)/df)*d
    ave <- sum(du)/m
    var <- d^2*(v1^2 + v2^2 + 2*cor^2*v1*v2)/(8*df*s^4) + vd/(df*s^2)
    se <- sqrt(sum(var)/m^2)
  } 
  else if (stdzr == 1) { 
    d <- (m1 - m2)/sd1
    du <- (1 - 3/(4*df - 1))*d
    ave <- sum(du)/m
    var <- d^2/(2*df) + vd/(df*v1)
    se <- sqrt(sum(var)/m^2)
  } 
  else {
    d <- (m1 - m2)/sd2
    du <- (1 - 3/(4*df - 1))*d
    ave <- sum(du)/m
    var <- d^2/(2*df) + vd/(df*v2)
    se <- sqrt(sum(var)/m^2)
  }  
  ll <- ave - z*se
  ul <- ave + z*se
  out <- cbind(round(ave, 4), round(se, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    if (stdzr == 0) {
      se <- sqrt(d^2*(v1^2 + v2^2 + 2*cor^2*v1*v2)/(8*df*s^4) + vd/(df*s^2))
    } else if (stdzr == 1) { 
      se <- sqrt(d^2/(2*df) + vd/(df*v1))
    } else {
      se <- sqrt(d^2/(2*df) + vd/(df*v2))
    }
    ll <- d - z*se
    ul <- d + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(d, 4), round(se, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return(out)
}


#  meta.ave.meanratio2 ==========================================================
#' Confidence interval for an average mean ratio from 2-group studies 
#'
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' geometric average mean ratio from two or more 2-group studies. A Satterthwaite 
#' adjustment to the degrees of freedom is used to improve the accuracy of the
#' confidence intervals. Equality of variances within or across studies is not assumed.
#'
#'
#' @param   alpha  	 alpha level for 1-alpha confidence
#' @param   m1     	 vector of estimated means for group 1
#' @param   m2     	 vector of estimated means for group 2
#' @param   sd1    	 vector of estimated SDs for group 1
#' @param   sd2    	 vector of estimated SDs for group 2
#' @param   n1     	 vector of group 1 sample sizes
#' @param   n2     	 vector of group 2 sample sizes
#' @param   bystudy  logical to also return each study estimate (TRUE) or not
#'
#'
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy 
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns: 
#'  * Estimate - estimated effect size 
#'  * SE - standard error 
#'  * LL - lower limit of the confidence interval
#'  * UL - upper limit of the confidence interval 
#'  * exp(Estimate) - exponentiated estimate 
#'  * exp(LL) - lower limit of the exponentiated confidence interval
#'  * exp(UL) - upper limit of the exponentiated confidence interval
#'  * df - degrees of freedom
#'
#'
#' @references 
#' \insertRef{Bonett2020}{vcmeta}
#'
#'
#' @examples
#' m1 <- c(7.4, 6.9)
#' m2 <- c(6.3, 5.7)
#' sd1 <- c(1.7, 1.5)
#' sd2 <- c(2.3, 2.0)
#' n1 <- c(40, 20)
#' n2 <- c(40, 20)
#' meta.ave.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
#'
#' # Should return:
#' #          Estimate         SE          LL        UL  exp(Estimate)
#' # Average 0.1759928 0.05738065 0.061437186 0.2905484       1.192429
#' # Study 1 0.1609304 0.06820167 0.024749712 0.2971110       1.174603
#' # Study 2 0.1910552 0.09229675 0.002986265 0.3791242       1.210526
#' #          exp(LL)  exp(UL)    df
#' # Average 1.063364 1.337161 66.27
#' # Study 1 1.025059 1.345965 65.70
#' # Study 2 1.002991 1.461004 31.71
#'
#'
#' @importFrom stats qt
#' @export
meta.ave.meanratio2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE) {
  m <- length(m1)
  v <- rep(1, m)*(1/m)
  logratio <- log(m1/m2)
  var1 <- sd1^2/(n1*m1^2) 
  var2 <- sd2^2/(n2*m2^2)
  est <- t(v)%*%logratio
  se <- sqrt(t(v)%*%(diag(var1 + var2))%*%v)
  df <- se^4/sum(v^4*var1^2/(n1 - 1) + v^4*var2^2/(n2 - 1))
  df <- round(df, 2)
  t <- qt(1 - alpha/2, df)
  ll <- est - t*se
  ul <- est + t*se
  out <- cbind(est, se, ll, ul, exp(est), exp(ll), exp(ul), df)
  row <- "Average"
  if (bystudy) {
    se <- sqrt(var1 + var2)
    u1 <- se^4
    u2 <- var1^2/(n1 - 1) + var2^2/(n2 - 1)
    df <- u1/u2
    t <- qt(1 - alpha/2, df)
    ll <- logratio - t*se
    ul <- logratio + t*se
    row2 <- t(t(paste(rep("Study", m), seq(1, m))))
    row <- rbind(row, row2)
    out2 <- cbind(logratio, se, ll, ul, exp(logratio), exp(ll), exp(ul), df)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)", "df")
  rownames(out) <- row 
  return(out)
}


#  meta.ave.meanratio.ps ==========================================================
#' Confidence interval for an average mean ratio from paired-samples studies
#'
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' geometric average mean ratio from two or more paired-samples studies. A 
#' Satterthwaite adjustment to the degrees of freedom is used to improve the 
#' accuracy of the confidence interval for the average effect size. Equality 
#' of variances within or across studies is not assumed.
#'
#'
#' @param   alpha		 alpha level for 1-alpha confidence
#' @param   m1		   vector of estimated means for measurement 1
#' @param   m2		   vector of estimated means for measurement 2
#' @param   sd1		   vector of estimated SDs for measurement 1
#' @param   sd2		   vector of estimated SDs for measurement 2
#' @param   cor		   vector of estimated correlations for paired measurements
#' @param   n		     vector of sample sizes
#' @param   bystudy  logical to also return each study estimate (TRUE) or not
#'
#'
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'  * exp(Estimate) - exponentiated estimate 
#'  * exp(LL) - lower limit of the exponentiated confidence interval
#'  * exp(UL) - upper limit of the exponentiated confidence interval
#'  * df - degrees of freedom
#'
#'
#' @references 
#' \insertRef{Bonett2020}{vcmeta}
#'
#'
#' @examples
#' m1 <- c(53, 60, 53, 57)
#' m2 <- c(55, 62, 58, 61)
#' sd1 <- c(4.1, 4.2, 4.5, 4.0)
#' sd2 <- c(4.2, 4.7, 4.9, 4.8)
#' cor <- c(.7, .7, .8, .85)
#' n <- c(30, 50, 30, 70)
#' meta.ave.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
#'
#' # Should return:
#' #            Estimate          SE          LL          UL
#' # Average -0.05695120 0.004350863 -0.06558008 -0.04832231
#' # Study 1 -0.03704127 0.010871086 -0.05927514 -0.01480740
#' # Study 2 -0.03278982 0.008021952 -0.04891054 -0.01666911
#' # Study 3 -0.09015110 0.009779919 -0.11015328 -0.07014892
#' # Study 4 -0.06782260 0.004970015 -0.07773750 -0.05790769
#' #         exp(Estimate)   exp(LL)   exp(UL)     df
#' # Average     0.9446402 0.9365240 0.9528266 103.03
#' # Study 1     0.9636364 0.9424474 0.9853017  29.00
#' # Study 2     0.9677419 0.9522663 0.9834691  49.00
#' # Study 3     0.9137931 0.8956968 0.9322550  29.00
#' # Study 4     0.9344262 0.9252073 0.9437371  69.00
#' 
#'
#' @importFrom stats qt
#'@export
meta.ave.meanratio.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, bystudy = TRUE) {
  m <- length(m1)
  v <- rep(1, m)*(1/m)
  logratio <- log(m1/m2)
  var <- (sd1^2/m1^2 + sd2^2/m2^2 - 2*cor*sd1*sd2/(m1*m2))/n
  est <- t(v)%*%logratio
  se <- sqrt(t(v)%*%(diag(var))%*%v)
  df <- se^4/sum(v^4*var^2/(n - 1))
  df <- round(df, 2)
  t <- qt(1 - alpha/2, df)
  ll <- est - t*se
  ul <- est + t*se
  out <- cbind(est, se, ll, ul, exp(est), exp(ll), exp(ul), df)
  row <- "Average"
  if (bystudy) {
    se <- sqrt(var)
    df <- n - 1
    t <- qt(1 - alpha/2, df)
    ll <- logratio - t*se
    ul <- logratio + t*se
    row2 <- t(t(paste(rep("Study", m), seq(1, m))))
    row <- rbind(row, row2)
    out2 <- cbind(logratio, se, ll, ul, exp(logratio), exp(ll), exp(ul), df)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)", "df")
  rownames(out) <- row 
  return(out)
}


#  meta.ave.cor ==========================================================
#' Confidence interval for an average Pearson or partial correlation
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average Pearson or partial correlation from two or more studies. The 
#' sample correlations must be all Pearson correlations or all partial
#' correlations. Use the meta.ave.cor.gen function to meta-analyze any 
#' combination of Pearson, partial, or Spearman correlations.
#' 
#' 
#' @param alpha	   alpha level for 1-alpha confidence
#' @param n     	 vector of sample sizes 
#' @param cor   	 vector of estimated correlations 
#' @param s     	 number of control variables (set to 0 for Pearson)
#' @param bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n <- c(55, 190, 65, 35)
#' cor <- c(.40, .65, .60, .45)
#' meta.ave.cor(.05, n, cor, 0, bystudy = TRUE)
#' 
#' # Should return:
#' #         Estimate      SE     LL     UL
#' # Average    0.525 0.05113 0.4177 0.6179
#' # Study 1    0.400 0.11431 0.1507 0.6015
#' # Study 2    0.650 0.04201 0.5594 0.7252
#' # Study 3    0.600 0.08000 0.4171 0.7362
#' # Study 4    0.450 0.13677 0.1374 0.6811
#' 
#' 
#' @references
#' \insertRef{Bonett2008a}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.cor <- function(alpha, n, cor, s, bystudy = TRUE) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  var.cor <- (1 - cor^2)^2/ (n - 3 - s)
  ave.cor <- sum(cor)/m
  se.ave <- sqrt(sum(var.cor)/m^2)
  z.ave <- log((1 + ave.cor)/(1 - ave.cor))/2
  ll0 <- z.ave - z*se.ave/(1 - ave.cor^2)
  ul0 <- z.ave + z*se.ave/(1 - ave.cor^2)
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- cbind(round(ave.cor, 4), round(se.ave, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    se.cor <- sqrt((1 - cor^2)^2/ (n - 1 - s))
    se.z <- sqrt(1/(n - 3 - s))
    z.cor <- log((1 + cor)/(1 - cor))/2
    ll0 <- z.cor - z*se.z
    ul0 <- z.cor + z*se.z
    ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
    ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(cor, 4), round(se.cor, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.ave.slope ==========================================================
#' Confidence interval for an average slope coefficient
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average slope coefficient in a simple linear regression model from two
#' or more studies. A Satterthwaite adjustment to the degrees of freedom
#' is used to improve the accuracy of the confidence interval.
#' 
#' 
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    n     	vector of sample sizes 
#' @param    cor   	vector of estimated correlations 
#' @param    sdy   	vector of estimated SDs of y
#' @param    sdx   	vector of estimated SDs of x
#' @param bystudy   logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @examples
#' n <- c(45, 85, 50, 60)
#' cor <- c(.24, .35, .16, .20)
#' sdy <- c(12.2, 14.1, 11.7, 15.9)
#' sdx <- c(1.34, 1.87, 2.02, 2.37)
#' meta.ave.slope(.05, n, cor, sdy, sdx, bystudy = TRUE)
#'
#' # Should return:
#' #          Estimate        SE         LL       UL     df
#' # Average 1.7731542 0.4755417  0.8335021 2.712806 149.48
#' # Study 1 2.1850746 1.3084468 -0.4536599 4.823809  43.00
#' # Study 2 2.6390374 0.7262491  1.1945573 4.083518  83.00
#' # Study 3 0.9267327 0.8146126 -0.7111558 2.564621  48.00
#' # Study 4 1.3417722 0.8456799 -0.3510401 3.034584  58.00
#' 
#' 
#' @importFrom stats qt
#' @export
meta.ave.slope <- function(alpha, n, cor, sdy, sdx, bystudy = TRUE) {
  m <- length(n)
  b <- cor*(sdy/sdx)
  var.b <- (sdy^2*(1 - cor^2)^2*(n - 1))/(sdx^2*(n - 1)*(n - 2))
  ave.b <- sum(b)/m
  se.ave <- sqrt(sum(var.b)/m^2)
  u1 <- sum(var.b)^2
  u2 <- sum(var.b^2/(n - 1))
  df <- round(u1/u2, 2)
  t <- qt(1 - alpha/2, df)
  ll <- ave.b - t*se.ave
  ul <- ave.b + t*se.ave
  out <- cbind(ave.b, se.ave, ll, ul, df)
  row <- "Average"
  if (bystudy) {
    se.b <- sqrt(var.b)
    df <- n - 2
    t <- qt(1 - alpha/2, df)
    ll <- b - t*se.b
    ul <- b + t*se.b
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(b, se.b, ll, ul, df)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) <- row
  return (out)
}


#  meta.ave.path ==========================================================
#' Confidence interval for an average slope coefficient in a general 
#' linear model or a path model. 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average slope coefficient in a general linear model (ANOVA, ANCOVA,
#' multiple regression) or a path model from two or more studies.
#'
#'
#' @param alpha    alpha level for 1-alpha confidence
#' @param n        vector of sample sizes 
#' @param slope    vector of slope estimates 
#' @param se       vector of slope standard errors
#' @param s        number of predictors of the response variable
#' @param bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy 
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'    
#'    
#' @examples  
#' n <- c(75, 85, 250, 160)
#' slope <- c(1.57, 1.38, 1.08, 1.25)
#' se <- c(.658, .724, .307, .493)
#' meta.ave.path(.05, n, slope, se, 2, bystudy = TRUE)
#'
#' #  Should return:
#' #          Estimate        SE          LL       UL     df
#' #  Average     1.32 0.2844334  0.75994528 1.880055 263.18
#' #  Study 1     1.57 0.6580000  0.25830097 2.881699  72.00
#' #  Study 2     1.38 0.7240000 -0.06026664 2.820267  82.00
#' #  Study 3     1.08 0.3070000  0.47532827 1.684672 247.00
#' #  Study 4     1.25 0.4930000  0.27623174 2.223768 157.00
#'
#'
#' @importFrom stats qt
#' @export
meta.ave.path <- function(alpha, n, slope, se, s, bystudy = TRUE) {
  m <- length(n)
  var.b <- se^2
  ave.b <- sum(slope)/m
  se.ave <- sqrt(sum(var.b)/m^2)
  u1 <- sum(var.b)^2
  u2 <- sum(var.b^2/(n - s - 1))
  df <- round(u1/u2, 2)
  t <- qt(1 - alpha/2, df)
  ll <- ave.b - t*se.ave
  ul <- ave.b + t*se.ave
  out <- cbind(ave.b, se.ave, ll, ul, df)
  row <- "Average"
  if (bystudy) {
    df <- n - s - 1
    t <- qt(1 - alpha/2, df)
    ll <- slope - t*se
    ul <- slope + t*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(slope, se, ll, ul, df)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) <- row
  return (out)
}


#  meta.ave.spear ==========================================================
#' Confidence interval for an average Spearman correlation 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average Spearman correlation from two or more studies. The Spearman 
#' correlation is preferred to the Pearson correlation if the relation 
#' between the two quantitative variables is monotonic rather than linear
#' or if the bivariate normality assumption is not plausible.
#'
#'
#' @param    alpha	  alpha level for 1-alpha confidence
#' @param    n     	  vector of sample sizes 
#' @param    cor   	  vector of estimated Spearman correlations 
#' @param    bystudy  logical to also return each study estimate (TRUE) or not
#' 
#'   
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n <- c(150, 200, 300, 200, 350)
#' cor <- c(.14, .29, .16, .21, .23)
#' meta.ave.spear(.05, n, cor, bystudy = TRUE)
#'
#' # Should return:
#' #         Estimate      SE      LL     UL
#' # Average    0.206 0.02944  0.1476 0.2629
#' # Study 1    0.140 0.08071 -0.0215 0.2944
#' # Study 2    0.290 0.06628  0.1548 0.4146
#' # Study 3    0.160 0.05671  0.0469 0.2691
#' # Study 4    0.210 0.06850  0.0719 0.3402
#' # Study 5    0.230 0.05136  0.1269 0.3282
#' 
#' 
#' @references
#' \insertRef{Bonett2008a}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.spear <- function(alpha, n, cor, bystudy = TRUE) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  var.cor <- (1 + cor^2/2)*(1 - cor^2)^2/(n - 3)
  ave.cor <- sum(cor)/m
  se.ave <- sqrt(sum(var.cor)/m^2)
  z.ave <- log((1 + ave.cor)/(1 - ave.cor))/2
  ll0 <- z.ave - z*se.ave/(1 - ave.cor^2)
  ul0 <- z.ave + z*se.ave/(1 - ave.cor^2)
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- cbind(round(ave.cor, 4), round(se.ave, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    se.cor <- sqrt((1 + cor^2/2)*(1 - cor^2)^2/(n - 1))
    se.z <- sqrt((1 + cor^2/2)/(n - 3))
    z.cor <- log((1 + cor)/(1 - cor))/2
    ll0 <- z.cor - z*se.z
    ul0 <- z.cor + z*se.z
    ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
    ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(cor, 4), round(se.cor, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.ave.pbcor ==========================================================
#' Confidence interval for an average point-biserial correlation
#'
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average point-biserial correlation from two or more studies. Two types
#' of point-biserial correlations can be meta-analyzed. One type uses
#' an unweighted variance and is appropriate in 2-group experimental
#' designs. The other type uses a weighted variance and is appropriate in
#' 2-group nonexperimental designs with simple random sampling (but not
#' stratified random sampling) within each study. This function requires 
#' all point-biserial correlations to be of the same type.  Use the
#' meta.ave.cor.gen function to meta-analyze any combination of biserial
#' correlation types. 
#'
#' 
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    m1	  	vector of estimated means for group 1 
#' @param    m2     vector of estimated means for group 2 
#' @param    sd1		vector of estimated SDs for group 1
#' @param    sd2		vector of estimated SDs for group 2
#' @param    n1		  vector of group 1 sample sizes
#' @param    n2		  vector of group 2 sample sizes
#' @param    type		
#' * set to 1 for weighted variance
#' * set to 2 for unweighted variance
#' @param bystudy   logical to also return each study estimate (TRUE) or not
#' 
#'    
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' m1 <- c(21.9, 23.1, 19.8)
#' m2 <- c(16.1, 17.4, 15.0)
#' sd1 <- c(3.82, 3.95, 3.67)
#' sd2 <- c(3.21, 3.30, 3.02)
#' n1 <- c(40, 30, 24)
#' n2 <- c(40, 28, 25)
#' meta.ave.pbcor(.05, m1, m2, sd1, sd2, n1, n2, 2, bystudy = TRUE)
#'
#' # Should return:
#' #         Estimate      SE     LL     UL
#' # Average   0.6109 0.04532 0.5144 0.6921
#' # Study 1   0.6350 0.06317 0.4842 0.7370
#' # Study 2   0.6165 0.07774 0.4260 0.7385
#' # Study 3   0.5812 0.09192 0.3551 0.7236
#' 
#' 
#' @references
#' \insertRef{Bonett2020b}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.pbcor <- function(alpha, m1, m2, sd1, sd2, n1, n2, type, bystudy = TRUE) {
  m <- length(m1)
  z <- qnorm(1 - alpha/2)
  n <- n1 + n2
  df1 <- n1 - 1
  df2 <- n2 - 1
  if (type == 1) {
    p <- n1/n
    b <- (n - 2)/(n*p*(1 - p))
    s <- sqrt((df1*sd1^2 + df2*sd2^2)/(df1 + df2))
    d <- (m1 - m2)/s
    se.d <- sqrt(d^2*(1/df1 + 1/df2)/8 + 1/n1 + 1/n2)
    se <- sqrt((b^2*se.d^2)/(d^2 + b)^3)
    cor <- d/sqrt(d^2 + b)
  }
  else {
    s <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s
    a1 <- d^2*(sd1^4/df1 + sd1^4/df2)/(8*s^4)
    a2 <- sd1^2/(s^2*df1) + sd2^2/(s^2*df2)
    se.d <- sqrt(a1 + a2)
    se <- sqrt((16*se.d^2)/(d^2 + 4)^3)
    cor <- d/sqrt(d^2 + 4)
  }
  ll.d <- d - z*se.d
  ul.d <- d + z*se.d
  ave <- sum(cor)/m
  var.ave <- sum(se^2)/m^2
  se.ave <- sqrt(var.ave)
  cor.f <- log((1 + ave)/(1 - ave))/2
  ll0 <- cor.f - z*se.ave/(1 - ave^2)
  ul0 <- cor.f + z*se.ave/(1 - ave^2)
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- cbind(round(ave, 4), round(se.ave, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    if (type == 1) {
      ll <- ll.d/sqrt(ll.d^2 + b)
      ul <- ul.d/sqrt(ul.d^2 + b)
    }
    else {    
      ll <- ll.d/sqrt(ll.d^2 + 4)
      ul <- ul.d/sqrt(ul.d^2 + 4)
    }
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(cor, 4), round(se, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return(out)
}


#  meta.ave.semipart ==========================================================
#' Confidence interval for an average semipartial correlation 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average semipartial correlation from two or more studies. 
#'
#'
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    n     	vector of sample sizes 
#' @param    cor   	vector of estimated semipartial correlations 
#' @param    r2  	  vector of squared multiple correlations for a model that
#' includes the IV and all control variables
#' @param bystudy   logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n <- c(128, 97, 210, 217)
#' cor <- c(.35, .41, .44, .39)
#' r2 <- c(.29, .33, .36, .39)
#' meta.ave.semipart(.05, n, cor, r2, bystudy = TRUE)
#'
#' # Should return:
#' #         Estimate      SE     LL     UL
#' # Average   0.3975 0.03221 0.3326 0.4587
#' # Study 1   0.3500 0.07175 0.2023 0.4821
#' # Study 2   0.4100 0.07886 0.2447 0.5521
#' # Study 3   0.4400 0.05147 0.3338 0.5351
#' # Study 4   0.3900 0.05085 0.2860 0.4849
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.semipart <- function(alpha, n, cor, r2, bystudy = TRUE) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  r0 <- r2 - cor^2
  var.cor <- (r2^2 - 2*r2 + r0 - r0^2 + 1)/(n - 3)
  ave.cor <- sum(cor)/m
  se.ave <- sqrt(sum(var.cor)/m^2)
  z.ave <- log((1 + ave.cor)/(1 - ave.cor))/2
  ll0 <- z.ave - z*se.ave/(1 - ave.cor^2)
  ul0 <- z.ave + z*se.ave/(1 - ave.cor^2)
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- cbind(round(ave.cor, 4), round(se.ave, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    se.cor = sqrt(var.cor)
    se.z <- se.cor/(1 - cor^2)
    z.cor <- log((1 + cor)/(1 - cor))/2
    ll0 <- z.cor - z*se.z
    ul0 <- z.cor + z*se.z
    ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
    ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(cor, 4), round(se.cor, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.ave.cronbach ==========================================================
#' Confidence interval for an average Cronbach alpha reliability 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average Cronbach reliability coefficient from two or more studies. 
#'
#'
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    n     	vector of sample sizes 
#' @param    rel   	vector of sample reliabilities 
#' @param    r     	number of measurements (e.g., items) used to compute each reliability
#' @param bystudy   logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n <- c(583, 470, 546, 680)
#' rel <- c(.91, .89, .90, .89)
#' meta.ave.cronbach(.05, n, rel, 10, bystudy = TRUE)
#'
#' # Should return:
#' #         Estimate      SE     LL     UL
#' # Average   0.8975 0.00326 0.8911 0.9039
#' # Study 1   0.9100 0.00557 0.8986 0.9204
#' # Study 2   0.8900 0.00758 0.8744 0.9041
#' # Study 3   0.9000 0.00639 0.8869 0.9119
#' # Study 4   0.8900 0.00630 0.8771 0.9018
#' 
#' 
#' @references
#' * \insertRef{Bonett2010}{vcmeta}
#' * \insertRef{Bonett2015b}{vcmeta}
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.cronbach <- function(alpha, n, rel, r, bystudy = TRUE) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  hn <- m/sum(1/n)
  a <- ((r - 2)*(m - 1))^.25
  var.rel <- 2*r*(1 - rel)^2/((r - 1)*(n - 2 - a))
  ave.rel <- sum(rel)/m
  se.ave <- sqrt(sum(var.rel)/m^2)
  log.ave <- log(1 - ave.rel) - log(hn/(hn - 1))
  ul <- 1 - exp(log.ave - z*se.ave/(1 - ave.rel))
  ll <- 1 - exp(log.ave + z*se.ave/(1 - ave.rel))
  out <- cbind(round(ave.rel, 4), round(se.ave, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    se.rel <- sqrt(2*r*(1 - rel)^2/((r - 1)*(n - 2)))
    log.rel <- log(1 - rel) - log(n/(n - 1))
    ul <- 1 - exp(log.rel - z*se.rel/(1 - rel))
    ll <- 1 - exp(log.rel + z*se.rel/(1 - rel))
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(rel, 4), round(se.rel, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.ave.oddsratio =====================================================
#' Confidence interval for average odds ratio from 2-group studies
#'  
#'  
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' geometric average odds ratio from two or more studies. 
#'
#'
#' @param    alpha  	alpha level for 1-alpha confidence
#' @param    f1     	vector of group 1 frequency counts 
#' @param    f2     	vector of group 2 frequency counts
#' @param    n1     	vector of group 1 sample sizes 
#' @param    n2     	vector of group 2 sample sizes
#' @param    bystudy  logical to also return each study estimate (TRUE) or not
#'
#'  
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
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
#' meta.ave.oddsratio(.05, f1, f2, n1, n2, bystudy = TRUE)
#' 
#' # Should return:
#' #           Estimate        SE          LL        UL 
#' # Average 0.86211102 0.2512852  0.36960107 1.3546210
#' # Study 1 0.02581353 0.3700520 -0.69947512 0.7511022
#' # Study 2 0.91410487 0.3830515  0.16333766 1.6648721
#' # Study 3 0.41496672 0.2226089 -0.02133877 0.8512722
#' # Study 4 1.52717529 0.6090858  0.33338907 2.7209615
#' # Study 5 1.42849472 0.9350931 -0.40425414 3.2612436
#' #         exp(Estimate)   exp(LL)   exp(UL)
#' # Average      2.368155 1.4471572  3.875292
#' # Study 1      1.026150 0.4968460  2.119335
#' # Study 2      2.494541 1.1774342  5.284997
#' # Study 3      1.514320 0.9788873  2.342625
#' # Study 4      4.605150 1.3956902 15.194925
#' # Study 5      4.172414 0.6674745 26.081952
#'
#'
#' @references 
#' \insertRef{Bonett2015}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.oddsratio <- function(alpha, f1, f2, n1, n2, bystudy = TRUE) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  lor <- log((f1 + .5)*(n2 - f2 + .5)/((f2 + .5)*(n1 - f1 + .5)))
  var.lor <- 1/(f1 + .5) + 1/(f2 + .5) + 1/(n1 - f1 + .5) + 1/(n2 - f2 + .5)
  ave.lor <- sum(lor)/m
  se.ave <- sqrt(sum(var.lor)/m^2)
  ll <- ave.lor - z*se.ave
  ul <- ave.lor + z*se.ave
  out <- cbind(ave.lor, se.ave, ll, ul, exp(ave.lor), exp(ll), exp(ul))
  row <- "Average"
  if (bystudy) {
    se <- sqrt(var.lor)
    ll <- lor - z*se
    ul <- lor + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(lor, se, ll, ul, exp(lor), exp(ll), exp(ul))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- row
  return (out)
}


#  meta.ave.propratio2 ==========================================================
#' Confidence interval for an average proportion ratio from 2-group studies 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' geometric average proportion ratio from two or more studies. 
#'
#'
#' @param    alpha  	alpha level for 1-alpha confidence
#' @param    f1     	vector of group 1 frequency counts
#' @param    f2     	vector of group 2 frequency counts
#' @param    n1     	vector of group 1 sample sizes 
#' @param    n2     	vector of group 2 sample sizes
#' @param    bystudy  logical to also return each study estimate (TRUE) or not
#'    
#'    
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * exp(Estimate) - exponentiated estimate
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' n1 <- c(204, 201, 932, 130, 77)
#' n2 <- c(106, 103, 415, 132, 83)
#' f1 <- c(24, 40, 93, 14, 5)
#' f2 <- c(12, 9, 28, 3, 1)
#' meta.ave.propratio2(.05, f1, f2, n1, n2, bystudy = TRUE)
#' 
#' #  Should return:
#' #           Estimate        SE          LL        UL 
#' # Average 0.84705608 0.2528742  0.35143178 1.3426804
#' # Study 1 0.03604257 0.3297404 -0.61023681 0.6823220
#' # Study 2 0.81008932 0.3442007  0.13546839 1.4847103
#' # Study 3 0.38746839 0.2065227 -0.01730864 0.7922454
#' # Study 4 1.49316811 0.6023296  0.31262374 2.6737125
#' # Study 5 1.50851199 0.9828420 -0.41782290 3.4348469
#' #     exp(Estimate)   exp(LL)   exp(UL)
#' # Average  2.332769 1.4211008  3.829294
#' # Study 1  1.036700 0.5432222  1.978466
#' # Study 2  2.248109 1.1450730  4.413686
#' # Study 3  1.473246 0.9828403  2.208350
#' # Study 4  4.451175 1.3670071 14.493677
#' # Study 5  4.520000 0.6584788 31.026662
#' 
#' 
#' @references 
#' \insertRef{Price2008}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.propratio2 <- function(alpha, f1, f2, n1, n2, bystudy = TRUE) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  p1 <- (f1 + 1/4)/(n1 + 7/4) 
  p2 <- (f2 + 1/4)/(n2 + 7/4)
  lrr <- log(p1/p2)
  v1 <- 1/(f1 + 1/4 + (f1 + 1/4)^2/(n1 - f1 + 3/2))
  v2 <- 1/(f2 + 1/4 + (f2 + 1/4)^2/(n2 - f2 + 3/2))
  var.lrr <- v1 + v2
  ave.lrr <- sum(lrr)/m
  se.ave <- sqrt(sum(var.lrr)/m^2)
  ll <- ave.lrr - z*se.ave
  ul <- ave.lrr + z*se.ave
  out <- cbind(ave.lrr, se.ave, ll, ul, exp(ave.lrr), exp(ll), exp(ul))
  row <- "Average"
  if (bystudy) {
    se <- sqrt(var.lrr)
    ll <- lrr - z*se
    ul <- lrr + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(lrr, se, ll, ul, exp(lrr), exp(ll), exp(ul))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- row
  return (out)
}


#  meta.ave.prop2 ==========================================================
#' Confidence interval for an average proportion difference in 
#' 2-group studies 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average proportion difference from two or more studies. 
#'
#'
#' @param    alpha  	alpha level for 1-alpha confidence
#' @param    f1     	vector of group 1 frequency counts
#' @param    f2     	vector of group 2 frequency counts
#' @param    n1     	vector of group 1 sample sizes 
#' @param    n2     	vector of group 2 sample sizes
#' @param bystudy     logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n1 <- c(204, 201, 932, 130, 77)
#' n2 <- c(106, 103, 415, 132, 83)
#' f1 <- c(24, 40, 93, 14, 5)
#' f2 <- c(12, 9, 28, 3, 1)
#' meta.ave.prop2(.05, f1, f2, n1, n2, bystudy = TRUE)
#' 
#' # Should return:
#' #             Estimate         SE            LL         UL
#' # Average 0.0567907589 0.01441216  2.854345e-02 0.08503807
#' # Study 1 0.0009888529 0.03870413 -7.486985e-02 0.07684756
#' # Study 2 0.1067323481 0.04018243  2.797623e-02 0.18548847
#' # Study 3 0.0310980338 0.01587717 -2.064379e-05 0.06221671
#' # Study 4 0.0837856174 0.03129171  2.245499e-02 0.14511624
#' # Study 5 0.0524199553 0.03403926 -1.429577e-02 0.11913568
#' 
#' 
#' @references
#' \insertRef{Bonett2014}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.prop2 <- function(alpha, f1, f2, n1, n2, bystudy = TRUE) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  p1 <- (f1 + 1/m)/(n1 + 2/m) 
  p2 <- (f2 + 1/m)/(n2 + 2/m)
  rd <- p1 - p2
  v1 <- p1*(1 - p1)/(n1 + 2/m)
  v2 <- p2*(1 - p2)/(n2 + 2/m) 
  var.rd <- v1 + v2
  ave.rd <- sum(rd)/m
  se.ave <- sqrt(sum(var.rd)/m^2)
  ll <- ave.rd - z*se.ave
  ul <- ave.rd + z*se.ave
  out <- cbind(ave.rd, se.ave, ll, ul)
  row <- "Average"
  if (bystudy) {
    p1 <- (f1 + 1)/(n1 + 2)   
    p2 <- (f2 + 1)/(n2 + 2)
    rd <- p1 - p2
    v1 <- p1*(1 - p1)/(n1 + 2)
    v2 <- p2*(1 - p2)/(n2 + 2) 
    se <- sqrt(v1 + v2)
    ll <- rd - z*se
    ul <- rd + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(rd, se, ll, ul)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row 
  return (out)
}


#  meta.ave.prop.ps ==========================================================
#' Confidence interval for an average proportion difference in 
#' paired-samples studies  
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average proportion difference from two or more studies. 
#'
#'
#' @param    alpha  	alpha level for 1-alpha confidence
#' @param    f11    	vector of frequency counts in cell 1,1
#' @param    f12    	vector of frequency counts in cell 1,2
#' @param    f21    	vector of frequency counts in cell 2,1
#' @param    f22    	vector of frequency counts in cell 2,2
#' @param bystudy     logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' f11 <- c(17, 28, 19)
#' f12 <- c(43, 56, 49)
#' f21 <- c(3, 5, 5)
#' f22 <- c(37, 54, 39)
#' meta.ave.prop.ps(.05, f11, f12, f21, f22, bystudy = TRUE)
#' 
#' # Should return:
#' #          Estimate         SE        LL        UL
#' # Average 0.3809573 0.03000016 0.3221581 0.4397565
#' # Study 1 0.3921569 0.05573055 0.2829270 0.5013867
#' # Study 2 0.3517241 0.04629537 0.2609869 0.4424614
#' # Study 3 0.3859649 0.05479300 0.2785726 0.4933572
#' 
#' 
#' @references 
#' \insertRef{Bonett2012}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.prop.ps <- function(alpha, f11, f12, f21, f22, bystudy = TRUE) {
  m <- length(f11)
  z <- qnorm(1 - alpha/2)
  n <- f11 + f12 + f21 + f22
  p12 <- (f12 + 1/m)/(n + 2/m) 
  p21 <- (f21 + 1/m)/(n + 2/m)
  rd <- p12 - p21
  var.rd <- (p12 + p21 - rd^2)/(n + 2/m)
  ave.rd <- sum(rd)/m
  se.ave <- sqrt(sum(var.rd)/m^2)
  ll <- ave.rd - z*se.ave
  ul <- ave.rd + z*se.ave
  out <- cbind(ave.rd, se.ave, ll, ul)
  row <- "Average"
  if (bystudy) {
    p12 <- (f12 + 1)/(n + 2) 
    p21 <- (f21 + 1)/(n + 2)
    rd <- p12 - p21
    se <- sqrt((p12 + p21 - rd^2)/(n + 2))
    ll <- rd - z*se
    ul <- rd + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(rd, se, ll, ul)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.ave.agree ==========================================================
#' Confidence interval for an average G-index agreement coefficient 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average G-index of agreement from two or more studies. This function 
#' assumes that two raters each provide a dichotomous rating to a sample
#' of objects. As a measure of agreement, the G-index is usually preferred
#' to Cohen's kappa.
#'
#'
#' @param    alpha  	alpha level for 1-alpha confidence
#' @param    f11    	vector of frequency counts in cell 1,1
#' @param    f12    	vector of frequency counts in cell 1,2
#' @param    f21    	vector of frequency counts in cell 2,1
#' @param    f22    	vector of frequency counts in cell 2,2
#' @param    bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' f11 <- c(43, 56, 49)
#' f12 <- c(7, 2, 9)
#' f21 <- c(3, 5, 5)
#' f22 <- c(37, 54, 39)
#' meta.ave.agree(.05, f11, f12, f21, f22, bystudy = TRUE)
#' 
#' # Should return:
#' #        Estimate      SE     LL     UL
#' # Average  0.7843 0.03540 0.7149 0.8537
#' # Study 1  0.7447 0.06884 0.6098 0.8796
#' # Study 2  0.8512 0.04771 0.7577 0.9447
#' # Study 3  0.6981 0.06954 0.5618 0.8344
#' 
#' 
#' @references 
#' \insertRef{Bonett2022}{vcmeta}
#'
#'
#' @importFrom stats qnorm
#' @export
meta.ave.agree <- function(alpha, f11, f12, f21, f22, bystudy = TRUE) {
  m <- length(f11)
  z <- qnorm(1 - alpha/2)
  n <- f11 + f12 + f21 + f22
  p0 <- (f11 + f22 + 2/m)/(n + 4/m)
  g <- 2*p0 - 1 
  ave.g <- sum(g)/m
  var.g <- 4*p0*(1 - p0)/(n + 4/m)
  se.ave <- sqrt(sum(var.g)/m^2)
  ll <- ave.g - z*se.ave
  ul <- ave.g + z*se.ave
  out <- cbind(round(ave.g, 4), round(se.ave, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    p0 <- (f11 + f22 + 2)/(n + 4)
    g <- 2*p0 - 1 
    se <- sqrt(4*p0*(1 - p0)/(n + 4))
    ll <- g - z*se 
    ul <- g + z*se 
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(g, 4), round(se, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return (out)
}


# meta.ave.var ==========================================================
#' Confidence interval for an average variance
#'
#'
#' @description
#' Computes the estimate and confidence interval for an average variance 
#' from two or more studies. The estimated average variance or the upper
#' confidence limit could be used as a variance planning value in sample
#' size planning. The confidence intervals assume normality, and this
#' function is not recommended if the variances have been estimated
#' from leptokurtic data.
#'
#'  
#' @param alpha  	 alpha level for 1-alpha confidence
#' @param var   	 vector of sample variances 
#' @param n     	 vector of sample sizes
#' @param bystudy  logical to also return each study estimate (TRUE) or not
#'
#'
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated variance
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' var <- c(26.63, 22.45, 34.12)
#' n <- c(40, 30, 50)
#' meta.ave.var(.05, var, n, bystudy = TRUE)
#'
#' # Should return:
#' #         Estimate       LL       UL
#' # Average 27.73333 21.45679 35.84589
#' # Study 1 26.63000 17.86939 43.90614
#' # Study 2 22.45000 14.23923 40.57127
#' # Study 3 34.12000 23.80835 52.98319
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats qchisq
#' @export
meta.ave.var <- function(alpha, var, n, bystudy = TRUE) {
 m <- length(n)
 z <- qnorm(1 - alpha/2)
 var.var <- 2*var^2/(n - 1)
 ave.var <- sum(var)/m
 se.ave <- sqrt(sum(var.var)/m^2)
 ln.ave <- log(ave.var)
 ll <- exp(ln.ave - z*se.ave/ave.var)
 ul <- exp(ln.ave + z*se.ave/ave.var)
 out <- cbind(ave.var, ll, ul)
 row <- "Average"
 if (bystudy) {
   chi.U <- qchisq(1 - alpha/2, (n - 1))
   ll <- (n - 1)*var/chi.U
   chi.L <- qchisq(alpha/2, (n - 1))
   ul <- (n - 1)*var/chi.L
   row2 <- t(t(paste(rep("Study", m), seq(1,m))))
   row <- rbind(row, row2)
   out2 <- cbind(var, ll, ul)
   out <- rbind(out, out2)
 }
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- row
 return (out)
}


#  meta.ave.gen ==========================================================
#' Confidence interval for an average of any parameter 
#' 
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average of any type of parameter from two or more studies. 
#'
#' 
#' @param    alpha  	alpha level for 1-alpha confidence
#' @param    est     	vector of parameter estimates 
#' @param    se      	vector of standard errors
#' @param    bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' 
#' est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
#' se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
#' meta.ave.gen(.05, est, se, bystudy = TRUE)
#' 
#' # Should return:
#' #          Estimate        SE          LL        UL
#' # Average  0.393125 0.1561622  0.08705266 0.6991973
#' # Study 1  0.022000 0.1240000 -0.22103553 0.2650355
#' # Study 2  0.751000 0.4640000 -0.15842329 1.6604233
#' # Study 3  0.421000 0.1020000  0.22108367 0.6209163
#' # Study 4  0.287000 0.5920000 -0.87329868 1.4472987
#' # Study 5  0.052000 0.8640000 -1.64140888 1.7454089
#' # Study 6  0.146000 0.2410000 -0.32635132 0.6183513
#' # Study 7  0.562000 0.2520000  0.06808908 1.0559109
#' # Study 8  0.904000 0.3180000  0.28073145 1.5272685
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.gen <- function(alpha, est, se, bystudy = TRUE) {
  m <- length(est)
  z <- qnorm(1 - alpha/2)
  ave.est <- sum(est)/m
  se.ave <- sqrt(sum(se^2)/m^2)
  ll <- ave.est - z*se.ave
  ul <- ave.est + z*se.ave
  out <- cbind(ave.est, se.ave, ll, ul)
  row <- "Average"
  if (bystudy) {
    ll <- est - z*se
    ul <- est + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(est, se, ll, ul)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row 
  return (out)
}


#  meta.ave.gen.cc ==========================================================
#' Confidence interval for an average effect size using a constant 
#' coefficient model 
#'
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' weighted average effect from two or more studies using the constant
#' coefficient (fixed-effect) meta-analysis model. 
#'
#'
#' @details
#' When the population effect sizes are not identical across studies, the
#' weighted average estimate will be biased regardless of the number of 
#' studies or the sample size in each study. The actual confidence interval 
#' coverage probability can be much smaller than the specified confidence
#' level when the population effect sizes are not identical across studies. 
#'
#' The constant coefficient model should be used with caution, and the varying
#' coefficient methods in this package are the recommended alternatives. The 
#' varying coefficient methods do not require effect-size homogeneity across 
#' the selected studies. This constant coefficient meta-analysis function is 
#' included in the vcmeta package primarily for classroom  demonstrations to 
#' illustrate the problematic characteristics of the constant coefficient 
#' meta-analysis model.
#'    
#' 
#' @param    alpha   	  alpha level for 1-alpha confidence
#' @param    est     	  vector of parameter estimates 
#' @param    se      	  vector of standard errors
#' @param    bystudy    logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
#' se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
#' meta.ave.gen.cc(.05, est, se, bystudy = TRUE)
#' 
#' # Should return:
#' #             Estimate         SE          LL        UL
#' # Average    0.3127916 0.06854394  0.17844794 0.4471352
#' # Study 1    0.0220000 0.12400000 -0.22103553 0.2650355
#' # Study 2    0.7510000 0.46400000 -0.15842329 1.6604233
#' # Study 3    0.4210000 0.10200000  0.22108367 0.6209163
#' # Study 4    0.2870000 0.59200000 -0.87329868 1.4472987
#' # Study 5    0.0520000 0.86400000 -1.64140888 1.7454089
#' # Study 6    0.1460000 0.24100000 -0.32635132 0.6183513
#' # Study 7    0.5620000 0.25200000  0.06808908 1.0559109
#' # Study 8    0.9040000 0.31800000  0.28073145 1.5272685
#' 
#' 
#' @references
#' \insertRef{Hedges1985}{vcmeta}
#'
#' \insertRef{Borenstein2009}{vcmeta}
#'
#'
#' @seealso \link[vcmeta]{meta.ave.gen}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.gen.cc <- function(alpha, est, se, bystudy = TRUE) {
  m <- length(est)
  z <- qnorm(1 - alpha/2)
  w <- 1/se^2
  weighted.est <- sum(w*est)/sum(w)
  se.w <- sqrt(1/sum(w))
  ll <- weighted.est - z*se.w
  ul <- weighted.est + z*se.w
  out <- cbind(weighted.est, se.w, ll, ul)
  row <- "Average"
  if (bystudy) {
    ll <- est - z*se
    ul <- est + z*se
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(est, se, ll, ul)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row 
  return (out)
}


#  meta.ave.gen.rc ==========================================================
#' Confidence interval for an average effect size using a random coefficient 
#' model 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' weighted average effect from multiple studies using the random 
#' coefficient (random-effects) meta-analysis model. An estimate of 
#' effect-size heterogeneity (tau-squared) is also computed. 
#'
#'
#' @details
#' The random coefficient model assumes that the studies in the meta-analysis
#' are a random sample from some definable superpopulation of studies. This
#' assumption is very difficult to justify. When the effect sizes are 
#' correlated with the weights (a common occurance), the weighted average 
#' estimate will be biased regardless of the number of studies or the sample
#' size in each study. The actual confidence interval coverage probability  
#' can much smaller than the specified confidence level if the effect sizes 
#' are correlated with the weights. 
#'
#' The confidence interval for tau-squared assumes that the true effect sizes
#' in the superpopulation of studies have a normal distribution. A large number 
#' of studies, each with a large sample size, is required to assess the 
#' superpopulation normality assumption and to accurately estimate 
#' tau-squared. The confidence interval for the population tau-squared is
#' hypersensitive to very minor and difficult-to-detect violations of the 
#' superpopulation normality assumption. 
#'
#' The random coefficient model should be used with caution, and the varying 
#' coefficient methods in this package are the recommended alternatives. The 
#' varying coefficient methods allows the effect sizes to differ across studies,
#' but do not require the studies to be a random sample from a definable 
#' superpopoulation of studies. This random coefficient function is included 
#' in the vcmeta package primarily for classroom demonstrations to illustrate
#' the problimatic characteristics of the random coefficient meta-analysis
#' model.
#' 
#' 
#' @param    alpha    alpha level for 1-alpha confidence
#' @param    est      vector of parameter estimates 
#' @param    se       vector of standard errors
#' @param    bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is true, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
#' se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
#' meta.ave.gen.rc(.05, est, se, bystudy = TRUE)
#' 
#' # Should return: 
#' #               Estimate        SE          LL        UL
#' # Tau-squared 0.03772628 0.0518109  0.00000000 0.1392738
#' # Average     0.35394806 0.1155239  0.12752528 0.5803708
#' # Study 1     0.02200000 0.1240000 -0.22103553 0.2650355
#' # Study 2     0.75100000 0.4640000 -0.15842329 1.6604233
#' # Study 3     0.42100000 0.1020000  0.22108367 0.6209163
#' # Study 4     0.28700000 0.5920000 -0.87329868 1.4472987
#' # Study 5     0.05200000 0.8640000 -1.64140888 1.7454089
#' # Study 6     0.14600000 0.2410000 -0.32635132 0.6183513
#' # Study 7     0.56200000 0.2520000  0.06808908 1.0559109
#' # Study 8     0.90400000 0.3180000  0.28073145 1.5272685
#' 
#' @references
#' \insertRef{Hedges1985}{vcmeta}
#'
#' \insertRef{Borenstein2009}{vcmeta}
#'
#'
#' @seealso \link[vcmeta]{meta.ave.gen}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @export
meta.ave.gen.rc <- function(alpha, est, se, bystudy = TRUE) {
 m <- length(est)
 z <- qnorm(1 - alpha/2)
 w1 <- 1/se^2
 sw1 <- sum(w1)
 sw2 <- sum(w1*w1)
 sw3 <- sum(w1*w1*w1)
 C <- sw1 - sw2/sw1
 w1.est <- sum(w1*est)/sw1
 Q <- sum(w1*(est - w1.est)*(est - w1.est))
 v <- Q - m + 1
 t2 <- v/C 
 if (t2 < 0) t2 = 0
 A <- (m - 1 + 2*C*t2 + (sw2 - 2*(sw3/sw1) + sw2^2/sw1^2)*t2^2)
 se.t2 <- sqrt(2*A/C^2)
 w2 <- 1/(se^2 + t2)
 w2.est <- sum(w2*est)/sum(w2)
 se.w2 <- sqrt(1/sum(w2))
 ll.t2 <- t2 - z*se.t2
 ul.t2 <- t2 + z*se.t2
 if (ll.t2 < 0) {ll.t2 = 0}
 ll <- w2.est - z*se.w2
 ul <- w2.est + z*se.w2
 out1 <- cbind(t2, se.t2, ll.t2, ul.t2)
 out2 <- cbind(w2.est, se.w2, ll, ul)
 out <- rbind(out1, out2)
 row1 <- "Tau-squared"
 row2 <- "Average"
 row <- rbind(row1, row2)
 if (bystudy) {
   ll <- est - z*se
   ul <- est + z*se
   row3 <- t(t(paste(rep("Study", m), seq(1,m))))
   row <- rbind(row1, row2, row3)
   out3 <- cbind(est, se, ll, ul)
   out <- rbind(out1, out2, out3)
 }
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- row 
 return (out)
}


#  meta.ave.cor.gen ==========================================================
#' Confidence interval for an average correlation of any type
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average correlation. Any type of correlation can be used (e.g., Pearson,
#' Spearman, semipartial, factor correlation, gamma coefficient, Somers d
#' coefficient, tetrachoric, point-biserial, biserial, correlation between
#' latent factors, etc.).
#' 
#' 
#' @param alpha	   alpha level for 1-alpha confidence
#' @param cor      vector of estimated correlations 
#' @param se   	   vector of standard errors 
#' @param bystudy  logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated effect size
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' cor <- c(.396, .454, .409, .502, .350)
#' se <- c(.104, .064, .058, .107, .086)
#' meta.ave.cor.gen(.05, cor, se, bystudy = TRUE)
#' 
#' # Should return:
#' #         Estimate      SE     LL     UL
#' # Average   0.4222 0.03853 0.3439 0.4947
#' # Study 1   0.3960 0.10400 0.1753 0.5788
#' # Study 2   0.4540 0.06400 0.3201 0.5701
#' # Study 3   0.4090 0.05800 0.2894 0.5160
#' # Study 4   0.5020 0.10700 0.2651 0.6817
#' # Study 5   0.3500 0.08600 0.1716 0.5061
#' 
#' 
#' @references
#' \insertRef{Bonett2008a}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.cor.gen <- function(alpha, cor, se, bystudy = TRUE) {
  m <- length(cor)
  z <- qnorm(1 - alpha/2)
  ave.cor <- sum(cor)/m
  se.ave <- sqrt(sum(se^2)/m^2)
  z.ave <- log((1 + ave.cor)/(1 - ave.cor))/2
  ll0 <- z.ave - z*se.ave/(1 - ave.cor^2)
  ul0 <- z.ave + z*se.ave/(1 - ave.cor^2)
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- cbind(round(ave.cor, 4), round(se.ave, 5), round(ll, 4), round(ul, 4))
  row <- "Average"
  if (bystudy) {
    se.z <- se/(1 - cor^2)
    z.cor <- log((1 + cor)/(1 - cor))/2
    ll0 <- z.cor - z*se.z
    ul0 <- z.cor + z*se.z
    ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
    ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(round(cor, 4), round(se, 5), round(ll, 4), round(ul, 4))
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- row
  return (out)
}


#  meta.ave.gen.log ==========================================================
#' Exponentiated confidence interval for an average of log-transformed
#' parameters 
#'                        
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' average of any type of log-transformed parameter (e.g., log mean ratio,
#' log proportion ratio, log odds ratio) from two or more studies. 
#'
#' 
#' @param    alpha  	alpha level for 1-alpha confidence
#' @param    est     	vector of log-transformed parameter estimates 
#' @param    se      	vector of standard errors for log-transformed estimates
#' @param    bystudy    logical to also return each study estimate (TRUE) or not
#' 
#' 
#' @return 
#' Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
#' is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
#' * Estimate - estimated log effect size (from input)
#' * SE - standard error of log effect size (from input)
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * exp(Estimate) - exponentiated estimate
#' * LL - lower limit of the exponentiated confidence interval
#' * UL - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' 
#' est <- c(.165, .193, .218)
#' se <- c(.0684, .0921, .0882)
#' meta.ave.gen.log(.05, est, se, bystudy = TRUE)
#' 
#' # Should return:
#' #         Estimate         SE         LL        UL exp(Estimate)
#' # Average    0.192 0.04823578 0.09745962 0.2865404      1.211671
#' # Study 1    0.165 0.06840000 0.03093846 0.2990615      1.179393
#' # Study 2    0.193 0.09210000 0.01248732 0.3735127      1.212883
#' # Study 3    0.218 0.08820000 0.04513118 0.3908688      1.243587
#' #          exp(LL)  exp(UL)
#' # Average 1.102367 1.331812
#' # Study 1 1.031422 1.348593
#' # Study 2 1.012566 1.452829
#' # Study 3 1.046165 1.478265
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.ave.gen.log <- function(alpha, est, se, bystudy = TRUE) {
  m <- length(est)
  z <- qnorm(1 - alpha/2)
  ave.est <- sum(est)/m
  se.ave <- sqrt(sum(se^2)/m^2)
  ll1 <- ave.est - z*se.ave
  ul1 <- ave.est + z*se.ave
  ll2 <- exp(ave.est - z*se.ave)
  ul2 <- exp(ave.est + z*se.ave)
  out <- cbind(ave.est, se.ave, ll1, ul1, exp(ave.est), ll2, ul2)
  row <- "Average"
  if (bystudy) {
    ll1 <- est - z*se
    ul1 <- est + z*se
    ll2 <- exp(est - z*se)
    ul2 <- exp(est + z*se)
    row2 <- t(t(paste(rep("Study", m), seq(1,m))))
    row <- rbind(row, row2)
    out2 <- cbind(est, se, ll1, ul1, exp(est), ll2, ul2)
    out <- rbind(out, out2)
  }
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- row 
  return (out)
}


use_imports <- function() {
  mathjaxr::preview_rd()
}
