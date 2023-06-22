#  replicate.mean2 ============================================================
#' Compares and combines 2-group mean differences in original and follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals for a 2-group mean difference
#' from an original study and a follow-up study. Confidence intervals for the
#' difference and average effect size also are computed. A Satterthwaite
#' adjustment to the degrees of freedom is used to improve the accuracy of the
#' confidence intervals. The same results can be obtained using the \link[vcmeta]{meta.lc.mean2}
#' function with appropriate contrast coefficients. The confidence level for the
#' difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m11		   estimated mean for group 1 in original study 
#' @param    m12		   estimated mean for group 2 in original study
#' @param    sd11   	 estimated SD for group 1 in original study
#' @param    sd12   	 estimated SD for group 2 in original study
#' @param    n11    	 sample size for group 1 in original study
#' @param    n12    	 sample size for group 2 in original study
#' @param    m21    	 estimated mean for group 1 in follow-up study 
#' @param    m22    	 estimated mean for group 2 in follow-up study
#' @param    sd21   	 estimated SD for group 1 in follow-up study
#' @param    sd22   	 estimated SD for group 2 in follow-up study
#' @param    n21    	 sample size for group 1 in follow-up study
#' @param    n22    	 sample size for group 2 in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#'
#'
#' The columns are:
#' * Estimate - mean difference estimate (single study, difference, average)
#' * SE - standard error
#' * t - t-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#'    
#' 
#' @examples
#' replicate.mean2(.05, 21.9, 16.1, 3.82, 3.21, 40, 40, 
#'                      25.2, 19.1, 3.98, 3.79, 75, 75)
#'
#' # Should return:
#' #                     Estimate        SE          t            p 
#' # Original:               5.80 0.7889312  7.3517180 1.927969e-10  
#' # Follow-up:              6.10 0.6346075  9.6122408 0.000000e+00  
#' # Original - Follow-up:  -0.30 1.0124916 -0.2962988 7.673654e-01 
#' # Average:                5.95 0.5062458 11.7531843 0.000000e+00 
#' #                               LL       UL        df
#' # Original:               4.228624 7.371376  75.75255
#' # Follow-up:              4.845913 7.354087 147.64728
#' # Original - Follow-up:  -1.974571 1.374571 169.16137
#' # Average:                4.950627 6.949373 169.16137
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
replicate.mean2 <- function(alpha, m11, m12, sd11, sd12, n11, n12, m21, m22, sd21, sd22, n21, n22){
  v11 <- sd11^2; v12 <- sd12^2 
  v21 <- sd21^2; v22 <- sd22^2
  est1 <- m11 - m12
  est2 <- m21 - m22
  est3 <- est1 - est2
  est4 <- (est1 + est2)/2
  se1 <- sqrt(v11/n11 + v12/n12)
  se2 <- sqrt(v21/n21 + v22/n22)
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  v1 <- v11^2/(n11^3 - n11^2)
  v2 <- v12^2/(n12^3 - n12^2)
  v3 <- v21^2/(n21^3 - n21^2)
  v4 <- v22^2/(n22^3 - n22^2)
  df1 <- (se1^4)/(v1 + v2)
  df2 <- (se2^4)/(v3 + v4)
  df3 <- (se3^4)/(v1 + v2 + v3 + v4)
  t1 <- est1/se1
  t2 <- est2/se2
  t3 <- est3/se3
  t4 <- est4/se4
  pval1 <- 2*(1 - pt(abs(t1),df1))
  pval2 <- 2*(1 - pt(abs(t2),df2)) 
  pval3 <- 2*(1 - pt(abs(t3),df3)) 
  pval4 <- 2*(1 - pt(abs(t4),df3))
  tcrit1 <- qt(1 - alpha/2, df1)
  tcrit2 <- qt(1 - alpha/2, df2)
  tcrit3 <- qt(1 - alpha, df3)
  tcrit4 <- qt(1 - alpha/2, df3)
  ll1 <- est1 - tcrit1*se1;  ul1 <- est1 + tcrit1*se1
  ll2 <- est2 - tcrit2*se2;  ul2 <- est2 + tcrit2*se2
  ll3 <- est3 - tcrit3*se3;  ul3 <- est3 + tcrit3*se3
  ll4 <- est4 - tcrit4*se4;  ul4 <- est4 + tcrit4*se4
  out1 <- t(c(est1, se1, t1, pval1, ll1, ul1, df1))
  out2 <- t(c(est2, se2, t2, pval2, ll2, ul2, df2))
  out3 <- t(c(est3, se3, t3, pval3, ll3, ul3, df3))
  out4 <- t(c(est4, se4, t4, pval4, ll4, ul4, df3))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.mean.ps ============================================================
#' Compares and combines paired-samples mean differences in original and follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals for a paired-samples mean difference
#' from an original study and a follow-up study. Confidence intervals for the
#' difference and average effect size also are computed. A Satterthwaite
#' adjustment to the degrees of freedom is used to improve the accuracy of the
#' confidence intervals for the difference and average. The same results can be 
#' obtained using the \link[vcmeta]{meta.lc.mean.ps} function with appropriate 
#' contrast coefficients. The confidence level for the difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m11		   estimated mean for group 1 in original study 
#' @param    m12		   estimated mean for group 2 in original study
#' @param    sd11   	 estimated SD for group 1 in original study
#' @param    sd12   	 estimated SD for group 2 in original study
#' @param    n1    	   sample size in original study
#' @param    cor1    	 estimated correlation of paired observations in orginal study
#' @param    m21    	 estimated mean for group 1 in follow-up study 
#' @param    m22    	 estimated mean for group 2 in follow-up study
#' @param    sd21   	 estimated SD for group 1 in follow-up study
#' @param    sd22   	 estimated SD for group 2 in follow-up study
#' @param    n2    	   sample size in follow-up study
#' @param    cor2    	 estimated correlation of paired observations in follow-up study
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#' 
#' 
#' The columns are:
#' * Estimate - mean difference estimate (single study, difference, average)
#' * SE - standard error
#' * df - degrees of freedom
#' * t - t-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.mean.ps(.05, 86.22, 70.93, 14.89, 12.32, .765, 20, 
#'                        84.81, 77.24, 15.68, 16.95, .702, 75)
#'
#' #  Should return:
#' #                     Estimate       SE        t            p  
#' # Original:              15.29 2.154344 7.097288 9.457592e-07 
#' # Follow-up:              7.57 1.460664 5.182575 1.831197e-06 
#' # Original - Follow-up:   7.72 2.602832 2.966000 5.166213e-03 
#' # Average:               11.43 1.301416 8.782740 1.010232e-10 
#' #                              LL       UL       df
#' # Original:             10.780906 19.79909 19.00000
#' # Follow-up:             4.659564 10.48044 74.00000
#' # Original - Follow-up:  3.332885 12.10712 38.40002
#' # Average:               8.796322 14.06368 38.40002
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
replicate.mean.ps <- function(alpha, m11, m12, sd11, sd12, cor1, n1, m21, m22, sd21, sd22, cor2, n2) {
  v11 <- sd11^2;  v12 <- sd12^2
  v21 <- sd21^2;  v22 <- sd22^2
  vd1 <- v11 + v12 - 2*cor1*sd11*sd12
  vd2 <- v21 + v22 - 2*cor2*sd21*sd22
  est1 <- m11 - m12
  est2 <- m21 - m22
  est3 <- est1 - est2
  est4 <- (est1 + est2)/2
  se1 <- sqrt(vd1/n1)
  se2 <- sqrt(vd2/n2)
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  df1 <- n1 - 1
  df2 <- n2 - 1
  df3 <- se3^4/(se1^4/df1 + se2^4/df2)
  t1 <- est1/se1
  t2 <- est2/se2
  t3 <- est3/se3
  t4 <- est4/se4
  pval1 <- 2*(1 - pt(t1, df1))
  pval2 <- 2*(1 - pt(t2, df2))
  pval3 <- 2*(1 - pt(t3, df3))
  pval4 <- 2*(1 - pt(t4, df3))
  tcrit1 <- qt(1 - alpha/2, df1)
  tcrit2 <- qt(1 - alpha/2, df2)
  tcrit3 <- qt(1 - alpha, df3)
  tcrit4 <- qt(1 - alpha/2, df3)
  ll1 <- est1 - tcrit1*se1;  ul1 <- est1 + tcrit1*se1
  ll2 <- est2 - tcrit2*se2; ul2 <- est2 + tcrit2*se2
  ll3 <- est3 - tcrit3*se3; ul3 <- est3 + tcrit3*se3
  ll4 <- est4 - tcrit4*se4; ul4 <- est4 + tcrit4*se4
  out1 <- t(c(est1, se1, t1, pval1, ll1, ul1, df1))
  out2 <- t(c(est2, se2, t2, pval2, ll2, ul2, df2))
  out3 <- t(c(est3, se3, t3, pval3, ll3, ul3, df3))
  out4 <- t(c(est4, se4, t4, pval4, ll4, ul4, df3))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.stdmean2 ============================================================
#' Compares and combines 2-group standardized mean differences in original and 
#' follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals for a 2-group standardized mean 
#' difference from an original study and a follow-up study. Confidence intervals
#' for the difference and average effect size also are computed. The same results 
#' can be obtained using the \link[vcmeta]{meta.lc.stdmean2} function with 
#' appropriate contrast coefficients. The confidence level for the difference is
#' 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m11		   estimated mean for group 1 in original study 
#' @param    m12		   estimated mean for group 2 in original study
#' @param    sd11   	 estimated SD for group 1 in original study
#' @param    sd12   	 estimated SD for group 2 in original study
#' @param    n11    	 sample size for group 1 in original study
#' @param    n12    	 sample size for group 2 in original study
#' @param    m21    	 estimated mean for group 1 in follow-up study 
#' @param    m22    	 estimated mean for group 2 in follow-up study
#' @param    sd21   	 estimated SD for group 1 in follow-up study
#' @param    sd22   	 estimated SD for group 2 in follow-up study
#' @param    n21    	 sample size for group 1 in follow-up study
#' @param    n22    	 sample size for group 2 in follow-up study
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#' 
#' 
#' The columns are:
#' * Estimate - standardized mean difference estimate (single study, difference, average)
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.stdmean2(.05, 21.9, 16.1, 3.82, 3.21, 40, 40, 
#'                         25.2, 19.1, 3.98, 3.79, 75, 75)
#'
#' # Should return: 
#' #                           Estimate        SE         LL        UL
#' #  Original:              1.62803662 0.2594668  1.1353486 2.1524396
#' #  Follow-up:             1.56170447 0.1870576  1.2030461 1.9362986
#' #  Original - Follow-up:  0.07422178 0.3198649 -0.4519092 0.6003527
#' #  Average:               1.59487055 0.1599325  1.2814087 1.9083324
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
replicate.stdmean2 <- function(alpha, m11, m12, sd11, sd12, n11, n12, m21, m22, sd21, sd22, n21, n22) {
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  v11 <- sd11^2
  v12 <- sd12^2
  v21 <- sd21^2
  v22 <- sd22^2
  df11 <- n11 - 1
  df12 <- n12 - 1
  df21 <- n21 - 1
  df22 <- n22 - 1
  s1 <- sqrt((v11 + v12)/2)
  s2 <- sqrt((v21 + v22)/2)
  a1 <- 1 - 3/(4*(n11 + n12) - 9)
  a2 <- 1 - 3/(4*(n21 + n22) - 9)
  est1 <- (m11 - m12)/s1
  est2 <- (m21 - m22)/s2
  est3 <- est1 - est2
  est4 <- (a1*est1 + a2*est2)/2
  se1 <- sqrt(est1^2*(v11^2/df11 + v12^2/df12)/(8*s1^4) + (v11/df11 + v12/df12)/s1^2)
  se2 <- sqrt(est2^2*(v21^2/df21 + v22^2/df22)/(8*s2^4) + (v21/df21 + v22/df22)/s2^2)
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  ll1 <- est1 - zcrit1*se1;  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2;  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3;  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4;  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(a1*est1, se1, ll1, ul1))
  out2 <- t(c(a2*est2, se2, ll2, ul2))
  out3 <- t(c(est3, se3, ll3, ul3))
  out4 <- t(c(est4, se4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.stdmean.ps ============================================================
#' Compares and combines paired-samples standardized mean differences in original and 
#' follow-up studies
#' 
#' 
#' @description 
#' This function computes confidence intervals for a paired-samples standardized
#' mean difference from an original study and a follow-up study. Confidence intervals
#' for the difference and average effect size also are computed. The same results can
#' be obtained using the \link[vcmeta]{meta.lc.stdmean.ps} function with appropriate 
#' contrast coefficients. The confidence level for the difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m11		   estimated mean for group 1 in original study 
#' @param    m12		   estimated mean for group 2 in original study
#' @param    sd11   	 estimated SD for group 1 in original study
#' @param    sd12   	 estimated SD for group 2 in original study
#' @param    cor1    	 estimated correlation of paired observations in orginal study
#' @param    n1        sample size in original study
#' @param    m21    	 estimated mean for group 1 in follow-up study 
#' @param    m22    	 estimated mean for group 2 in follow-up study
#' @param    sd21   	 estimated SD for group 1 in follow-up study
#' @param    sd22   	 estimated SD for group 2 in follow-up study
#' @param    cor2    	 estimated correlation of paired observations in follow-up study
#' @param    n2        sample size in follow-up study
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#' 
#' 
#' The columns are:
#' * Estimate - standardized mean difference estimate (single study, difference, average)
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.stdmean.ps(alpha = .05, 
#'   m11 = 86.22, m12 = 70.93, sd11 = 14.89, sd12 = 12.32, cor1 = .765, n1 = 20, 
#'   m21 = 84.81, m22 = 77.24, sd21 = 15.68, sd22 = 16.95, cor2 = .702, n2 = 75
#' )
#'
#' # Should return:
#' #                         Estimate         SE        LL        UL
#' #  Orginal:              1.0890300 0.22915553 0.6697353 1.5680085
#' #  Follow-up:            0.4604958 0.09590506 0.2756687 0.6516096
#' #  Original - Follow-up: 0.6552328 0.24841505 0.2466264 1.0638392
#' #  Average:              0.7747629 0.12420752 0.5313206 1.0182052
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
replicate.stdmean.ps <- function(alpha, m11, m12, sd11, sd12, cor1, n1, m21, m22, sd21, sd22, cor2, n2) {
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  v11 <- sd11^2
  v12 <- sd12^2
  v21 <- sd21^2
  v22 <- sd22^2
  df1 <- n1 - 1
  df2 <- n2 - 1
  s1 <- sqrt((v11 + v12)/2)
  s2 <- sqrt((v21 + v22)/2)
  vd1 <- v11 + v12 - 2*cor1*sd11*sd12
  vd2 <- v21 + v22 - 2*cor2*sd21*sd22
  a1 <- sqrt((n1 - 2)/df1)
  a2 <- sqrt((n2 - 2)/df2)
  est1 <- (m11 - m12)/s1
  est2 <- (m21 - m22)/s2
  est3 <- est1 - est2
  est4 <- (a1*est1 + a2*est2)/2
  se1 <- sqrt(est1^2*(v11^2 + v12^2 + 2*cor1^2*v11*v12)/(8*df1*s1^4) + vd1/(df1*s1^2))
  se2 <- sqrt(est2^2*(v21^2 + v22^2 + 2*cor2^2*v21*v22)/(8*df2*s2^4) + vd2/(df2*s2^2))
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  ll1 <- est1 - zcrit1*se1;  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2;  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3;  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4;  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(a1*est1, se1, ll1, ul1))
  out2 <- t(c(a2*est2, se2, ll2, ul2))
  out3 <- t(c(est3, se3, ll3, ul3))
  out4 <- t(c(est4, se4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Orginal:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.cor ============================================================
#' Compares and combines Pearson or partial correlations in original and follow-up studies
#' 
#'
#' @description 
#' This function can be used to compare and combine Pearson or partial 
#' correlations from the original study and the follow-up study. The 
#' confidence level for the difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha	 alpha level for 1-alpha confidence
#' @param    cor1  	 estimated Pearson correlation in original study
#' @param    n1    	 sample size in original study
#' @param    cor2  	 estimated Pearson correlation in follow-up study
#' @param    n2    	 sample size in follow-up study
#' @param    s     	 number of control variables in each study (0 for Pearson)
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#' 
#' 
#' The columns are:
#' * Estimate - Pearson or partial correlation estimate (single study, difference, average)
#' * SE - standard error
#' * z - t-value for rows 1 and 2; z-value for rows 3 and 4
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.cor(.05, .598, 80, .324, 200, 0)
#'
#' # Should return:
#' #                      Estimate         SE        z            p        LL        UL
#' # Original:                0.598 0.07320782 6.589418 4.708045e-09 0.4355043 0.7227538
#' # Follow-up:               0.324 0.06376782 4.819037 2.865955e-06 0.1939787 0.4428347
#' # Original - Follow-up:    0.274 0.09708614 2.633335 8.455096e-03 0.1065496 0.4265016
#' # Average:                 0.461 0.04854307 7.634998 2.264855e-14 0.3725367 0.5411607
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats pt
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
replicate.cor <- function(alpha, cor1, n1, cor2, n2, s) {
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  zr1 <- log((1 + cor1)/(1 - cor1))/2
  zr2 <- log((1 + cor2)/(1 - cor2))/2
  se1 <- sqrt((1 - cor1^2)^2/(n1 - 3 - s))
  se2 <- sqrt((1 - cor2^2)^2/(n2 - 3 - s))
  dif <- cor1 - cor2
  ave <- (cor1 + cor2)/2
  ave.z <- log((1 + ave)/(1 - ave))/2
  se1.z <- sqrt(1/((n1 - 3 - s)))
  se2.z <- sqrt(1/((n2 - 3 - s)))
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- sqrt(se1^2 + se2^2)/2
  se4.z <- sqrt(((se1^2 + se2^2)/4)/(1 - ave^2))
  t1 <- cor1*sqrt(n1 - 2)/sqrt(1 - cor1^2) 
  t2 <- cor2*sqrt(n2 - 2)/sqrt(1 - cor2^2) 
  t3 <- (zr1 - zr2)/sqrt(se1.z^2 + se2.z^2)
  t4 <- (zr1 + zr2)/sqrt(se1.z^2 + se2.z^2)
  pval1 <- 2*(1 - pt(abs(t1), n1 - 2 - s))
  pval2 <- 2*(1 - pt(abs(t2), n2 - 2 - s))
  pval3 <- 2*(1 - pnorm(abs(t3)))
  pval4 <- 2*(1 - pnorm(abs(t4)))
  ll0a <- zr1 - zcrit1*se1.z;  ul0a <- zr1 + zcrit1*se1.z
  ll1a <- (exp(2*ll0a) - 1)/(exp(2*ll0a) + 1)
  ul1a <- (exp(2*ul0a) - 1)/(exp(2*ul0a) + 1)
  ll0b <- zr1 - zcrit2*se1.z;  ul0b <- zr1 + zcrit2*se1.z
  ll1b <- (exp(2*ll0b) - 1)/(exp(2*ll0b) + 1)
  ul1b <- (exp(2*ul0b) - 1)/(exp(2*ul0b) + 1)
  ll0a <- zr2 - zcrit1*se2.z;  ul0a <- zr2 + zcrit1*se2.z
  ll2a <- (exp(2*ll0a) - 1)/(exp(2*ll0a) + 1)
  ul2a <- (exp(2*ul0a) - 1)/(exp(2*ul0a) + 1)
  ll0b <- zr2 - zcrit2*se2.z;  ul0b <- zr2 + zcrit2*se2.z
  ll2b <- (exp(2*ll0b) - 1)/(exp(2*ll0b) + 1)
  ul2b <- (exp(2*ul0b) - 1)/(exp(2*ul0b) + 1)
  ll3 <- dif - sqrt((cor1 - ll1b)^2 + (ul2b - cor2)^2)
  ul3 <- dif + sqrt((ul1b - cor1)^2 + (cor2 - ll2b)^2)
  ll0 <- ave.z - zcrit1*se4.z
  ul0 <- ave.z + zcrit1*se4.z
  ll4 <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul4 <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out1 <- t(c(cor1, se1, t1, pval1, ll1a, ul1a))
  out2 <- t(c(cor2, se2, t2, pval2, ll2a, ul2a))
  out3 <- t(c(dif, se3, t3, pval3, ll3, ul3))
  out4 <- t(c(ave, se4, t4, pval4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


# replicate.prop2 ============================================================
#' Compares and combines 2-group proportion differences in original and follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals for a 2-group proportion difference
#' from an original study and a follow-up study. Confidence intervals for the
#' difference and average effect size also are computed. The same results can be 
#' obtained using the \link[vcmeta]{meta.lc.prop2} function with appropriate 
#' contrast coefficients. The confidence level for the difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    f11		   frequency count for group 1 in original study 
#' @param    f12		   frequency count for group 2 in original study
#' @param    n11    	 sample size for group 1 in original study
#' @param    n12    	 sample size for group 2 in original study
#' @param    f21    	 frequency count for group 1 in follow-up study 
#' @param    f22    	 frequency count for group 2 in follow-up study
#' @param    n21    	 sample size for group 1 in follow-up study
#' @param    n22    	 sample size for group 2 in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#'
#'
#' The columns are:
#' * Estimate - proportion difference estimate (single study, difference, average)
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'    
#' 
#' @examples
#' replicate.prop2(.05, 21, 16, 40, 40, 19, 13, 60, 60)
#'
#' # Should return:
#' #                         Estimate         SE         z         p
#' # Original:             0.11904762 0.10805233 1.1017590 0.2705665
#' # Follow-up:            0.09677419 0.07965047 1.2149858 0.2243715
#' # Original - Follow-up: 0.02359056 0.13542107 0.1742016 0.8617070
#' # Average:              0.11015594 0.06771053 1.6268656 0.1037656
#' #                                LL        UL
#' # Original:             -0.09273105 0.3308263
#' # Follow-up:            -0.05933787 0.2528863
#' # Original - Follow-up: -0.19915727 0.2463384
#' # Average:              -0.02255427 0.2428661
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
replicate.prop2 <- function(alpha, f11, f12, n11, n12, f21, f22, n21, n22){
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  p11.o <- (f11 + 1)/(n11 + 2)
  p12.o <- (f12 + 1)/(n12 + 2)
  p21.f <- (f21 + 1)/(n21 + 2)
  p22.f <- (f22 + 1)/(n22 + 2)
  est1 <- p11.o - p12.o
  est2 <- p21.f - p22.f
  p11 <- (f11 + .5)/(n11 + 1)
  p12 <- (f12 + .5)/(n12 + 1)
  p21 <- (f21 + .5)/(n21 + 1)
  p22 <- (f22 + .5)/(n22 + 1)
  est3 <- (p11 - p12) - (p21 - p22)
  est4 <- ((p11 - p12) + (p21 - p22))/2
  v11 <- p11.o*(1 - p11.o)/(n11 + 2)
  v12 <- p12.o*(1 - p12.o)/(n12 + 2)
  v21 <- p21.f*(1 - p21.f)/(n21 + 2)
  v22 <- p22.f*(1 - p22.f)/(n22 + 2)
  se1 <- sqrt(v11 + v12)
  se2 <- sqrt(v21 + v22)
  v11 <- p11*(1 - p11)/(n11 + 1)
  v12 <- p12*(1 - p12)/(n12 + 1)
  v21 <- p21*(1 - p21)/(n21 + 1)
  v22 <- p22*(1 - p22)/(n22 + 1)
  se3 <- sqrt(v11 + v12 + v21 + v22)
  se4 <- se3/2
  z1 <- est1/se1
  z2 <- est2/se2
  z3 <- est3/se3
  z4 <- est4/se4
  p1 <- 2*(1 - pnorm(abs(z1)))
  p2 <- 2*(1 - pnorm(abs(z2)))
  p3 <- 2*(1 - pnorm(abs(z3))) 
  p4 <- 2*(1 - pnorm(abs(z4)))
  ll1 <- est1 - zcrit1*se1;  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2;  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3;  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4;  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(est1, se1, z1, p1, ll1, ul1))
  out2 <- t(c(est2, se2, z2, p2, ll2, ul2))
  out3 <- t(c(est3, se3, z3, p3, ll3, ul3))
  out4 <- t(c(est4, se4, z4, p4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


# replicate.oddsratio ============================================================
#' Compares and combines odds ratios in original and follow-up studies
#' 
#' @description 
#' This function computes confidence intervals for an odds ratio from an
#' original study and a follow-up study. Confidence intervals for the
#' ratio of odds ratios and average odds ratio size also are computed. 
#' The confidence level for the difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    est1		   estimate of log odds ratio in original study 
#' @param    se1		   standard error of log odds ratio in original study
#' @param    est2   	 estimate of log odds ratio in follow-up study 
#' @param    se2    	 standard error of log odds ratio in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the ratio of odds ratios between studies
#' * Row 4 estimates the average effect size for the two studies
#'
#'
#' The columns are:
#' * Estimate - odds ratio estimate (single study, difference, average)
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - exponentiated lower limit of the confidence interval
#' * UL - exponentiated upper limit of the confidence interval
#'    
#' 
#' @examples
#' replicate.oddsratio(.05, 1.39, .302, 1.48, .206)
#'
#' # Should return:
#' #                          Estimate        SE          z            p
#' # Original:              1.39000000 0.3020000  4.6026490 4.171509e-06
#' # Follow-up:             1.48000000 0.2060000  7.1844660 6.747936e-13
#' # Original - Follow-up: -0.06273834 0.3655681 -0.1716188 8.637372e-01
#' # Average:               0.36067292 0.1827840  1.9732190 4.847061e-02
#' #                         exp(LL)  exp(UL)
#' # Original:             2.2212961 7.256583
#' # Follow-up:            2.9336501 6.578144
#' # Original/Fllow-up:    0.5147653 1.713551
#' # Average:              1.0024257 2.052222
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
replicate.oddsratio <- function(alpha, est1, se1, est2, se2){
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  est3 <- log(est1) - log(est2)
  est4 <- (log(est1) + log(est2))/2
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  z1 <- est1/se1
  z2 <- est2/se2
  z3 <- est3/se3
  z4 <- est4/se4
  p1 <- 2*(1 - pnorm(abs(z1)))
  p2 <- 2*(1 - pnorm(abs(z2)))
  p3 <- 2*(1 - pnorm(abs(z3))) 
  p4 <- 2*(1 - pnorm(abs(z4)))
  ll1 <- exp(est1 - zcrit1*se1);  ul1 <- exp(est1 + zcrit1*se1)
  ll2 <- exp(est2 - zcrit1*se2);  ul2 <- exp(est2 + zcrit1*se2)
  ll3 <- exp(est3 - zcrit2*se3);  ul3 <- exp(est3 + zcrit2*se3)
  ll4 <- exp(est4 - zcrit1*se4);  ul4 <- exp(est4 + zcrit1*se4)
  out1 <- t(c(est1, se1, z1, p1, ll1, ul1))
  out2 <- t(c(est2, se2, z2, p2, ll2, ul2))
  out3 <- t(c(est3, se3, z3, p3, ll3, ul3))
  out4 <- t(c(est4, se4, z4, p4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "exp(LL)", "exp(UL)")
  rownames(out) <- c("Original:", "Follow-up:", "Original/Follow-up:", "Average:")
  return(out)
}


# replicate.slope ============================================================ 
#' Compares and combines slope coefficients in original and follow-up studies
#'
#' @description 
#' Computes confidence intervals for a slope in original and follow-up studies,
#' the difference in slopes, and the average of the slopes. Equality of error 
#' variances between studies is not assumed. The confidence interval for the
#' difference uses a 1 - 2alpha confidence level. Use the replicate.gen 
#' function for slopes in other types of models (e.g., binary logistic, ordinal 
#' logistic, SEM). 
#'
#'
#' @param    alpha	alpha level for 1-alpha or 1 - 2alpha confidence
#' @param    b1     sample slope in original study 
#' @param    se1    standard error of slope in original study
#' @param    n1     sample size in original study
#' @param    b2     sample slope in follow-up study
#' @param    se2    standard error of slope in follow-up study
#' @param    n2     sample size in follow-up study
#' @param    s      number of predictor variables in model
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#'
#'
#' The columns are:
#' * Estimate - slope estimate (single study, difference, average)
#' * SE - standard error
#' * t - t-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#'    
#' 
#' @examples
#' replicate.slope(.05, 23.4, 5.16, 50, 18.5, 4.48, 90, 4)
#'
#' # Should return: 
#' #                       Estimate       SE         t            p
#' # Original:                23.40 5.160000 4.5348837 4.250869e-05
#' # Follow-up:               18.50 4.480000 4.1294643 8.465891e-05
#' # Original - Follow-up:     4.90 6.833447 0.7170612 4.749075e-01
#' # Average:                 20.95 3.416724 6.1316052 1.504129e-08
#' #                              LL       UL       df
#' # Original:             13.007227 33.79277  45.0000
#' # Follow-up:             9.592560 27.40744  85.0000
#' # Original - Follow-up: -6.438743 16.23874 106.4035
#' # Average:              14.176310 27.72369 106.4035
#'
#'
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
replicate.slope <- function(alpha, b1, se1, n1, b2, se2, n2, s) {
 df1 <- n1 - s - 1
 df2 <- n2 - s - 1
 est1 <- b1
 est2 <- b2
 est3 <- est1 - est2
 est4 <- (est1 + est2)/2
 se3 <- sqrt(se1^2 + se2^2)
 se4 <- se3/2
 v1 <- se1^4/df1
 v2 <- se2^4/df2
 df3 <- (se3^4)/(v1 + v2)
 t1 <- est1/se1
 t2 <- est2/se2
 t3 <- est3/se3
 t4 <- est4/se4
 pval1 <- 2*(1 - pt(abs(t1),df1))
 pval2 <- 2*(1 - pt(abs(t2),df2)) 
 pval3 <- 2*(1 - pt(abs(t3),df3)) 
 pval4 <- 2*(1 - pt(abs(t4),df3))
 tcrit1 <- qt(1 - alpha/2, df1)
 tcrit2 <- qt(1 - alpha/2, df2)
 tcrit3 <- qt(1 - alpha, df3)
 tcrit4 <- qt(1 - alpha/2, df3)
 ll1 <- est1 - tcrit1*se1;  ul1 <- est1 + tcrit1*se1
 ll2 <- est2 - tcrit2*se2;  ul2 <- est2 + tcrit2*se2
 ll3 <- est3 - tcrit3*se3;  ul3 <- est3 + tcrit3*se3
 ll4 <- est4 - tcrit4*se4;  ul4 <- est4 + tcrit4*se4
 out1 <- t(c(est1, se1, t1, pval1, ll1, ul1, df1))
 out2 <- t(c(est2, se2, t2, pval2, ll2, ul2, df2))
 out3 <- t(c(est3, se3, t3, pval3, ll3, ul3, df3))
 out4 <- t(c(est4, se4, t4, pval4, ll4, ul4, df3))
 out <- rbind(out1, out2, out3, out4)
 colnames(out) <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
 rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
 return(out)
}


#  replicate.gen ============================================================
#' Compares and combines effect sizes in original and follow-up studies
#'
#'
#' @description 
#' This function can be used to compare and combine any effect size (e.g.,odds
#' ratio, proportion ratio, proportion difference, slope coefficient, etc.)
#' using the effect size estimate and its standard error from the original study
#' and the follow-up study. The same results can be obtained using the
#' \link[vcmeta]{meta.lc.gen} function with appropriate contrast coefficients. 
#' The confidence level for the difference is 1 – 2alpha.
#' 
#'
#' @param    alpha		 alpha level for 1-alpha confidence 
#' @param    est1  	   estimated effect size in original study
#' @param    se1    	 effect size standard error in original study
#' @param    est2  	   estimated effect size in follow-up study
#' @param    se2    	 effect size standard error in follow-up study
#'   
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#' 
#' 
#' Columns are:
#' * Estimate - effect size estimate (single study, difference, average)
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.gen(.05, .782, .210, .650, .154)
#'
#' # Should return: 
#' #                      Estimate        SE         z            p         LL        UL
#' #  Original:              0.782 0.2100000 3.7238095 1.962390e-04  0.3704076 1.1935924
#' #  Follow-up:             0.650 0.1540000 4.2207792 2.434593e-05  0.3481655 0.9518345
#' #  Original - Follow-up:  0.132 0.2604151 0.5068831 6.122368e-01 -0.2963446 0.5603446
#' #  Average:               0.716 0.1302075 5.4989141 3.821373e-08  0.4607979 0.9712021
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
replicate.gen <- function(alpha, est1, se1, est2, se2) {
  est3 <- est1 - est2
  est4 <- (est1 + est2)/2
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  z1 <- est1/se1
  z2 <- est2/se2
  z3 <- est3/se3
  z4 <- est4/se4
  pval1 <- 2*(1 - pnorm(abs(z1)))
  pval2 <- 2*(1 - pnorm(abs(z2))) 
  pval3 <- 2*(1 - pnorm(abs(z3))) 
  pval4 <- 2*(1 - pnorm(abs(z4)))
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  ll1 <- est1 - zcrit1*se1;  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2;  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3;  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4;  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(est1, se1, z1, pval1, ll1, ul1))
  out2 <- t(c(est2, se2, z2, pval2, ll2, ul2))
  out3 <- t(c(est3, se3, z3, pval3, ll3, ul3))
  out4 <- t(c(est4, se4, z4, pval4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


# replicate.spear ===============================================================
#' Compares and combines Spearman correlations in original and follow-up studies
#' 
#'                           
#' @description 
#' This function can be used to compare and combine Spearman correlations from
#' the original study and the follow-up study. The confidence level for the 
#' difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha	 alpha level for 1-alpha confidence
#' @param    cor1  	 estimated Spearman correlation in original study
#' @param    n1    	 sample size in original study
#' @param    cor2  	 estimated Spearman correlation in follow-up study
#' @param    n2    	 sample size in follow-up study
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average effect size for the two studies
#' 
#' 
#' The columns are:
#' * Estimate - Spearman correlation estimate (single study, difference, average)
#' * SE - standard error
#' * z - z-value
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.spear(.05, .598, 80, .324, 200)
#'
#' # Should return:
#' #                      Estimate         SE        z            p         LL        UL
#' # Original:                0.598 0.07948367 5.315140 1.065752e-07 0.41985966 0.7317733
#' # Follow-up:               0.324 0.06541994 4.570582 4.863705e-06 0.19049455 0.4457384
#' # Original - Follow-up:    0.274 0.10294378 3.437975 5.860809e-04 0.09481418 0.4342171
#' # Average:                 0.461 0.05147189 9.967944 0.000000e+00 0.36695230 0.5457190
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
replicate.spear <- function(alpha, cor1, n1, cor2, n2) {
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  dif <- cor1 - cor2
  ave <- (cor1 + cor2)/2
  ave.z <- log((1 + ave)/(1 - ave))/2
  zr1 <- log((1 + cor1)/(1 - cor1))/2
  zr2 <- log((1 + cor2)/(1 - cor2))/2
  se1 <- sqrt((1 + cor1^2/2)*(1 - cor1^2)^2/(n1 - 3))
  se2 <- sqrt((1 + cor2^2/2)*(1 - cor2^2)^2/(n2 - 3))
  se1.z <- sqrt((1 + cor1^2/2)/((n1 - 3)))
  se2.z <- sqrt((1 + cor2^2/2)/((n2 - 3)))
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- sqrt(se1^2 + se2^2)/2
  se4.z <- sqrt(((se1^2 + se2^2)/4)/(1 - ave^2))
  t1 <- cor1*sqrt(n1 - 1) 
  t2 <- cor2*sqrt(n2 - 1)
  t3 <- (zr1 - zr2)/sqrt(se1^2 + se2^2)
  t4 <- (zr1 + zr2)/sqrt(se1^2 + se2^2)
  pval1 <- 2*(1 - pnorm(abs(t1)))
  pval2 <- 2*(1 - pnorm(abs(t2)))
  pval3 <- 2*(1 - pnorm(abs(t3)))
  pval4 <- 2*(1 - pnorm(abs(t4)))
  ll0a <- zr1 - zcrit1*se1.z;  ul0a <- zr1 + zcrit1*se1.z
  ll1a <- (exp(2*ll0a) - 1)/(exp(2*ll0a) + 1)
  ul1a <- (exp(2*ul0a) - 1)/(exp(2*ul0a) + 1)
  ll0a <- zr2 - zcrit1*se2.z;  ul0a <- zr2 + zcrit1*se2.z
  ll2a <- (exp(2*ll0a) - 1)/(exp(2*ll0a) + 1)
  ul2a <- (exp(2*ul0a) - 1)/(exp(2*ul0a) + 1)
  ll0b <- zr1 - zcrit2*se1.z;  ul0b <- zr1 + zcrit2*se1.z
  ll1b <- (exp(2*ll0b) - 1)/(exp(2*ll0b) + 1)
  ul1b <- (exp(2*ul0b) - 1)/(exp(2*ul0b) + 1)
  ll0b <- zr2 - zcrit2*se2.z;  ul0b <- zr2 + zcrit2*se2.z
  ll2b <- (exp(2*ll0b) - 1)/(exp(2*ll0b) + 1)
  ul2b <- (exp(2*ul0b) - 1)/(exp(2*ul0b) + 1)
  ll3 <- dif - sqrt((cor1 - ll1b)^2 + (ul2b - cor2)^2)
  ul3 <- dif + sqrt((ul1b - cor1)^2 + (cor2 - ll2b)^2)
  ll0 <- ave.z - zcrit1*se4.z
  ul0 <- ave.z + zcrit1*se4.z
  ll4 <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul4 <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out1 <- t(c(cor1, se1, t1, pval1, ll1a, ul1a))
  out2 <- t(c(cor2, se2, t2, pval2, ll2a, ul2a))
  out3 <- t(c(dif, se3, t3, pval3, ll3, ul3))
  out4 <- t(c(ave, se4, t4, pval4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.prop1 ============================================================
#' Compares and combines single proportion in original and follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals for a single proportion from an 
#' original study and a follow-up study. Confidence intervals for the
#' difference between the two proportions and average of the two proportions 
#' also are computed. The same results can be obtained using the \link[vcmeta]{meta.lc.prop1}
#' function with appropriate contrast coefficients. The confidence level for the
#' difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    f1	  	   frequency count in original study 
#' @param    n1     	 sample size in original study
#' @param    f2 	   	 frequency count in follow-up study 
#' @param    n2    	 	 sample size for in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average mean for the two studies
#'
#'
#' The columns are:
#' * Estimate - proportion estimate (single study, difference, average)
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'    
#' 
#' @examples
#' replicate.prop1(.05, 21, 300, 35, 400)
#'
#' # Should return:
#' #                          Estimate         SE          LL         UL
#' # Original:              0.07565789 0.01516725  0.04593064 0.10538515
#' # Follow-up:             0.09158416 0.01435033  0.06345803 0.11971029
#' # Original - Follow-up: -0.01670456 0.02065098 -0.05067239 0.01726328
#' # Average:               0.08119996 0.01032549  0.06096237 0.10143755
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
replicate.prop1 <- function(alpha, f1, n1, f2, n2){
  est1 <- (f1 + 2)/(n1 + 4)
  est2 <- (f2 + 2)/(n2 + 4)
  est1.d <- (f1 + 1)/(n1 + 2) 
  est2.d <- (f2 + 1)/(n2 + 2)
  est3 <- est1.d - est2.d
  est4 <- (est1.d + est2.d)/2
  se1 <- sqrt(est1*(1 - est1)/(n1 + 4))
  se2 <- sqrt(est2*(1 - est2)/(n2 + 4))
  se3 <- sqrt(est1.d*(1 - est1.d)/(n1 + 2) + est2.d*(1 - est2.d)/(n2 + 2))
  se4 <- se3/2
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit3 <- qnorm(1 - alpha)
  ll1 <- est1 - zcrit1*se1;  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2;  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit3*se3;  ul3 <- est3 + zcrit3*se3
  ll4 <- est4 - zcrit1*se4;  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(est1, se1, ll1, ul1))
  out2 <- t(c(est2, se2, ll2, ul2))
  out3 <- t(c(est3, se3, ll3, ul3))
  out4 <- t(c(est4, se4, ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.mean1 ============================================================
#' Compares and combines single mean in original and follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals for a single mean from an 
#' original study and a follow-up study. Confidence intervals for the
#' difference between the two means and average of the two means also are
#' computed. A Satterthwaite adjustment to the degrees of freedom is used to
#' improve the accuracy of the confidence intervals for the difference and
#' average. The same results can be obtained using the \link[vcmeta]{meta.lc.mean1}
#' function with appropriate contrast coefficients. The confidence level for the
#' difference is 1 – 2alpha.
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m1	  	   estimated mean in original study 
#' @param    sd1   	   estimated SD in original study
#' @param    n1     	 sample size in original study
#' @param    m2 	   	 estimated mean in follow-up study 
#' @param    sd2   		 estimated SD in follow-up study
#' @param    n2    	 	 sample size for in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference between studies
#' * Row 4 estimates the average mean for the two studies
#'
#'
#' The columns are:
#' * Estimate - mean estimate (single study, difference, average)
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#'    
#' 
#' @examples
#' replicate.mean1(.05, 21.9, 3.82, 40, 25.2, 3.98, 75)
#'
#' # Should return:
#' #                       Estimate        SE        LL        UL       df
#' # Original:                21.90 0.6039950 20.678305 23.121695 39.00000
#' # Follow-up:               25.20 0.4595708 24.284285 26.115715 74.00000
#' # Original - Follow-up:    -3.30 0.7589567 -4.562527 -2.037473 82.63282
#' # Average:                 23.55 0.3794784 22.795183 24.304817 82.63282
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @export
replicate.mean1 <- function(alpha, m1, sd1, n1, m2, sd2, n2){
  v1 <- sd1^2
  v2 <- sd2^2
  est1 <- m1
  est2 <- m2
  est3 <- est1 - est2
  est4 <- (est1 + est2)/2
  se1 <- sqrt(v1/n1)
  se2 <- sqrt(v2/n2)
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  v1 <- v1^2/(n1^3 - n1^2)
  v2 <- v2^2/(n2^3 - n2^2)
  df1 <- n1 - 1
  df2 <- n2 - 1
  df3 <- (se3^4)/(v1 + v2)
  tcrit1 <- qt(1 - alpha/2, df1)
  tcrit2 <- qt(1 - alpha/2, df2)
  tcrit3 <- qt(1 - alpha, df3)
  tcrit4 <- qt(1 - alpha/2, df3)
  ll1 <- est1 - tcrit1*se1;  ul1 <- est1 + tcrit1*se1
  ll2 <- est2 - tcrit2*se2;  ul2 <- est2 + tcrit2*se2
  ll3 <- est3 - tcrit3*se3;  ul3 <- est3 + tcrit3*se3
  ll4 <- est4 - tcrit4*se4;  ul4 <- est4 + tcrit4*se4
  out1 <- t(c(est1, se1, ll1, ul1, df1))
  out2 <- t(c(est2, se2, ll2, ul2, df2))
  out3 <- t(c(est3, se3, ll3, ul3, df3))
  out4 <- t(c(est4, se4, ll4, ul4, df3))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}
