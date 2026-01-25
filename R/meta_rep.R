#  replicate.mean2 ============================================================
#' Compares and combines 2-group mean differences in original and 
#' follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals from an original study and a 
#' follow-up study where the effect size is a 2-group mean difference. Confidence
#' intervals for the difference and average effect size are also computed. 
#' Equality of variances within or across studies is not assumed. A
#' Satterthwaite adjustment to the degrees of freedom is used to improve the 
#' accuracy of the confidence intervals. The confidence level for the difference
#' is 1 – 2*alpha, which is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m11		 estimated mean for group 1 in original study 
#' @param    m12		 estimated mean for group 2 in original study
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
#' * Row 3 estimates the difference in mean differences
#' * Row 4 estimates the average mean difference
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
#' #                       Estimate        SE      t     p        LL       UL     df
#' # Original:                 5.80 0.7889312  7.352 0.000  4.228624 7.371376  75.75
#' # Follow-up:                6.10 0.6346075  9.612 0.000  4.845913 7.354087 147.65
#' # Original - Follow-up:    -0.30 1.0124916 -0.296 0.767 -1.974571 1.374571 169.16
#' # Average:                  5.95 0.5062458 11.753 0.000  4.950627 6.949373 169.16
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
  out1 <- t(c(est1, se1, round(t1, 3), round(pval1, 3), ll1, ul1, round(df1, 2)))
  out2 <- t(c(est2, se2, round(t2, 3), round(pval2, 3), ll2, ul2, round(df2, 2)))
  out3 <- t(c(est3, se3, round(t3, 3), round(pval3, 3), ll3, ul3, round(df3, 2)))
  out4 <- t(c(est4, se4, round(t4, 3), round(pval4, 3), ll4, ul4, round(df3, 2)))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.mean.ps ============================================================
#' Compares and combines paired-samples mean differences in original and 
#' follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals from an original study and a 
#' follow-up study where the effect size is a paired-samples mean difference. 
#' Confidence intervals for the difference and average effect size are also
#' computed. Equality of variances within or across studies is not assumed. 
#' A Satterthwaite adjustment to the degrees of freedom is used to 
#' improve the accuracy of the confidence intervals for the difference and 
#' average. The confidence level for the difference is 1 – 2*alpha, which is
#' recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m11		 estimated mean for measurement 1 in original study 
#' @param    m12		 estimated mean for measurement 2 in original study
#' @param    sd11   	 estimated SD for measurement 1 in original study
#' @param    sd12   	 estimated SD for measurement 2 in original study
#' @param    n1    	     sample size in original study
#' @param    cor1    	 estimated correlation of paired measurements in orginal study
#' @param    m21    	 estimated mean for measurement 1 in follow-up study 
#' @param    m22    	 estimated mean for measurement 2 in follow-up study
#' @param    sd21   	 estimated SD for measurement 1 in follow-up study
#' @param    sd22   	 estimated SD for measurement 2 in follow-up study
#' @param    n2    	     sample size in follow-up study
#' @param    cor2    	 estimated correlation of paired measurements in follow-up study
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in mean differences
#' * Row 4 estimates the average mean difference
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
#' #                        Estimate       SE     t     p        LL       UL   df
#' #  Original:                15.29 2.154344 7.097 0.000 10.780906 19.79909 19.0
#' #  Follow-up:                7.57 1.460664 5.183 0.000  4.659564 10.48044 74.0
#' #  Original - Follow-up:     7.72 2.602832 2.966 0.005  3.332885 12.10712 38.4
#' #  Average:                 11.43 1.301416 8.783 0.000  8.796322 14.06368 38.4
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
  out1 <- t(c(est1, se1, round(t1, 3), round(pval1, 3), ll1, ul1, round(df1, 2)))
  out2 <- t(c(est2, se2, round(t2, 3), round(pval2, 3), ll2, ul2, round(df2, 2)))
  out3 <- t(c(est3, se3, round(t3, 3), round(pval3, 3), ll3, ul3, round(df3, 2)))
  out4 <- t(c(est4, se4, round(t4, 3), round(pval4, 3), ll4, ul4, round(df3, 2)))
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
#' This function computes confidence intervals from an original study and a 
#' follow-up study where the effect size is a 2-group standardized mean 
#' difference. Confidence intervals for the difference and average effect 
#' size are also computed. Equality of variances within or across studies
#' is not assumed. The confidence level for the difference is 1 – 2*alpha, 
#' which is recommended for equivalence testing. Square root unweighted 
#' variances, square root weighted variances, and single-group standard 
#' deviation are options for the standardizer.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	     alpha level for 1-alpha confidence
#' @param    m11	     estimated mean for group 1 in original study 
#' @param    m12	     estimated mean for group 2 in original study
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
#' @param    stdzr
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#' * set to 3 for square root weighted average variance standardizer 
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in standardized mean differences
#' * Row 4 estimates the average standardized mean difference
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
#'                         25.2, 19.1, 3.98, 3.79, 75, 75, 0)
#'
#' # Should return: 
#' #                        Estimate     SE      LL     UL
#' #  Original:               1.6280 0.2595  1.1353 2.1524
#' #  Follow-up:              1.5617 0.1871  1.2030 1.9363
#' #  Original - Follow-up:   0.0742 0.3199 -0.4519 0.6004
#' #  Average:                1.5949 0.1599  1.2814 1.9083
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
replicate.stdmean2 <- function(alpha, m11, m12, sd11, sd12, n11, n12, m21, m22, sd21, sd22, n21, n22, stdzr) {
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
  if (stdzr == 0) {
    s1 <- sqrt((v11 + v12)/2)
    s2 <- sqrt((v21 + v22)/2)
    a1 <- 1 - 3/(4*(n11 + n12) - 9)
    a2 <- 1 - 3/(4*(n21 + n22) - 9)
    est1 <- (m11 - m12)/s1
    est2 <- (m21 - m22)/s2
    se1 <- sqrt(est1^2*(v11^2/df11 + v12^2/df12)/(8*s1^4) + (v11/df11 + v12/df12)/s1^2)
    se2 <- sqrt(est2^2*(v21^2/df21 + v22^2/df22)/(8*s2^4) + (v21/df21 + v22/df22)/s2^2)
  } else if (stdzr == 1) { 
    a1 <- (1 - 3/(4*n11 - 5))
	a2 <- (1 - 3/(4*n21 - 5))
    est1 <- (m11 - m12)/sd11
    est2 <- (m21 - m22)/sd21
	se1 <- sqrt(est1^2/(2*df11) + 1/df11 + v12/(df12*v11))
	se2 <- sqrt(est2^2/(2*df21) + 1/df21 + v22/(df22*v21))
  } else if (stdzr == 2) {
    a1 <- (1 - 3/(4*n12 - 5))
	a2 <- (1 - 3/(4*n22 - 5))
    est1 <- (m11 - m12)/sd12
    est2 <- (m21 - m22)/sd22
	se1 <- sqrt(est1^2/(2*df12) + 1/df12 + v11/(df11*v11))
	se2 <- sqrt(est2^2/(2*df22) + 1/df22 + v21/(df21*v22))
  } else {
    s1 <- sqrt((df11*v11 + df12*v12)/(df11 + df12))
    s2 <- sqrt((df21*v21 + df22*v22)/(df21 + df22))
    a1 <- 1 - 3/(4*(n11 + n12) - 9)
    a2 <- 1 - 3/(4*(n21 + n22) - 9)
    est1 <- (m11 - m12)/s1
    est2 <- (m21 - m22)/s2
    se1 <- sqrt(est1^2*(1/df11 + 1/df12)/8 + (v11/n11 + v12/n12)/s1^2)
    se2 <- sqrt(est2^2*(1/df12 + 1/df22)/8 + (v12/n12 + v22/n22)/s2^2)
  }
  est3 <- est1 - est2
  est4 <- (a1*est1 + a2*est2)/2
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  ll1 <- est1 - zcrit1*se1
  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2
  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3
  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4
  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(a1*est1, se1, ll1, ul1))
  out2 <- t(c(a2*est2, se2, ll2, ul2))
  out3 <- t(c(est3, se3, ll3, ul3))
  out4 <- t(c(est4, se4, ll4, ul4))
  out <- rbind(round(out1, 4), round(out2, 4), round(out3, 4), round(out4, 4))
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
#' This function computes confidence intervals from an original study and a follow-up
#' study where the effect size is a paired-samples standardized mean difference. 
#' Confidence intervals for the difference and average effect size are also computed.
#' Equality of variances within or across studies is not assumed. The confidence level
#' for the difference is 1 – 2*alpha, which is recommended for equivalence testing. 
#' Square root unweighted variances and single-condition standard deviation are options
#' for the standardizer.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	 alpha level for 1-alpha confidence
#' @param    m11	 estimated mean for group 1 in original study 
#' @param    m12	 estimated mean for group 2 in original study
#' @param    sd11   	 estimated SD for group 1 in original study
#' @param    sd12   	 estimated SD for group 2 in original study
#' @param    cor1    	 estimated correlation of paired observations in orginal study
#' @param    n1          sample size in original study
#' @param    m21    	 estimated mean for group 1 in follow-up study 
#' @param    m22    	 estimated mean for group 2 in follow-up study
#' @param    sd21   	 estimated SD for group 1 in follow-up study
#' @param    sd22   	 estimated SD for group 2 in follow-up study
#' @param    cor2    	 estimated correlation of paired observations in follow-up study
#' @param    n2          sample size in follow-up study
#' @param    stdzr
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for measurement 1 SD standardizer 
#' * set to 2 for measurement 2 SD standardizer 
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in standardized mean differences
#' * Row 4 estimates the average standardized mean difference
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
#' replicate.stdmean.ps(alpha = .05, 86.22, 70.93, 14.89, 12.32, .765, 20, 
#'                                   84.81, 77.24, 15.68, 16.95, .702, 75, 0)
#'
#' # Should return:
#' #                       Estimate     SE     LL     UL
#' # Orginal:                1.0890 0.2292 0.6697 1.5680
#' # Follow-up:              0.4605 0.0959 0.2757 0.6516
#' # Original - Follow-up:   0.6552 0.2484 0.2466 1.0638
#' # Average:                0.7748 0.1242 0.5313 1.0182
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
replicate.stdmean.ps <- function(alpha, m11, m12, sd11, sd12, cor1, n1, m21, m22, sd21, sd22, cor2, n2, stdzr) {
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  v11 <- sd11^2
  v12 <- sd12^2
  v21 <- sd21^2
  v22 <- sd22^2
  df1 <- n1 - 1
  df2 <- n2 - 1
  vd1 <- v11 + v12 - 2*cor1*sd11*sd12
  vd2 <- v21 + v22 - 2*cor2*sd21*sd22
  if (stdzr == 0) {
    s1 <- sqrt((v11 + v12)/2)
    s2 <- sqrt((v21 + v22)/2)
    a1 <- sqrt((n1 - 2)/df1)
    a2 <- sqrt((n2 - 2)/df2)
    est1 <- (m11 - m12)/s1
    est2 <- (m21 - m22)/s2
    est3 <- est1 - est2
    est4 <- (a1*est1 + a2*est2)/2
    se1 <- sqrt(est1^2*(v11^2 + v12^2 + 2*cor1^2*v11*v12)/(8*df1*s1^4) + vd1/(df1*s1^2))
    se2 <- sqrt(est2^2*(v21^2 + v22^2 + 2*cor2^2*v21*v22)/(8*df2*s2^4) + vd2/(df2*s2^2))
  } else if (stdzr == 1){
    a1 <- 1 - 3/(4*df1 - 1)
    a2 <- 1 - 3/(4*df2 - 1)
    est1 <- (m11 - m12)/sd11
    est2 <- (m21 - m22)/sd21
    est3 <- est1 - est2
    est4 <- (a1*est1 + a2*est2)/2
	se1 <- sqrt(est1^2/(2*df1) + vd1/(df1*v11))
	se2 <- sqrt(est2^2/(2*df2) + vd2/(df2*v12))
  } else {
    a1 <- 1 - 3/(4*df1 - 1)
	a2 <- 1 - 3/(4*df2 - 1)
    est1 <- (m11 - m12)/sd12
    est2 <- (m21 - m22)/sd22
    est3 <- est1 - est2
    est4 <- (a1*est1 + a2*est2)/2
    se1 <- sqrt(est1^2/(2*df1) + vd1/(df1*v12))
	se2 <- sqrt(est2^2/(2*df2) + vd2/(df2*v22))
  }
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  ll1 <- est1 - zcrit1*se1
  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2
  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3
  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4
  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(a1*est1, se1, ll1, ul1))
  out2 <- t(c(a2*est2, se2, ll2, ul2))
  out3 <- t(c(est3, se3, ll3, ul3))
  out4 <- t(c(est4, se4, ll4, ul4))
  out <- rbind(round(out1, 4), round(out2, 4), round(out3, 4), round(out4, 4))
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Orginal:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.cor ============================================================
#' Compares and combines Pearson or partial correlations in original and 
#' follow-up studies
#' 
#'
#' @description 
#' This function can be used to compare and combine Pearson or partial 
#' correlations from an original study and a follow-up study. The 
#' confidence level for the difference is 1 – 2*alpha, which is recommended 
#' for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	 alpha level for 1-alpha confidence
#' @param    cor1  	 estimated correlation in original study
#' @param    n1    	 sample size in original study
#' @param    cor2  	 estimated correlation in follow-up study
#' @param    n2    	 sample size in follow-up study
#' @param    s     	 number of control variables in each study (0 for Pearson)
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in correlations
#' * Row 4 estimates the average correlation
#' 
#' 
#' The columns are:
#' * Estimate -correlation estimate (single study, difference, average)
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
#' #                       Estimate      SE     z     p     LL     UL
#' # Original:                0.598 0.07321 6.589 0.000 0.4355 0.7228
#' # Follow-up:               0.324 0.06377 4.819 0.000 0.1940 0.4428
#' # Original - Follow-up:    0.274 0.09709 2.633 0.008 0.1065 0.4265
#' # Average:                 0.461 0.04854 7.635 0.000 0.3725 0.5412
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
  out1 <- t(c(round(cor1, 4), round(se1, 5), round(t1, 3), round(pval1, 3), round(ll1a, 4), round(ul1a, 4)))
  out2 <- t(c(round(cor2, 4), round(se2, 5), round(t2, 3), round(pval2, 3), round(ll2a, 4), round(ul2a, 4)))
  out3 <- t(c(round(dif, 4), round(se3, 5), round(t3, 3), round(pval3, 3), round(ll3, 4), round(ul3, 4)))
  out4 <- t(c(round(ave, 4), round(se4, 5), round(t4, 3), round(pval4, 3), round(ll4, 4), round(ul4, 4)))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


# replicate.prop2 ============================================================
#' Compares and combines 2-group proportion differences in original and 
#' follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals from an original study and a
#' follow-up study where the effect size is a 2-group proportion difference. 
#' Confidence intervals for the difference and average effect size are also 
#' computed. The confidence level for the difference is 1 – 2*alpha, which 
#' is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    f11		 frequency count for group 1 in original study 
#' @param    f12		 frequency count for group 2 in original study
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
#' * Row 3 estimates the difference in proportion differences
#' * Row 4 estimates the average proportion difference
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
#' #                         Estimate         SE     z     p          LL        UL
#' # Original:             0.11904762 0.10805233 1.102 0.271 -0.09273105 0.3308263
#' # Follow-up:            0.09677419 0.07965047 1.215 0.224 -0.05933787 0.2528863
#' # Original - Follow-up: 0.02359056 0.13542107 0.174 0.862 -0.19915727 0.2463384
#' # Average:              0.11015594 0.06771053 1.627 0.104 -0.02255427 0.2428661
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
  pval1 <- 2*(1 - pnorm(abs(z1)))
  pval2 <- 2*(1 - pnorm(abs(z2)))
  pval3 <- 2*(1 - pnorm(abs(z3))) 
  pval4 <- 2*(1 - pnorm(abs(z4)))
  ll1 <- est1 - zcrit1*se1;  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2;  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3;  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4;  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(est1, se1, round(z1, 3), round(pval1, 3), ll1, ul1))
  out2 <- t(c(est2, se2, round(z2, 3), round(pval2, 3), ll2, ul2))
  out3 <- t(c(est3, se3, round(z3, 3), round(pval3, 3), ll3, ul3))
  out4 <- t(c(est4, se4, round(z4, 3), round(pval4, 3), ll4, ul4))
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
#' ratio of odds ratios and geometric average odds ratio are also  
#' computed. The confidence level for the ratio of ratios is 1 – 2*alpha, which
#' is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
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
#' * Row 3 estimates the ratio of odds ratios 
#' * Row 4 estimates the geometric average odds ratio
#'
#'
#' The columns are:
#' * Estimate - odds ratio estimate (single study, ratio, average)
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
#' #                        Estimate        SE      z     p   exp(LL)  exp(UL)
#' # Original:            1.39000000 0.3020000  4.603 0.000 2.2212961 7.256583
#' # Follow-up:           1.48000000 0.2060000  7.184 0.000 2.9336501 6.578144
#' # Original/Follow-up: -0.06273834 0.3655681 -0.172 0.864 0.5147653 1.713551
#' # Average:             0.36067292 0.1827840  1.973 0.048 1.0024257 2.052222
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
  pval1 <- 2*(1 - pnorm(abs(z1)))
  pval2 <- 2*(1 - pnorm(abs(z2)))
  pval3 <- 2*(1 - pnorm(abs(z3))) 
  pval4 <- 2*(1 - pnorm(abs(z4)))
  ll1 <- exp(est1 - zcrit1*se1);  ul1 <- exp(est1 + zcrit1*se1)
  ll2 <- exp(est2 - zcrit1*se2);  ul2 <- exp(est2 + zcrit1*se2)
  ll3 <- exp(est3 - zcrit2*se3);  ul3 <- exp(est3 + zcrit2*se3)
  ll4 <- exp(est4 - zcrit1*se4);  ul4 <- exp(est4 + zcrit1*se4)
  out1 <- t(c(est1, se1, round(z1, 3), round(pval1, 3), ll1, ul1))
  out2 <- t(c(est2, se2, round(z2, 3), round(pval2, 3), ll2, ul2))
  out3 <- t(c(est3, se3, round(z3, 3), round(pval3, 3), ll3, ul3))
  out4 <- t(c(est4, se4, round(z4, 3), round(pval4, 3), ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "exp(LL)", "exp(UL)")
  rownames(out) <- c("Original:", "Follow-up:", "Original/Follow-up:", "Average:")
  return(out)
}


# replicate.slope ============================================================ 
#' Compares and combines slope coefficients in original and follow-up studies
#'
#' @description 
#' This function computes confidence intervals for an OLS slope in a GLM from 
#' the original and follow-up studies, the difference in slopes, and the average
#' of the slopes. Equality of error variances across studies is not assumed. The 
#' confidence interval for the difference uses a 1 - 2*alpha confidence level, 
#' which is recommended for equivalence testing. Use the \link[vcmeta]{replicate.gen} 
#' function for slopes in other types of models (e.g., binary logistic, ordinal 
#' logistic, SEM). A Satterthwaite adjustment to the degrees of freedom is used
#' to improve the accuracy of the confidence intervals for the average and the
#' difference.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
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
#' * Row 3 estimates the difference in slopes
#' * Row 4 estimates the average slope
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
#' #                       Estimate       SE     t     p        LL       UL    df
#' # Original:                23.40 5.160000 4.535 0.000 13.007227 33.79277  45.0
#' # Follow-up:               18.50 4.480000 4.129 0.000  9.592560 27.40744  85.0
#' # Original - Follow-up:     4.90 6.833447 0.717 0.475 -6.438743 16.23874 106.4
#' # Average:                 20.95 3.416724 6.132 0.000 14.176310 27.72369 106.4
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
 out1 <- t(c(est1, se1, round(t1, 3), round(pval1, 3), ll1, ul1, round(df1, 2)))
 out2 <- t(c(est2, se2, round(t2, 3), round(pval2, 3), ll2, ul2, round(df2, 2)))
 out3 <- t(c(est3, se3, round(t3, 3), round(pval3, 3), ll3, ul3, round(df3, 2)))
 out4 <- t(c(est4, se4, round(t4, 3), round(pval4, 3), ll4, ul4, round(df3, 2)))
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
#' This function can be used to compare and combine any effect size using the 
#' effect size estimate and its standard error from the original study and 
#' the follow-up study. The confidence level for the difference is 1 – 2*alpha,
#' which is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#'
#' @param    alpha		 alpha level for 1-alpha confidence 
#' @param    est1  	     estimated effect size in original study
#' @param    se1    	 effect size standard error in original study
#' @param    est2  	     estimated effect size in follow-up study
#' @param    se2    	 effect size standard error in follow-up study
#'   
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in effect sizes
#' * Row 4 estimates the average effect size
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
#' #                       Estimate        SE     z     p         LL        UL
#' # Original:                0.782 0.2100000 3.724 0.000  0.3704076 1.1935924
#' # Follow-up:               0.650 0.1540000 4.221 0.000  0.3481655 0.9518345
#' # Original - Follow-up:    0.132 0.2604151 0.507 0.612 -0.2963446 0.5603446
#' # Average:                 0.716 0.1302075 5.499 0.000  0.4607979 0.9712021
#' 
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
  out1 <- t(c(est1, se1, round(z1, 3), round(pval1, 3), ll1, ul1))
  out2 <- t(c(est2, se2, round(z2, 3), round(pval2, 3), ll2, ul2))
  out3 <- t(c(est3, se3, round(z3, 3), round(pval3, 3), ll3, ul3))
  out4 <- t(c(est4, se4, round(z4, 3), round(pval4, 3), ll4, ul4))
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
#' an original study and a follow-up study. The confidence level for the 
#' difference is 1 – 2*alpha, which is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
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
#' * Row 3 estimates the difference in correlations
#' * Row 4 estimates the average correlation 
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
#' #                       Estimate      SE     z     p     LL     UL
#' # Original:                0.598 0.07948 5.315 0.000 0.4199 0.7318
#' # Follow-up:               0.324 0.06542 4.571 0.000 0.1905 0.4457
#' # Original - Follow-up:    0.274 0.10294 3.438 0.001 0.0948 0.4342
#' # Average:                 0.461 0.05147 9.968 0.000 0.3670 0.5457
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
  z1 <- cor1*sqrt(n1 - 1) 
  z2 <- cor2*sqrt(n2 - 1)
  z3 <- (zr1 - zr2)/sqrt(se1^2 + se2^2)
  z4 <- (zr1 + zr2)/sqrt(se1^2 + se2^2)
  pval1 <- 2*(1 - pnorm(abs(z1)))
  pval2 <- 2*(1 - pnorm(abs(z2)))
  pval3 <- 2*(1 - pnorm(abs(z3)))
  pval4 <- 2*(1 - pnorm(abs(z4)))
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
  out1 <- t(c(round(cor1, 4), round(se1, 5), round(z1, 3), round(pval1, 3), round(ll1a, 4), round(ul1a, 4)))
  out2 <- t(c(round(cor2, 4), round(se2, 5), round(z2, 3), round(pval2, 3), round(ll2a, 4), round(ul2a, 4)))
  out3 <- t(c(round(dif, 4), round(se3, 5), round(z3, 3), round(pval3, 3), round(ll3, 4), round(ul3, 4)))
  out4 <- t(c(round(ave, 4), round(se4, 5), round(z4, 3), round(pval4, 3), round(ll4, 4), round(ul4, 4)))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.prop1 ============================================================
#' Compares and combines single proportions in original and follow-up studies
#' 
#'
#' @description 
#' This function computes confidence intervals for a single proportion from an 
#' original study and a follow-up study. Confidence intervals for the
#' difference between the two proportions and average of the two proportions 
#' are also computed. The confidence level for the difference is 1 – 2*alpha, 
#' which is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	 alpha level for 1-alpha confidence
#' @param    f1	  	 frequency count in original study 
#' @param    n1      sample size in original study
#' @param    f2 	 frequency count in follow-up study 
#' @param    n2    	 sample size for in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in proportions
#' * Row 4 estimates the average proportion
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
#' difference between the two means and average of the two means are also
#' computed. Equality of variances across studies is not assumed. A 
#' Satterthwaite adjustment to the degrees of freedom is used to improve 
#' the accuracy of the confidence intervals for the difference and average. 
#' The confidence level for the difference is 1 – 2*alpha, which is 
#' recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha		 alpha level for 1-alpha confidence
#' @param    m1	  	     estimated mean in original study 
#' @param    sd1   	     estimated SD in original study
#' @param    n1     	 sample size in original study
#' @param    m2 	   	 estimated mean in follow-up study 
#' @param    sd2   		 estimated SD in follow-up study
#' @param    n2    	 	 sample size for in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in means
#' * Row 4 estimates the average mean
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
#' #                       Estimate        SE        LL        UL    df
#' # Original:                21.90 0.6039950 20.678305 23.121695 39.00
#' # Follow-up:               25.20 0.4595708 24.284285 26.115715 74.00
#' # Original - Follow-up:    -3.30 0.7589567 -4.562527 -2.037473 82.63
#' # Average:                 23.55 0.3794784 22.795183 24.304817 82.63
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
  out1 <- t(c(est1, se1, ll1, ul1, round(df1, 2)))
  out2 <- t(c(est2, se2, ll2, ul2, round(df2, 2)))
  out3 <- t(c(est3, se3, ll3, ul3, round(df3, 2)))
  out4 <- t(c(est4, se4, ll4, ul4, round(df3, 3)))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


# replicate.propratio2 =======================================================
#' Compares and combines 2-group proportion ratios in original and follow-up 
#' studies
#' 
#'
#' @description 
#' This function computes confidence intervals from an original study and a 
#' follow-up study where the effect size is a 2-group proportion ratio. 
#' Confidence intervals for the ratio and geometric average of effect sizes
#' are also computed. The confidence level for the ratio of ratios is 1 – 2*alpha, 
#' which is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	     alpha level for 1-alpha confidence																																												
#' @param    f11		 frequency count for group 1 in original study 
#' @param    f12		 frequency count for group 2 in original study
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
#' * Row 3 estimates the ratio of proportion ratios
#' * Row 4 estimates the geometric average proportion ratio
#'
#'
#' The columns are:
#' * Estimate - proportion ratio estimate (single study, ratio, average)
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'    
#' 
#' @examples
#' replicate.propratio2(.05, 21, 16, 40, 40, 19, 13, 60, 60)
#'
#' # Should return:
#' #                      Estimate        LL       UL
#' # Original:           1.3076923 0.8068705 2.119373
#' # Follow-up:          1.4528302 0.7939881 2.658372
#' # Original/Follow-up: 0.9000999 0.4703209 1.722611
#' # Average:            1.3783522 0.9362893 2.029132
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
replicate.propratio2 <- function(alpha, f11, f12, n11, n12, f21, f22, n21, n22){
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  p11 <- (f11 + 1/4)/(n11 + 7/4)
  p12 <- (f12 + 1/4)/(n12 + 7/4)
  v11 <- 1/(f11 + 1/4 + (f11 + 1/4)^2/(n11 - f11 + 3/2))
  v12 <- 1/(f12 + 1/4 + (f12 + 1/4)^2/(n12 - f12 + 3/2))
  se1 <- sqrt(v11 + v12)
  est1 <- log(p11/p12)
  p21 <- (f21 + 1/4)/(n21 + 7/4)
  p22 <- (f22 + 1/4)/(n22 + 7/4)
  v21 <- 1/(f21 + 1/4 + (f21 + 1/4)^2/(n21 - f21 + 3/2))
  v22 <- 1/(f22 + 1/4 + (f22 + 1/4)^2/(n22 - f22 + 3/2))
  se2 <- sqrt(v21 + v22)
  est2 <- log(p21/p22)
  est3 <- est1 - est2
  est4 <- (est1 + est2)/2
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- se3/2
  ll1 <- exp(est1 - zcrit1*se1);  ul1 <- exp(est1 + zcrit1*se1)
  ll2 <- exp(est2 - zcrit1*se2);  ul2 <- exp(est2 + zcrit1*se2)
  ll3 <- exp(est3 - zcrit2*se3);  ul3 <- exp(est3 + zcrit2*se3)
  ll4 <- exp(est4 - zcrit1*se4);  ul4 <- exp(est4 + zcrit1*se4)
  out1 <- t(c(exp(est1), ll1, ul1))
  out2 <- t(c(exp(est2), ll2, ul2))
  out3 <- t(c(exp(est3), ll3, ul3))
  out4 <- t(c(exp(est4), ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original/Follow-up:", "Average:")
  return(out)
}


# replicate.prop.ps ===========================================================
#' Compares and combines paired-samples proportion differences in original and 
#' follow-up studies
#'                   
#'
#' @description 
#' This function computes confidence intervals from an original study and a 
#' follow-up study where the effect size is a paired-samples proportion 
#' difference. Confidence intervals for the difference and average of effect
#' sizes are also computed. The confidence level for the difference is
#' 1 – 2*alpha, which is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	   alpha level for 1-alpha confidence
#' @param    f1		   vector of frequency counts for 2x2 table in original study 
#' @param    f2		   vector of frequency counts for 2x2 table in follow-up study
#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in proportion differences
#' * Row 4 estimates the average proportion difference 
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
#' f1 <- c(42, 2, 15, 61)
#' f2 <- c(69, 5, 31, 145)
#' replicate.prop.ps(.05, f1, f2)
#'
#' # Should return:
#' #                          Estimate         SE     z     p          LL         UL
#' # Original:             0.106557377 0.03440159 3.097 0.002  0.03913151 0.17398325
#' # Follow-up:            0.103174603 0.02358274 4.375 0.000  0.05695329 0.14939592
#' # Original - Follow-up: 0.003852359 0.04097037 0.094 0.925 -0.06353791 0.07124263
#' # Average:              0.105511837 0.02048519 5.151 0.000  0.06536161 0.14566206
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#'
#' \insertRef{Bonett2012}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
replicate.prop.ps <- function(alpha, f1, f2){
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  n1 <- sum(f1)
  p01 <- (f1[2] + 1)/(n1 + 2)
  p10 <- (f1[3] + 1)/(n1 + 2)
  est1 <- p10 - p01
  se1 <- sqrt(((p01 + p10) - (p01 - p10)^2)/(n1 + 2))
  n2 <- sum(f2)
  p01 <- (f2[2] + 1)/(n2 + 2)
  p10 <- (f2[3] + 1)/(n2 + 2)
  est2 <- p10 - p01
  se2 <- sqrt(((p01 + p10) - (p01 - p10)^2)/(n2 + 2))
  p011 <- (f1[2] + .5)/(n1 + 1)
  p101 <- (f1[3] + .5)/(n1 + 1)
  p012 <- (f2[2] + .5)/(n2 + 1)
  p102 <- (f2[3] + .5)/(n2 + 1)
  est3 <- p101 - p011 - p102 + p012
  v1 = ((p101 + p011) - (p101 - p011)^2)/(n1 + 1)
  v2 = ((p102 + p012) - (p102 - p012)^2)/(n2 + 1)
  se3 <- sqrt(v1 + v2)
  est4 <- ((p101 - p011) + (p102 - p012))/2
  se4 <- se3/2
  z1 <- est1/se1
  z2 <- est2/se2
  z3 <- est3/se3
  z4 <- est4/se4
  pval1 <- 2*(1 - pnorm(abs(z1)))
  pval2 <- 2*(1 - pnorm(abs(z2)))
  pval3 <- 2*(1 - pnorm(abs(z3))) 
  pval4 <- 2*(1 - pnorm(abs(z4)))
  ll1 <- est1 - zcrit1*se1;  ul1 <- est1 + zcrit1*se1
  ll2 <- est2 - zcrit1*se2;  ul2 <- est2 + zcrit1*se2
  ll3 <- est3 - zcrit2*se3;  ul3 <- est3 + zcrit2*se3
  ll4 <- est4 - zcrit1*se4;  ul4 <- est4 + zcrit1*se4
  out1 <- t(c(est1, se1, round(z1, 3), round(pval1, 3), ll1, ul1))
  out2 <- t(c(est2, se2, round(z2, 3), round(pval2, 3), ll2, ul2))
  out3 <- t(c(est3, se3, round(z3, 3), round(pval3, 3), ll3, ul3))
  out4 <- t(c(est4, se4, round(z4, 3), round(pval4, 3), ll4, ul4))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.cor.gen ==========================================================
#' Compares and combines any type of correlation in original and 
#' follow-up studies
#' 
#'
#' @description 
#' This function can be used to compare and combine any type of correlation 
#' from an original study and a follow-up study. The confidence level for the 
#' difference is 1 – 2*alpha, which is recommended for equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	 alpha level for 1-alpha confidence
#' @param    cor1  	 estimated correlation in original study
#' @param    se1   	 standard error of correlation in original study
#' @param    cor2  	 estimated correlation in follow-up study
#' @param    se2   	 standard error of correlation in follow-up study
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in correlations
#' * Row 4 estimates the average correlation
#' 
#' 
#' The columns are:
#' * Estimate - correlation estimate (single study, difference, average)
#' * SE - standard error
#' * z - z-value 
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.cor.gen(.05, .454, .170, .318, .098)
#'
#' # Should return:
#' #                       Estimate      SE     z     p      LL     UL
#' # Original:                0.454 0.17000 2.287 0.022  0.0699 0.7209
#' # Follow-up:               0.318 0.09800 3.022 0.003  0.1152 0.4953
#' # Original - Follow-up:    0.136 0.19622 0.667 0.505 -0.2154 0.4237
#' # Average:                 0.386 0.09811 3.409 0.001  0.1961 0.5480
#' 
#' 
#' @references
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @export
replicate.cor.gen <- function(alpha, cor1, se1, cor2, se2) {
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  zr1 <- log((1 + cor1)/(1 - cor1))/2
  zr2 <- log((1 + cor2)/(1 - cor2))/2
  dif <- cor1 - cor2
  ave <- (cor1 + cor2)/2
  ave.z <- log((1 + ave)/(1 - ave))/2
  se1.z <- se1/(1 - cor1^2)
  se2.z <- se2/(1 - cor2^2)
  se3 <- sqrt(se1^2 + se2^2)
  se4 <- sqrt(se1^2 + se2^2)/2
  se4.z <- sqrt(((se1^2 + se2^2)/4)/(1 - ave^2))
  z1 <- zr1/se1.z 
  z2 <- zr2/se2.z 
  z3 <- (zr1 - zr2)/sqrt(se1.z^2 + se2.z^2)
  z4 <- (zr1 + zr2)/sqrt(se1.z^2 + se2.z^2)
  pval1 <- 2*(1 - pnorm(abs(z1)))
  pval2 <- 2*(1 - pnorm(abs(z2)))
  pval3 <- 2*(1 - pnorm(abs(z3)))
  pval4 <- 2*(1 - pnorm(abs(z4)))
  ll0a <- zr1 - zcrit1*se1.z
  ul0a <- zr1 + zcrit1*se1.z
  ll1a <- (exp(2*ll0a) - 1)/(exp(2*ll0a) + 1)
  ul1a <- (exp(2*ul0a) - 1)/(exp(2*ul0a) + 1)
  ll0b <- zr1 - zcrit2*se1.z
  ul0b <- zr1 + zcrit2*se1.z
  ll1b <- (exp(2*ll0b) - 1)/(exp(2*ll0b) + 1)
  ul1b <- (exp(2*ul0b) - 1)/(exp(2*ul0b) + 1)
  ll0a <- zr2 - zcrit1*se2.z
  ul0a <- zr2 + zcrit1*se2.z
  ll2a <- (exp(2*ll0a) - 1)/(exp(2*ll0a) + 1)
  ul2a <- (exp(2*ul0a) - 1)/(exp(2*ul0a) + 1)
  ll0b <- zr2 - zcrit2*se2.z
  ul0b <- zr2 + zcrit2*se2.z
  ll2b <- (exp(2*ll0b) - 1)/(exp(2*ll0b) + 1)
  ul2b <- (exp(2*ul0b) - 1)/(exp(2*ul0b) + 1)
  ll3 <- dif - sqrt((cor1 - ll1b)^2 + (ul2b - cor2)^2)
  ul3 <- dif + sqrt((ul1b - cor1)^2 + (cor2 - ll2b)^2)
  ll0 <- ave.z - zcrit1*se4.z
  ul0 <- ave.z + zcrit1*se4.z
  ll4 <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul4 <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out1 <- t(c(round(cor1, 4), round(se1, 5), round(z1, 3), round(pval1, 3), round(ll1a, 4), round(ul1a, 4)))
  out2 <- t(c(round(cor2, 4), round(se2, 5), round(z2, 3), round(pval2, 3), round(ll2a, 4), round(ul2a, 4)))
  out3 <- t(c(round(dif, 4), round(se3, 5), round(z3, 3), round(pval3, 3), round(ll3, 4), round(ul3, 4)))
  out4 <- t(c(round(ave, 4), round(se4, 5), round(z4, 3), round(pval4, 3), round(ll4, 4), round(ul4, 4)))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


# replicate.agree =============================================================
#' Compares and combines G-index of agreement in original and follow-up studies
#' 
#'
#' @description 
#' This function computes adjusted Wald confidence intervals from an original 
#' study and a follow-up study where the effect size is a G-index of agreement. 
#' Adjusted Wald confidence intervals for the difference and average effect 
#' size are also computed. The confidence level for the difference is 
#' 1 – 2*alpha, which is recommended for equivalence testing. As a measurement
#' of agreement, the G-index is usually preferred to Cohen's kappa.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	     alpha level for 1-alpha confidence
#' @param    f1		     number of objects rated in agreement in original study 
#' @param    n1		     sample size (number of objects) in original study
#' @param    f2		     number of objects rated in agreement in follow-up study 
#' @param    n2		     sample size (number of objects) in follow-up study
#' @param    k		     number of rating categories

#' 
#' 
#' @return A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in G-indicies
#' * Row 4 estimates the average G-index
#'
#'
#' The columns are:
#' * Estimate - MLE of G-index (single study, difference, average)
#' * SE - standard error of adjusted estimate
#' * LL - lower limit of the adjusted confidence interval
#' * UL - upper limit of the adjusted confidence interval
#'    
#' 
#' @examples
#' replicate.agree(.05, 85, 100, 160, 200, 2)
#'
#' # Should return:
#' #                       Estimate      SE      LL     UL
#' # Original:                 0.70 0.07252  0.5309 0.8152
#' # Follow-up:                0.60 0.05662  0.4773 0.6992
#' # Original - Follow-up:     0.10 0.09160 -0.0584 0.2429
#' # Average:                  0.65 0.04580  0.5504 0.7299
#' 
#' 
#' @references
#' \insertRef{Bonett2022}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
replicate.agree <- function(alpha, f1, n1, f2, n2, k){
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  a <- k/(k - 1)
  p1 <- f1/n1
  p2 <- f2/n2
  p1.adj1 <- (f1 + 2)/(n1 + 4)
  p2.adj1 <- (f2 + 2)/(n2 + 4)
  p1.adj2 <- (f1 + 1)/(n1 + 2)
  p2.adj2 <- (f2 + 1)/(n2 + 2)
  est1 <- a*p1 - 1/(k - 1)
  est2 <- a*p2 - 1/(k - 1)
  est3 <- est1 - est2
  est4 <- (est1 + est2)/2
  est1.adj <- a*p1.adj1 - 1/(k - 1)
  est2.adj <- a*p2.adj1 - 1/(k - 1)
  est3.adj <- a*p1.adj2 - a*p2.adj2
  est4.adj <- (a*p1.adj2 + a*p2.adj2 - 2/(k - 1))/2
  se1 <- a*sqrt(p1.adj1*(1 - p1.adj1)/(n1 + 4))
  se2 <- a*sqrt(p2.adj1*(1 - p2.adj1)/(n2 + 4))
  se3 <- a*sqrt(p1.adj2*(1 - p1.adj2)/(n1 + 2) + p2.adj2*(1 - p2.adj2)/(n2 + 2))
  se4 <- se3/2
  ll1 <- est1.adj - zcrit1*se1;  ul1 <- est1.adj + zcrit1*se1
  ll2 <- est2.adj - zcrit1*se2;  ul2 <- est2.adj + zcrit1*se2
  ll3 <- est3.adj - zcrit2*se3;  ul3 <- est3.adj + zcrit2*se3
  ll4 <- est4.adj - zcrit1*se4;  ul4 <- est4.adj + zcrit1*se4
  out1 <- t(c(round(est1, 4), round(se1, 5), round(ll1, 4), round(ul1, 4)))
  out2 <- t(c(round(est2, 4), round(se2, 5), round(ll2, 4), round(ul2, 4)))
  out3 <- t(c(round(est3, 4), round(se3, 5), round(ll3, 4), round(ul3, 4)))
  out4 <- t(c(round(est4, 4), round(se4, 5), round(ll4, 4), round(ul4, 4)))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}


#  replicate.cronbach =========================================================
#' Compares and combines Cronbach reliablity in original and follow-up studies
#' 
#'
#' @description 
#' This function can be used to compare and combine a Cronbach reliablity 
#' coefficient from an original study and a follow-up study. The confidence 
#' level for the difference is 1 – 2*alpha, which is recommended for 
#' equivalence testing.
#'
#' For more details, see Chapter 4 of Bonett (2021, Volume 5).
#' 
#' 
#' @param    alpha	 alpha level for 1-alpha confidence
#' @param    rel1  	 estimated reliability in original study
#' @param    n1   	 sample size in original study
#' @param    rel2  	 estimated reliability in follow-up study
#' @param    n2   	 sample size in follow-up study
#' @param    r   	 number of measurements (e.g., items) 
#' 
#' 
#' @return
#' A 4-row matrix. The rows are:
#' * Row 1 summarizes the original study
#' * Row 2 summarizes the follow-up study
#' * Row 3 estimates the difference in correlations
#' * Row 4 estimates the average correlation
#' 
#' 
#' The columns are:
#' * Estimate - correlation estimate (single study, difference, average)
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' replicate.cronbach(.05, .883, 100, .869, 200, 6)
#'
#' # Should return:
#' #                       Estimate      SE      LL     UL
#' # Original:                0.883 0.01831  0.8436 0.9152
#' # Follow-up:               0.869 0.01442  0.8387 0.8952
#' # Original - Follow-up:    0.014 0.02331 -0.0334 0.0582
#' # Average:                 0.876 0.01172  0.8519 0.8977
#' 
#' 
#' @references
#' \insertRef{Bonett2010}{vcmeta}
#'
#' \insertRef{Bonett2015}{vcmeta}
#'
#' \insertRef{Bonett2021}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @importFrom stats qf
#' @export
replicate.cronbach <- function(alpha, rel1, n1, rel2, n2, r) {
  zcrit1 <- qnorm(1 - alpha/2)
  zcrit2 <- qnorm(1 - alpha)
  se1 <- sqrt((2*r*(1 - rel1)^2)/((r - 1)*(n1 - 2)))
  df11 <- n1 - 1
  df12 <- n1*(r - 1)
  f1 <- qf(1 - alpha/2, df11, df12)
  f2 <- qf(1 - alpha/2, df12, df11)
  f0 <- 1/(1 - rel1)
  ll1 <- 1 - f1/f0
  ul1 <- 1 - 1/(f0*f2)
  se2 <- sqrt((2*r*(1 - rel2)^2)/((r - 1)*(n2 - 2)))
  df21 <- n2 - 1
  df22 <- n2*(r - 1)
  f1 <- qf(1 - alpha/2, df21, df22)
  f2 <- qf(1 - alpha/2, df22, df21)
  f0 <- 1/(1 - rel2)
  ll2 <- 1 - f1/f0
  ul2 <- 1 - 1/(f0*f2)
  se3 <- sqrt(se1^2 + se2^2)
  dif <- rel1 - rel2
  ll3 <- dif - sqrt((rel1 - ll1)^2 + (ul2 - rel2)^2)
  ul3 <- dif + sqrt((ul1 - rel1)^2 + (rel2 - ll2)^2)
  hn <- 2/(1/n1 + 1/n2)
  a <- ((r - 2)*(2 - 1))^.25
  se10 <- sqrt((2*r*(1 - rel1)^2)/((r - 1)*(n1 - 2 - a)))
  se20 <- sqrt((2*r*(1 - rel2)^2)/((r - 1)*(n2 - 2 - a)))
  se4 <- sqrt(se10^2 + se20^2)/2
  ave <- (rel1 + rel2)/2
  log.ave <- log(1 - ave) - log(hn/(hn - 1))
  ul4 <- 1 - exp(log.ave - zcrit1*se4/(1 - ave))
  ll4 <- 1 - exp(log.ave + zcrit1*se4/(1 - ave))
  out1 <- t(c(round(rel1, 4), round(se1, 5), round(ll1, 4), round(ul1, 4)))
  out2 <- t(c(round(rel2, 4), round(se2, 5), round(ll2, 4), round(ul2, 4)))
  out3 <- t(c(round(dif, 4), round(se3, 5), round(ll3, 4), round(ul3, 4)))
  out4 <- t(c(round(ave, 4), round(se4, 5), round(ll4, 4), round(ul4, 4)))
  out <- rbind(out1, out2, out3, out4)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Original:", "Follow-up:", "Original - Follow-up:", "Average:")
  return(out)
}

