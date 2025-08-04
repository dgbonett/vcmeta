# ================= Sub-group Comparison of Effect Sizes ============
#  meta.sub.cor =====================================================
#' Confidence interval for a subgroup difference in average Pearson 
#' or partial correlations
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' difference in average Pearson or partial correlations for two mutually
#' exclusive subgroups of studies. Each subgroup can have one or more 
#' studies. All of the correlations must be either Pearson correlations
#' or partial correlations.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     n      	vector of sample sizes 
#' @param     cor    	vector of estimated correlations 
#' @param     s      	number of control variables (set to 0 for Pearson)
#' @param     group  	vector of group indicators:
#' * 1 for set A
#' * 2 for set B
#' * 0 to ignore
#' 
#'
#' @return 
#' Returns a matrix with three rows:
#' * Row 1 - estimate for Set A
#' * Row 2 - estimate for Set B
#' * Row 3 - estimate for difference, Set A - Set B
#'
#' The columns are:
#' * Estimate - estimated average correlation or difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n <- c(55, 190, 65, 35)
#' cor <- c(.40, .65, .60, .45)
#' group <- c(1, 1, 2, 0)
#' meta.sub.cor(.05, n, cor, 0, group)
#' 
#' # Should return:
#' #                Estimate         SE         LL        UL
#' # Set A:            0.525 0.06195298  0.3932082 0.6356531
#' # Set B:            0.600 0.08128008  0.4171458 0.7361686
#' # Set A - Set B:   -0.075 0.10219894 -0.2645019 0.1387283
#' 
#' 
#' @references
#' \insertRef{Bonett2008a}{vcmeta}
#'
#'
#' @importFrom stats qnorm
#' @export
meta.sub.cor <- function(alpha, n, cor, s, group) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  var <- (1 - cor^2)^2/(n - 3 - s)
  g1 <- (group == rep(1, m))*1
  g2 <- (group == rep(2, m))*1
  m1 <- sum(g1)
  m2 <- sum(g2)
  ave.A <- sum(g1*cor)/m1
  se.ave.A <- sqrt(sum(g1*var)/m1^2)
  z.A <- log((1 + ave.A)/(1 - ave.A))/2
  ll0.A <- z.A - z*se.ave.A/(1 - ave.A^2)
  ul0.A <- z.A + z*se.ave.A/(1 - ave.A^2)
  ll.A <- (exp(2*ll0.A) - 1)/(exp(2*ll0.A) + 1)
  ul.A <- (exp(2*ul0.A) - 1)/(exp(2*ul0.A) + 1)
  ave.B <- sum(g2*cor)/m2
  se.ave.B <- sqrt(sum(g2*var)/m2^2)
  z.B <- log((1 + ave.B)/(1 - ave.B))/2
  ll0.B <- z.B - z*se.ave.B/(1 - ave.B^2)
  ul0.B <- z.B + z*se.ave.B/(1 - ave.B^2)
  ll.B <- (exp(2*ll0.B) - 1)/(exp(2*ll0.B) + 1)
  ul.B <- (exp(2*ul0.B) - 1)/(exp(2*ul0.B) + 1)
  diff <- ave.A - ave.B
  se.diff <- sqrt(se.ave.A^2 + se.ave.B^2)
  ll.diff <- diff - sqrt((ave.A - ll.A)^2 + (ul.B - ave.B)^2)
  ul.diff <- diff + sqrt((ul.A - ave.A)^2 + (ave.B - ll.B)^2)
  out1 <- t(c(ave.A, se.ave.A, ll.A, ul.A))
  out2 <- t(c(ave.B, se.ave.B, ll.B, ul.B))
  out3 <- t(c(diff, se.diff, ll.diff, ul.diff))
  out <- rbind(out1, out2, out3)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Set A:", "Set B:", "Set A - Set B:")
  return (out)
}


#  meta.sub.spear =============================================
#' Confidence interval for a subgroup difference in average
#' Spearman correlations
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval 
#' for a  difference in average Spearman correlations for two 
#' mutually exclusive subgroups of studies. Each subgroup can have
#' one or more studies. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     n      	vector of sample sizes 
#' @param     cor    	vector of estimated Spearman correlations 
#' @param     group  	vector of group indicators:
#' * 1 for set A
#' * 2 for set B
#' * 0 to ignore
#' 
#' 
#' @return 
#' Returns a matrix with three rows:
#' * Row 1 - estimate for Set A
#' * Row 2 - estimate for Set B
#' * Row 3 - estimate for difference, Set A - Set B
#'
#' The columns are:
#' * Estimate - estimated average correlation or difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n <- c(55, 190, 65, 35)
#' cor <- c(.40, .65, .60, .45)
#' group <- c(1, 1, 2, 0)
#' meta.sub.spear(.05, n, cor, group)
#' 
#' # Should return:
#' #                Estimate         SE         LL        UL
#' # Set A:            0.525 0.06483629  0.3865928 0.6402793
#' # Set B:            0.600 0.08829277  0.3992493 0.7458512
#' # Set A - Set B:   -0.075 0.10954158 -0.2760700 0.1564955
#' 
#' 
#' @references
#' \insertRef{Bonett2008a}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.sub.spear <- function(alpha, n, cor, group) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  var <- (1 + cor^2/2)*(1 - cor^2)^2/(n - 3)
  g1 <- (group == rep(1, m))*1
  g2 <- (group == rep(2, m))*1
  m1 <- sum(g1)
  m2 <- sum(g2)
  ave.A <- sum(g1*cor)/m1
  se.ave.A <- sqrt(sum(g1*var)/m1^2)
  z.A <- log((1 + ave.A)/(1 - ave.A))/2
  ll0.A <- z.A - z*se.ave.A/(1 - ave.A^2)
  ul0.A <- z.A + z*se.ave.A/(1 - ave.A^2)
  ll.A <- (exp(2*ll0.A) - 1)/(exp(2*ll0.A) + 1)
  ul.A <- (exp(2*ul0.A) - 1)/(exp(2*ul0.A) + 1)
  ave.B <- sum(g2*cor)/m2
  se.ave.B <- sqrt(sum(g2*var)/m2^2)
  z.B <- log((1 + ave.B)/(1 - ave.B))/2
  ll0.B <- z.B - z*se.ave.B/(1 - ave.B^2)
  ul0.B <- z.B + z*se.ave.B/(1 - ave.B^2)
  ll.B <- (exp(2*ll0.B) - 1)/(exp(2*ll0.B) + 1)
  ul.B <- (exp(2*ul0.B) - 1)/(exp(2*ul0.B) + 1)
  diff <- ave.A - ave.B
  se.diff <- sqrt(se.ave.A^2 + se.ave.B^2)
  ll.diff <- diff - sqrt((ave.A - ll.A)^2 + (ul.B - ave.B)^2)
  ul.diff <- diff + sqrt((ul.A - ave.A)^2 + (ave.B - ll.B)^2)
  out1 <- t(c(ave.A, se.ave.A, ll.A, ul.A))
  out2 <- t(c(ave.B, se.ave.B, ll.B, ul.B))
  out3 <- t(c(diff, se.diff, ll.diff, ul.diff))
  out <- rbind(out1, out2, out3)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Set A:", "Set B:", "Set A - Set B:")
  return (out)
}


#  meta.sub.pbcor ===============================================
#' Confidence interval for a subgroup difference in average 
#' point-biserial correlations
#' 
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for  
#' a difference in average point-biserial correlations for two mutually
#' exclusive subgroups of studies. Each subgroup can have one or more 
#' studies. Two types of point-biserial correlations can be analyzed. 
#' One type uses an unweighted variance and is recommended for 2-group
#' experimental designs. The other type uses a weighted variance and 
#' is recommended for 2-group nonexperimental designs with simple random
#' sampling (but not stratified random sampling) within each study. 
#' Equality of variances within or across studies is not assumed.
#'    
#'    
#' @param     alpha   	alpha level for 1-alpha confidence
#' @param     m1    	  vector of estimated means for group 1 
#' @param     m2    	  vector of estimated means for group 2 
#' @param     sd1   	  vector of estimated SDs for group 1
#' @param     sd2   	  vector of estimated SDs for group 2
#' @param     n1    	  vector of group 1 sample sizes
#' @param     n2    	  vector of group 2 sample sizes
#' @param     type  	
#' * set to 1 for weighted variance
#' * set to 2 for unweighted variance
#' @param     group  	vector of group indicators:
#' * 1 for set A
#' * 2 for set B
#' * 0 to ignore
#' 
#' 
#' @return 
#' Returns a matrix with three rows:
#' * Row 1 - estimate for Set A
#' * Row 2 - estimate for Set B
#' * Row 3 - estimate for difference, Set A - Set B
#'
#' The columns are:
#' * Estimate - estimated average correlation or difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' m1 <- c(45.1, 39.2, 36.3, 34.5)
#' m2 <- c(30.0, 35.1, 35.3, 36.2)
#' sd1 <- c(10.7, 10.5, 9.4, 11.5)
#' sd2 <- c(12.3, 12.0, 10.4, 9.6)
#' n1 <- c(40, 20, 50, 25)
#' n2 <- c(40, 20, 48, 26)
#' group <- c(1, 1, 2, 2)
#' meta.sub.pbcor(.05,  m1, m2, sd1, sd2, n1, n2, 2, group)
#' 
#' # Should return:
#' #                   Estimate         SE         LL        UL
#' # Set A:          0.36338772 0.08552728  0.1854777 0.5182304
#' # Set B:         -0.01480511 0.08741322 -0.1840491 0.1552914
#' # Set A - Set B:  0.37819284 0.12229467  0.1320530 0.6075828
#' 
#' 
#' @references
#' \insertRef{Bonett2020b}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.sub.pbcor <- function(alpha,  m1, m2, sd1, sd2, n1, n2, type, group) {
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
    var <- (b^2*se.d^2)/(d^2 + b)^3
    cor <- d/sqrt(d^2 + b)
  }
  else {
    s <- sqrt((sd1^2 + sd2^2)/2)
    d <- (m1 - m2)/s
    a1 <- d^2*(sd1^4/df1 + sd1^4/df2)/(8*s^4)
    a2 <- sd1^2/(s^2*df1) + sd2^2/(s^2*df2)
    se.d <- sqrt(a1 + a2)
    var <- (16*se.d^2)/(d^2 + 4)^3
    cor <- d/sqrt(d^2 + 4)
  }
  g1 <- (group == rep(1, m))*1
  g2 <- (group == rep(2, m))*1
  m.A <- sum(g1)
  m.B <- sum(g2)
  ave.A <- sum(g1*cor)/m.A
  se.ave.A <- sqrt(sum(g1*var)/m.A^2)
  z.A <- log((1 + ave.A)/(1 - ave.A))/2
  ll0.A <- z.A - z*se.ave.A/(1 - ave.A^2)
  ul0.A <- z.A + z*se.ave.A/(1 - ave.A^2)
  ll.A <- (exp(2*ll0.A) - 1)/(exp(2*ll0.A) + 1)
  ul.A <- (exp(2*ul0.A) - 1)/(exp(2*ul0.A) + 1)
  ave.B <- sum(g2*cor)/m.B
  se.ave.B <- sqrt(sum(g2*var)/m.B^2)
  z.B <- log((1 + ave.B)/(1 - ave.B))/2
  ll0.B <- z.B - z*se.ave.B/(1 - ave.B^2)
  ul0.B <- z.B + z*se.ave.B/(1 - ave.B^2)
  ll.B <- (exp(2*ll0.B) - 1)/(exp(2*ll0.B) + 1)
  ul.B <- (exp(2*ul0.B) - 1)/(exp(2*ul0.B) + 1)
  diff <- ave.A - ave.B
  se.diff <- sqrt(se.ave.A^2 + se.ave.B^2)
  ll.diff <- diff - sqrt((ave.A - ll.A)^2 + (ul.B - ave.B)^2)
  ul.diff <- diff + sqrt((ul.A - ave.A)^2 + (ave.B - ll.B)^2)
  out1 <- t(c(ave.A, se.ave.A, ll.A, ul.A))
  out2 <- t(c(ave.B, se.ave.B, ll.B, ul.B))
  out3 <- t(c(diff, se.diff, ll.diff, ul.diff))
  out <- rbind(out1, out2, out3)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Set A:", "Set B:", "Set A - Set B:")
  return (out)
}


#  meta.sub.semipart =========================================
#' Confidence interval for a subgroup difference in average 
#' semipartial correlations
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval
#' for a difference in average semipartial correlations for two 
#' subgroups of mutually exclusive studies. Each subgroup can
#' have one or more studies. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     n      	vector of sample sizes 
#' @param     cor    	vector of estimated semi-partial correlations 
#' @param     r2   	  vector of squared multiple correlations for a model that
#' includes the IV and all control variables
#' @param     group  	vector of group indicators:
#' * 1 for set A
#' * 2 for set B
#' * 0 to ignore
#' 
#' 
#' @return 
#' Returns a matrix with three rows:
#' * Row 1 - estimate for Set A
#' * Row 2 - estimate for Set B
#' * Row 3 - estimate for difference, Set A - Set B
#'
#' The columns are:
#' * Estimate - estimated average correlation or difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' n <- c(55, 190, 65, 35)
#' cor <- c(.40, .65, .60, .45)
#' r2 <- c(.25, .41, .43, .39)
#' group <- c(1, 1, 2, 0)	
#' meta.sub.semipart(.05, n, cor, r2, group)
#' 
#' # Should return:
#' #                Estimate         SE         LL        UL
#' # Set A:            0.525 0.05955276  0.3986844 0.6317669
#' # Set B:            0.600 0.07931155  0.4221127 0.7333949
#' # Set A - Set B:   -0.075 0.09918091 -0.2587113 0.1324682
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.sub.semipart <- function(alpha, n, cor, r2, group) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  r0 <- r2 - cor^2
  var = (r2^2 - 2*r2 + r0 - r0^2 + 1)/(n - 3)
  g1 <- (group == rep(1, m))*1
  g2 <- (group == rep(2, m))*1
  m1 <- sum(g1)
  m2 <- sum(g2)
  ave.A <- sum(g1*cor)/m1
  se.ave.A <- sqrt(sum(g1*var)/m1^2)
  z.A <- log((1 + ave.A)/(1 - ave.A))/2
  ll0.A <- z.A - z*se.ave.A/(1 - ave.A^2)
  ul0.A <- z.A + z*se.ave.A/(1 - ave.A^2)
  ll.A <- (exp(2*ll0.A) - 1)/(exp(2*ll0.A) + 1)
  ul.A <- (exp(2*ul0.A) - 1)/(exp(2*ul0.A) + 1)
  ave.B <- sum(g2*cor)/m2
  se.ave.B <- sqrt(sum(g2*var)/m2^2)
  z.B <- log((1 + ave.B)/(1 - ave.B))/2
  ll0.B <- z.B - z*se.ave.B/(1 - ave.B^2)
  ul0.B <- z.B + z*se.ave.B/(1 - ave.B^2)
  ll.B <- (exp(2*ll0.B) - 1)/(exp(2*ll0.B) + 1)
  ul.B <- (exp(2*ul0.B) - 1)/(exp(2*ul0.B) + 1)
  diff <- ave.A - ave.B
  se.diff <- sqrt(se.ave.A^2 + se.ave.B^2)
  ll.diff <- diff - sqrt((ave.A - ll.A)^2 + (ul.B - ave.B)^2)
  ul.diff <- diff + sqrt((ul.A - ave.A)^2 + (ave.B - ll.B)^2)
  out1 <- t(c(ave.A, se.ave.A, ll.A, ul.A))
  out2 <- t(c(ave.B, se.ave.B, ll.B, ul.B))
  out3 <- t(c(diff, se.diff, ll.diff, ul.diff))
  out <- rbind(out1, out2, out3)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Set A:", "Set B:", "Set A - Set B:")
  return (out)
}


#  meta.sub.cronbach ==============================================
#' Confidence interval for a subgroup difference in average Cronbach 
#' reliabilities 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' difference in average Cronbach reliability coefficients for two 
#' mutually exclusive subgroups of studies. Each set can have one or
#' more studies. The number of measurements used to compute the sample
#' reliablity coefficient is assumed to be the same for all studies.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     n      	vector of sample sizes 
#' @param     rel    	vector of estimated Cronbach reliabilities 
#' @param     r      	number of measurements (e.g., items)
#' @param     group  	vector of group indicators:
#' * 1 for set A
#' * 2 for set B
#' * 0 to ignore
#' 
#' 
#' @return 
#' Returns a matrix with three rows:
#' * Row 1 - estimate for Set A
#' * Row 2 - estimate for Set B
#' * Row 3 - estimate for difference, Set A - Set B
#'
#' The columns are:
#' * Estimate - estimated average correlation or difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'        
#' @examples
#' n <- c(120, 170, 150, 135)
#' rel <- c(.89, .87, .73, .71)
#' group <- c(1, 1, 2, 2)
#' r <- 10
#' meta.sub.cronbach(.05, n, rel, r, group)
#' 
#' # Should return: 
#' #                Estimate         SE        LL        UL
#' # Set A:             0.88 0.01068845 0.8581268 0.8999386
#' # Set B:             0.72 0.02515130 0.6684484 0.7668524
#' # Set A - Set B:     0.16 0.02732821 0.1082933 0.2152731
#' 
#' 
#' @references
#' \insertRef{Bonett2010}{vcmeta}
#' \insertRef{Bonett2015b}{vcmeta]
#'
#'
#' @importFrom stats qnorm
#' @export
meta.sub.cronbach <- function(alpha, n, rel, r, group) {
  m <- length(n)
  z <- qnorm(1 - alpha/2)
  g1 <- (group == rep(1, m))*1
  g2 <- (group == rep(2, m))*1
  m1 <- sum(g1)
  m2 <- sum(g2)
  hn1 <- m1/sum(g1/n)
  hn2 <- m2/sum(g2/n)
  a1 <- ((r - 2)*(m1 - 1))^.25
  var1 <- 2*r*(1 - rel)^2/((r - 1)*(n - 2 - a1))
  a2 <- ((r - 2)*(m2 - 1))^.25
  var2 <- 2*r*(1 - rel)^2/((r - 1)*(n - 2 - a2))
  ave.A <- sum(g1*rel)/m1
  se.ave.A <- sqrt(sum(g1*var1)/m1^2)
  log.A <- log(1 - ave.A) - log(hn1/(hn1 - 1))
  ul.A <- 1 - exp(log.A - z*se.ave.A/(1 - ave.A))
  ll.A <- 1 - exp(log.A + z*se.ave.A/(1 - ave.A))
  ave.B <- sum(g2*rel)/m2
  se.ave.B <- sqrt(sum(g2*var2)/m2^2)
  log.B <- log(1 - ave.B) - log(hn2/(hn2 - 1))
  ul.B <- 1 - exp(log.B - z*se.ave.B/(1 - ave.B))
  ll.B <- 1 - exp(log.B + z*se.ave.B/(1 - ave.B))
  diff <- ave.A - ave.B
  se.diff <- sqrt(se.ave.A^2 + se.ave.B^2)
  ll.diff <- diff - sqrt((ave.A - ll.A)^2 + (ul.B - ave.B)^2)
  ul.diff <- diff + sqrt((ul.A - ave.A)^2 + (ave.B - ll.B)^2)
  out1 <- t(c(ave.A, se.ave.A, ll.A, ul.A))
  out2 <- t(c(ave.B, se.ave.B, ll.B, ul.B))
  out3 <- t(c(diff, se.diff, ll.diff, ul.diff))
  out <- rbind(out1, out2, out3)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Set A:", "Set B:", "Set A - Set B:")
  return (out)
}


#  meta.sub.gen ===============================================================
#' Confidence interval for a subgroup difference in average effect size
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' difference in the average effect size (any type of effect size) for
#' two mutually exclusive subgroups of studies. Each subgroup can have one
#' or more studies. All of the effects sizes should be compatible. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     est    	vector of estimated effect sizes 
#' @param     se    	vector of effect size standard errors 
#' @param     group  	vector of group indicators:
#' * 1 for set A
#' * 2 for set B
#' * 0 to ignore
#' 
#'
#' @return 
#' Returns a matrix with three rows:
#' * Row 1 - estimate for Set A
#' * Row 2 - estimate for Set B
#' * Row 3 - estimate for difference, Set A - Set B
#'
#' The columns are:
#' * Estimate - estimated average effect size or difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' est <- c(.920, .896, .760, .745)
#' se <- c(.098, .075, .069, .055) 
#' group <- c(1, 1, 2, 2)
#' meta.sub.gen(.05, est, se, group)
#' 
#' # Should return:
#' #                Estimate         SE          LL        UL
#' # Set A:           0.9080 0.06170292 0.787064504 1.0289355
#' # Set B:           0.7525 0.04411916 0.666028042 0.8389720
#' # Set A - Set B:   0.1555 0.07585348 0.006829917 0.3041701
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.sub.gen <- function(alpha, est, se, group) {
  m <- length(est)
  z <- qnorm(1 - alpha/2)
  var <- se^2
  g1 <- (group == rep(1, m))*1
  g2 <- (group == rep(2, m))*1
  m1 <- sum(g1)
  m2 <- sum(g2)
  ave.A <- sum(g1*est)/m1
  se.ave.A <- sqrt(sum(g1*var)/m1^2)
  ll.A <- ave.A - z*se.ave.A
  ul.A <- ave.A + z*se.ave.A
  ave.B <- sum(g2*est)/m2
  se.ave.B <- sqrt(sum(g2*var)/m2^2)
  ll.B <- ave.B - z*se.ave.B
  ul.B <- ave.B + z*se.ave.B
  diff <- ave.A - ave.B
  se.diff <- sqrt(se.ave.A^2 + se.ave.B^2)
  ll.diff <- diff - z*se.diff
  ul.diff <- diff + z*se.diff
  out1 <- t(c(ave.A, se.ave.A, ll.A, ul.A))
  out2 <- t(c(ave.B, se.ave.B, ll.B, ul.B))
  out3 <- t(c(diff, se.diff, ll.diff, ul.diff))
  out <- rbind(out1, out2, out3)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- c("Set A:", "Set B:", "Set A - Set B:")
  return (out)
}


# ================= Linear Contrasts of Effect Sizes ================
#  meta.lc.mean2 ====================================================
#' Confidence interval for a linear contrast of mean differences from
#' 2-group studies
#'
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' linear contrast of 2-group mean differences from two or more studies.
#' A Satterthwaite adjustment to the degrees of freedom is used to improve
#' the accuracy of the confidence interval. Equality of variances within or across
#' studies is not assumed. 
#'
#'
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    m1    	vector of estimated means for group 1
#' @param    m2    	vector of estimated means for group 2
#' @param    sd1   	vector of estimated SDs for group 1
#' @param    sd2   	vector of estimated SDs for group 2
#' @param    n1    	vector of group 1 sample sizes
#' @param    n2    	vector of group 2 sample sizes
#' @param    v     	vector of contrast coefficients
#'
#' @return 
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#'
#'
#' @examples
#' m1 <- c(45.1, 39.2, 36.3, 34.5)
#' m2 <- c(30.0, 35.1, 35.3, 36.2)
#' sd1 <- c(10.7, 10.5, 9.4, 11.5)
#' sd2 <- c(12.3, 12.0, 10.4, 9.6)
#' n1 <- c(40, 20, 50, 25)
#' n2 <- c(40, 20, 48, 26)
#' v <- c(.5, .5, -.5, -.5)
#' meta.lc.mean2(.05, m1, m2, sd1, sd2, n1, n2, v)
#'
#' # Should return:
#' #          Estimate       SE       LL       UL       df
#' # Contrast     9.95 2.837787 4.343938 15.55606 153.8362
#'
#'
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#'
#'
#' @importFrom stats qt
#' @export
meta.lc.mean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, v) {
  m <- length(m1)
  var1 <- sd1^2
  var2 <- sd2^2
  var <- var1/n1 + var2/n2
  d <- m1 - m2
  con <- t(v)%*%d
  var <- t(v)%*%(diag(var))%*%v
  se <- sqrt(var)
  u1 <- var^2*sum(v^2)^2
  u2 <- sum(v^4*var1^2/(n1^3 - n1^2) + v^4*var2^2/(n2^3 - n2^2))
  df <- u1/u2
  t <- qt(1 - alpha/2, df)
  ll <- con - t*se
  ul <- con + t*se
  out <- cbind(con, se, ll, ul, df)
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) = "Contrast" 
  return(out)
}


# meta.lc.stdmean2 ==================================================
#' Confidence interval for a linear contrast of standardized mean 
#' differences from 2-group studies  
#' 
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' linear contrast of 2-group standardized mean differences from two or 
#' more studies. Equality of variances within or across studies is not assumed. 
#' Use the square root average variance standardizer (stdzr = 0) for 2-group
#' experimental designs.  Use the square root weighted variance standardizer
#' (stdzr = 3) for 2-group nonexperimental designs with simple random sampling.
#' The stdzr = 1 and stdzr = 2 options can be used with either experimental
#' or nonexperimental designs.
#'
#'
#' @param    alpha  alpha level for 1-alpha confidence
#' @param    m1    	vector of estimated means for group 1 
#' @param    m2    	vector of estimated means for group 2 
#' @param    sd1   	vector of estimated SDs for group 1
#' @param    sd2   	vector of estimated SDs for group 2
#' @param    n1    	vector of group 1 sample sizes
#' @param    n2    	vector of group 2 sample sizes
#' @param    v     	vector of contrast coefficients
#' @param stdzr
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for group 1 SD standardizer 
#' * set to 2 for group 2 SD standardizer 
#' * set to 3 for square root weighted average variance standardizer
#' 
#' 
#' @return 
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples 
#' m1 <- c(45.1, 39.2, 36.3, 34.5)
#' m2 <- c(30.0, 35.1, 35.3, 36.2)
#' sd1 <- c(10.7, 10.5, 9.4, 11.5)
#' sd2 <- c(12.3, 12.0, 10.4, 9.6)
#' n1 <- c(40, 20, 50, 25)
#' n2 <- c(40, 20, 48, 26)
#' v <- c(.5, .5, -.5, -.5)
#' meta.lc.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, v, 0)
#' 
#' # Should return: 
#' #           Estimate        SE        LL       UL
#' # Contrast 0.8557914 0.2709192 0.3247995 1.386783
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.lc.stdmean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, v, stdzr) {
  df1 <- n1 - 1
  df2 <- n2 - 1
  m <- length(m1)
  z <- qnorm(1 - alpha/2)
  var1 <- sd1^2
  var2 <- sd2^2
  if (stdzr == 0) {
    s1 <- sqrt((var1 + var2)/2)
    d <- (m1 - m2)/s1
    du <- (1 - 3/(4*(n1 + n2) - 9))*d
    con <- t(v)%*%du
    var <- d^2*(var1^2/df1 + var2^2/df2)/(8*s1^4) + (var1/df1 + var2/df2)/s1^2
    se <- sqrt(t(v)%*%(diag(var))%*%v)
  } else if (stdzr == 1) { 
    d <- (m1 - m2)/sd1
    du <- (1 - 3/(4*n1 - 5))*d
    con <- t(v)%*%du
    var <- d^2/(2*df1) + 1/df1 + var2/(df2*var1)
    se <- sqrt(t(v)%*%(diag(var))%*%v)
  } else if (stdzr ==2) {
    d <- (m1 - m2)/sd2
    du <- (1 - 3/(4*n2 - 5))*d
    con <- t(v)%*%du
    var <- d^2/(2*df2) + 1/df2 + var1/(df1*var2)
    se <- sqrt(t(v)%*%(diag(var))%*%v)
  } else {
    s2 <- sqrt((df1*var1 + df2*var2)/(df1 + df2))
    d <- (m1 - m2)/s2
    du <- (1 - 3/(4*(n1 + n2) - 9))*d
    con <- t(v)%*%du
    var <- d^2*(1/df1 + 1/df2)/8 + (var1/n1 + var2/n2)/s2^2
    se <- sqrt(t(v)%*%(diag(var))%*%v)
  }
  ll <- con - z*se
  ul <- con + z*se
  out <- cbind(con, se, ll, ul)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- "Contrast"
  return(out)
}


#  meta.lc.mean.ps ====================================================
#' Confidence interval for a linear contrast of mean differences from 
#' paired-samples studies 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' linear contrast of paired-samples mean differences from two or more studies.
#' A Satterthwaite adjustment to the degrees of freedom is used to improve
#' the accuracy of the confidence interval. Equality of variances within or across
#' studies is not assumed. 
#'
#'
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    m1    	vector of estimated means for measurement 1 
#' @param    m2    	vector of estimated means for measurement 2 
#' @param    sd1   	vector of estimated SDs for measurement 1
#' @param    sd2   	vector of estimated SDs for measurement 2
#' @param    cor   	vector of estimated correlations for paired measurements
#' @param    n     	vector of sample sizes
#' @param    v     	vector of contrast coefficients
#' 
#' 
#' @return 
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @examples
#' m1 <- c(53, 60, 53, 57)
#' m2 <- c(55, 62, 58, 61)
#' sd1 <- c(4.1, 4.2, 4.5, 4.0)
#' sd2 <- c(4.2, 4.7, 4.9, 4.8)
#' cor <- c(.72, .78, .81, .85)
#' n <- c(30, 50, 30, 70)
#' v <- c(.5, .5, -.5, -.5)
#' meta.lc.mean.ps(.05, m1, m2, sd1, sd2, cor, n, v)
#' 
#' # Should return:
#' #          Estimate        SE      LL      UL      df
#' # Contrast      2.5 0.4681205 1.57207 3.42793 107.657
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @export
meta.lc.mean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, v) {
  m <- length(m1)
  var1 <- sd1^2
  var2 <- sd2^2
  var <- (var1 + var2 - 2*cor*sd1*sd2)/n
  d <- m1 - m2
  con <- t(v)%*%d
  se <- sqrt(t(v)%*%(diag(var))%*%v)
  u1 <- sum(var*v^2)^2
  u2 <- sum((var*v^2)^2/(n - 1))
  df <- u1/u2
  t <- qt(1 - alpha/2, df)
  ll <- con - t*se
  ul <- con + t*se
  out <- cbind(con, se, ll, ul, df)
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) = "Contrast" 
  return(out)
}


#  meta.lc.stdmean.ps ===================================================
#' Confidence interval for a linear contrast of standardized 
#' mean differences from paired-samples studies 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' linear contrast of paired-samples standardized mean differences from two or 
#' more studies. Equality of variances within or across studies is not assumed. 
#'
#'
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    m1    	vector of estimated means for measurement 1 
#' @param    m2    	vector of estimated means for measurement 2 
#' @param    sd1   	vector of estimated SDs for measurement 1
#' @param    sd2   	vector of estimated SDs for measurement 2
#' @param    cor	  vector of estimated correlations for paired measurements
#' @param    n     	vector of sample sizes
#' @param    v     	vector of contrast coefficients
#' @param stdzr
#' * set to 0 for square root unweighted average variance standardizer 
#' * set to 1 for measurement 1 SD standardizer 
#' * set to 2 for measurement 2 SD standardizer 
#' 
#' 
#' @return 
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' m1 <- c(53, 60, 53, 57)
#' m2 <- c(55, 62, 58, 61)
#' sd1 <- c(4.1, 4.2, 4.5, 4.0)
#' sd2 <- c(4.2, 4.7, 4.9, 4.8)
#' cor <- c(.72, .78, .81, .85)
#' n <- c(30, 50, 30, 70)
#' v <- c(.5, .5, -.5, -.5)
#' meta.lc.stdmean.ps(.05, m1, m2, sd1, sd2, cor, n, v, 0)
#' 
#' # Should return:
#' #            Estimate        SE        LL        UL
#' # Contrast  0.5127577 0.1346794 0.2487908 0.7767245
#' 
#' 
#' @references
#' \insertRef{Bonett2009a}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.lc.stdmean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, v, stdzr) {
  df <- n - 1
  m <- length(m1)
  z <- qnorm(1 - alpha/2)
  var1 <- sd1^2
  var2 <- sd2^2
  vd <- var1 + var2 - 2*cor*sd1*sd2
  if (stdzr == 0) {
    s <- sqrt((var1 + var2)/2)
    d <- (m1 - m2)/s
    du <- sqrt((n - 2)/df)*d
    var <- d^2*(var1^2 + var2^2 + 2*cor^2*var1*var2)/(8*df*s^4) + vd/(df*s^2)
    con <- t(v)%*%du
    se <- sqrt(t(v)%*%(diag(var))%*%v)
  } else if (stdzr == 1) { 
    d <- (m1 - m2)/sd1
    du <- (1 - 3/(4*df - 1))*d
    con <- t(v)%*%du
    var <- d^2/(2*df) + vd/(df*var1)
    se <- sqrt(t(v)%*%(diag(var))%*%v)
  } else {
    d <- (m1 - m2)/sd2
    du <- (1 - 3/(4*df - 1))*d
    con <- t(v)%*%du
    var <- d^2/(2*df) + vd/(df*var2)
    se <- sqrt(t(v)%*%(diag(var))%*%v)
  }
  ll <- con - z*se
  ul <- con + z*se
  out <- cbind(con, se, ll, ul)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- "Contrast"
  return(out)
}


#  meta.lc.meanratio2 =================================================
#' Confidence interval for a log-linear contrast of mean ratios from
#' 2-group studies  	
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' log-linear contrast of 2-group mean ratios from two or more studies. A
#' Satterthwaite adjustment to the degrees of freedom is used to improve
#' the accuracy of the confidence interval. Equality of variances within or across
#' studies is not assumed. 
#'
#'
#' @param    alpha 	alpha level for 1-alpha confidence
#' @param    m1    	vector of estimated means for group 1 
#' @param    m2    	vector of estimated means for group 2 
#' @param    sd1   	vector of estimated SDs for group 1
#' @param    sd2	  vector of estimated SDs for group 2
#' @param    n1    	vector of group 1 sample sizes
#' @param    n2    	vector of group 2 sample sizes
#' @param    v     	vector of contrast coefficients
#' 
#' 
#' @return 
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated log-linear contrast
#' * SE - standard error of log-linear contrast
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * exp(Estimate) - exponentiated log-linear contrast
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @examples
#' m1 <- c(45.1, 39.2, 36.3, 34.5)
#' m2 <- c(30.0, 35.1, 35.3, 36.2)
#' sd1 <- c(10.7, 10.5, 9.4, 11.5)
#' sd2 <- c(12.3, 12.0, 10.4, 9.6)
#' n1 <- c(40, 20, 50, 25)
#' n2 <- c(40, 20, 48, 26)
#' v <- c(.5, .5, -.5, -.5)
#' meta.lc.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, v)
#' 
#' # Should return:
#' #           Estimate         SE        LL        UL  exp(Estimate)
#' # Contrast 0.2691627 0.07959269 0.1119191 0.4264064      1.308868
#' #           exp(LL)  exp(UL)       df
#' # Contrast 1.118422 1.531743 152.8665
#' 
#' 
#' @references
#' \insertRef{Bonett2020}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @export
meta.lc.meanratio2 <- function(alpha, m1, m2, sd1, sd2, n1, n2, v) {
  logratio <- log(m1/m2)
  var1 <- sd1^2/(n1*m1^2) 
  var2 <- sd2^2/(n2*m2^2)
  est <- t(v)%*%logratio
  se <- sqrt(t(v)%*%(diag(var1 + var2))%*%v)
  df <- se^4/sum(v^4*var1^2/(n1 - 1) + v^4*var2^2/(n2 - 1))
  t <- qt(1 - alpha/2, df)
  ll <- est - t*se
  ul <- est + t*se
  out <- cbind(est, se, ll, ul, exp(est), exp(ll), exp(ul), df)
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "exp(Estimate)", 
                     "exp(LL)", "exp(UL)", "df")
  rownames(out) <- "Contrast"
  return(out)
}


#  meta.lc.meanratio.ps ================================================
#' Confidence interval for a log-linear contrast of mean ratios from 
#' paired-samples studies 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' log-linear contrast of paired-sample mean ratios from two or more studies.
#' A Satterthwaite adjustment to the degrees of freedom is used to improve
#' the accuracy of the confidence interval. Equality of variances within or across
#' studies is not assumed. 
#'
#'
#' @param    alpha	alpha level for 1-alpha confidence
#' @param    m1    	vector of estimated means for measurement 1 
#' @param    m2    	vector of estimated means for measurement 2 
#' @param    sd1   	vector of estimated SDs for measurement 1
#' @param    sd2   	vector of estimated SDs for measurement 2
#' @param    cor   	vector of estimated correlations for paired measurements
#' @param    n     	vector of sample sizes
#' @param    v     	vector of contrast coefficients
#' 
#' 
#' @return
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated log-linear contrast
#' * SE - standard error of log-linear contrast
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * exp(Estimate) - exponentiated log-linear contrast
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @examples
#' m1 <- c(53, 60, 53, 57)
#' m2 <- c(55, 62, 58, 61)
#' sd1 <- c(4.1, 4.2, 4.5, 4.0)
#' sd2 <- c(4.2, 4.7, 4.9, 4.8)
#' cor <- c(.72, .78, .81, .85)
#' n <- c(30, 50, 30, 70)
#' v <- c(.5, .5, -.5, -.5)
#' meta.lc.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, v)
#' 
#' # Should return:
#' #           Estimate       SE         LL         UL exp(Estimate)
#' # Contrast 0.0440713 0.008265 0.02767047 0.06047213      1.045057
#' #           exp(LL)  exp(UL)       df
#' # Contrast 1.028057 1.062338 98.38086
#'
#' 
#' @references
#' \insertRef{Bonett2020}{vcmeta}
#'
#'
#' @importFrom stats qt
#' @export
meta.lc.meanratio.ps <- function(alpha, m1, m2, sd1, sd2, cor, n, v) {
  logratio <- log(m1/m2)
  var <- (sd1^2/m1^2 + sd2^2/m2^2 - 2*cor*sd1*sd2/(m1*m2))/n
  est <- t(v)%*%logratio
  se <- sqrt(t(v)%*%(diag(var))%*%v)
  df <- se^4/sum(v^4*var^2/(n - 1))
  t <- qt(1 - alpha/2, df)
  ll <- est - t*se
  ul <- est + t*se
  out <- cbind(est, se, ll, ul, exp(est), exp(ll), exp(ul), df)
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "exp(Estimate)", 
                    "exp(LL)", "exp(UL)", "df")
  rownames(out) <- "Contrast"
  return(out)
}


#  meta.lc.oddsratio ==================================================
#' Confidence interval for a log-linear contrast of odds ratios 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' exponentiated log-linear contrast of odds ratios from two or more studies.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f1     	vector of group 1 frequency counts
#' @param     f2     	vector of group 2 frequency counts
#' @param     n1     	vector of group 1 sample sizes 
#' @param     n2     	vector of group 2 sample sizes
#' @param     v      	vector of contrast coefficients
#' 
#' 
#' @return 
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated log-linear contrast
#' * SE - standard error of log-linear contrast
#' * exp(Estimate) - exponentiated log-linear contrast
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' n1 <- c(50, 150, 150)
#' f1 <- c(16, 50, 25)
#' n2 <- c(50, 150, 150)
#' f2 <- c(7, 15, 20)
#' v <- c(1, -1, 0)
#' meta.lc.oddsratio(.05, f1, f2, n1, n2, v)
#' 
#' # Should return:
#' #            Estimate        SE  exp(Estimate)   exp(LL)  exp(UL)
#' # Contrast -0.4596883 0.5895438      0.6314805 0.1988563 2.005305
#' 
#' 
#' @references
#' \insertRef{Bonett2015}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.lc.oddsratio <- function(alpha, f1, f2, n1, n2, v) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  lor <- log((f1 + .5)*(n2 - f2 + .5)/((f2 + .5)*(n1 - f1 + .5)))
  var.lor <- 1/(f1 + .5) + 1/(f2 + .5) + 1/(n1 - f1 + .5) + 1/(n2 - f2 + .5)
  con.lor <- t(v)%*%lor
  se.lor <- sqrt(t(v)%*%(diag(var.lor))%*%v)
  ll <- exp(con.lor - z*se.lor)
  ul <- exp(con.lor + z*se.lor)
  con.or <- exp(con.lor)
  se.or <- con.or*se.lor
  out <- cbind(con.lor, se.lor, con.or, ll, ul)
  colnames(out) <- c("Estimate", "SE", "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- "Contrast"
  return (out)
}


#  meta.lc.propratio2 =====================================================
#' Confidence interval for a log-linear contrast of proportion ratios from 2-group studies
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for an 
#' exponentiated log-linear contrast of 2-group proportion ratios from
#' two or more studies.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f1     	vector of group 1 frequency counts
#' @param     f2     	vector of group 2 frequency counts
#' @param     n1     	vector of group 1 sample sizes 
#' @param     n2     	vector of group 2 sample sizes
#' @param     v      	vector of contrast coefficients
#' 
#' @return 
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated log-linear contrast
#' * SE - standard error of log-linear contrast
#' * exp(Estimate) - exponentiated log-linear contrast
#' * exp(LL) - lower limit of the exponentiated confidence interval
#' * exp(UL) - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' n1 <- c(50, 150, 150)
#' f1 <- c(16, 50, 25)
#' n2 <- c(50, 150, 150)
#' f2 <- c(7, 15, 20)
#' v <- c(1, -1, 0)
#' meta.lc.propratio2(.05, f1, f2, n1, n2, v)
#'
#' # Should return:
#' #            Estimate        SE  exp(Estimate)   exp(LL)  exp(UL)
#' # Contrast -0.3853396 0.4828218      0.6802196 0.2640405 1.752378
#' 
#' 
#' @references
#' \insertRef{Price2008}{vcmeta}
#' \insertRef{Bonett2015}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.lc.propratio2 <- function(alpha, f1, f2, n1, n2, v) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  p1 <- (f1 + 1/4)/(n1 + 7/4) 
  p2 <- (f2 + 1/4)/(n2 + 7/4)
  lrr <- log(p1/p2)
  var1 <- 1/(f1 + 1/4 + (f1 + 1/4)^2/(n1 - f1 + 3/2))
  var2 <- 1/(f2 + 1/4 + (f2 + 1/4)^2/(n2 - f2 + 3/2))
  var.lrr <- var1 + var2
  con.lrr <- t(v)%*%lrr
  se.lrr = sqrt(t(v)%*%(diag(var.lrr))%*%v)
  ll <- exp(con.lrr - z*se.lrr)
  ul <- exp(con.lrr + z*se.lrr)
  con.rr <- exp(con.lrr)
  se.rr <- con.rr*se.lrr
  out <- cbind(con.lrr, se.lrr, con.rr, ll, ul)
  colnames(out) <- c("Estimate", "SE", "exp(Estimate)", "exp(LL)", "exp(UL)")
  rownames(out) <- "Contrast"
  return (out)
} 


#  meta.lc.prop2 =============================================================
#' Confidence interval for a linear contrast of proportion differences in 2-group studies 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and adjusted Wald confidence interval for a 
#' linear contrast of 2-group proportion differences from two or more studies.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f1     	vector of group 1 frequency counts
#' @param     f2     	vector of group 2 frequency counts
#' @param     n1     	vector of group 1 sample sizes 
#' @param     n2     	vector of group 2 sample sizes
#' @param     v      	vector of contrast coefficients
#' 
#'
#' @return
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#' 
#' 
#' @examples
#' n1 <- c(50, 150, 150)
#' n2 <- c(50, 150, 150)
#' f1 <- c(16, 50, 25)
#' f2 <- c(7, 15, 20)
#' v <- c(1, -1, 0)
#' meta.lc.prop2(.05, f1, f2, n1, n2, v)
#' 
#' # Should return:
#' #             Estimate         SE         LL        UL
#' # Contrast -0.05466931 0.09401019 -0.2389259 0.1295873
#' 
#'
#' @references 
#' \insertRef{Bonett2014}{vcmeta}
#'
#'
#' @importFrom stats qnorm
#' @export
meta.lc.prop2 <- function(alpha, f1, f2, n1, n2, v) {
  m <- length(n1)
  z <- qnorm(1 - alpha/2)
  p1 <- (f1 + 1/m)/(n1 + 2/m) 
  p2 <- (f2 + 1/m)/(n2 + 2/m)
  rd <- p1 - p2
  var1 <- p1*(1 - p1)/(n1 + 2/m)
  var2 <- p2*(1 - p2)/(n2 + 2/m) 
  var <- var1 + var2
  con <- t(v)%*%rd
  se <- sqrt(t(v)%*%(diag(var))%*%v)
  ll <- con - z*se
  ul <- con + z*se
  out <- cbind(con, se, ll, ul)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) = "Contrast" 
  return (out)
}


# meta.lc.prop.ps =======================================================
#' Confidence interval for a linear contrast of proportion differences in 
#' paired-samples studies 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and adjusted Wald confidence interval
#' for a linear contrast of paired-samples proportion differences from two or
#' more studies.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f11    	vector of frequency counts in cell 1,1
#' @param     f12    	vector of frequency counts in cell 1,2
#' @param     f21    	vector of frequency counts in cell 2,1
#' @param     f22    	vector of frequency counts in cell 2,2
#' @param     v      	vector of contrast coefficients
#' 
#' 
#' @return
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#' 
#' 
#' @examples
#' f11 <- c(17, 28, 19)
#' f12 <- c(43, 56, 49)
#' f21 <- c(3, 5, 5)
#' f22 <- c(37, 54, 39)
#' v <- c(.5, .5, -1)
#' meta.lc.prop.ps(.05, f11, f12, f21, f22, v)
#'
#' # Should return:
#' #              Estimate         SE         LL       UL
#' #  Contrast -0.01436285 0.06511285 -0.1419817 0.113256
#' 
#' 
#' @references
#' \insertRef{Bonett2012}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.lc.prop.ps <- function(alpha, f11, f12, f21, f22, v) {
  m <- length(f11)
  z <- qnorm(1 - alpha/2)
  n <- f11 + f12 + f21 + f22
  p12 <- (f12 + 1/m)/(n + 2/m) 
  p21 <- (f21 + 1/m)/(n + 2/m)
  rd <- p12 - p21
  con <- t(v)%*%rd
  var <- (p12 + p21 - rd^2)/(n + 2/m)
  se <- sqrt(t(v)%*%(diag(var))%*%v)
  ll <- con - z*se
  ul <- con + z*se
  out <- cbind(con, se, ll, ul)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- "Contrast"
  return (out)
}


#  meta.lc.agree ===================================================
#' Confidence interval for a linear contrast of G-index coefficients  
#' 
#' 
#' @description
#' Computes the estimate, standard error, and adjusted Wald confidence 
#' interval for a linear contrast of G-index of agreement coefficients 
#' from two or more studies. This function assumes that two raters each
#' provide a dichotomous rating for a sample of objects.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f11    	vector of frequency counts in cell 1,1
#' @param     f12    	vector of frequency counts in cell 1,2
#' @param     f21    	vector of frequency counts in cell 2,1
#' @param     f22    	vector of frequency counts in cell 2,2
#' @param     v      	vector of contrast coefficients
#' 
#' @return
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#' 
#' 
#' @examples
#' f11 <- c(43, 56, 49)
#' f12 <- c(7, 2, 9)
#' f21 <- c(3, 5, 5)
#' f22 <- c(37, 54, 39)
#' v <- c(.5, .5, -1)
#' meta.lc.agree(.05, f11, f12, f21, f22, v)
#' 
#' # Should return:
#' #            Estimate        SE         LL        UL
#' # Contrast 0.1022939 0.07972357 -0.05396142 0.2585492
#' 
#' 
#' @references 
#' \insertRef{Bonett2022}{vcmeta}
#' 
#' 
#' @importFrom stats qnorm
#' @export
meta.lc.agree <- function(alpha, f11, f12, f21, f22, v) {
  m <- length(f11)
  z <- qnorm(1 - alpha/2)
  n <- f11 + f12 + f21 + f22
  p0 <- (f11 + f22 + 2/m)/(n + 4/m)
  g <- 2*p0 - 1 
  con <- t(v)%*%g
  var <- 4*p0*(1 - p0)/(n + 4/m)
  se <- sqrt(t(v)%*%(diag(var))%*%v)
  ll <- con - z*se
  ul <- con + z*se
  out <- cbind(con, se, ll, ul)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- "Contrast"
  return (out)
}


# meta.lc.mean1 ===================================================
#' Confidence interval for a linear contrast of means
#' 
#'
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' linear contrast of means from two or more studies. This function will
#' use either an unequal variance (recommended) or an equal variance method. 
#' A Satterthwaite adjustment to the degrees of freedom is used with the
#' unequal variance method. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     m     	vector of estimated means
#' @param     sd    	vector of estimated standard deviations
#' @param     n     	vector of sample sizes
#' @param     v     	vector of contrast coefficients
#' @param     eqvar 	
#' * FALSE for unequal variance method
#' * TRUE for equal variance method
#' 
#' @return
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' * df - degrees of freedom
#' 
#' 
#' @examples
#' m <- c(33.5, 37.9, 38.0, 44.1)
#' sd <- c(3.84, 3.84, 3.65, 4.98)
#' n <- c(10, 10, 10, 10)
#' v <- c(.5, .5, -.5, -.5)
#' meta.lc.mean1(.05, m, sd, n, v, eqvar = FALSE)
#'
#' # Should return:
#' #          Estimate       SE        LL        UL       df
#' # Contrast    -5.35 1.300136 -7.993583 -2.706417 33.52169
#' 
#' @references
#' \insertRef{Snedecor1980}{vcmeta}
#' 
#' 
#' @importFrom stats qt
#' @export
meta.lc.mean1 <- function(alpha, m, sd, n, v, eqvar = FALSE) {
  est <- t(v)%*%m 
  k <- length(m)
  if (eqvar){
    df <- sum(n) - k
    v1 <- sum((n - 1)*sd^2)/df
    se <- sqrt(v1*t(v)%*%solve(diag(n))%*%v)
    t1 <- qt(1 - alpha/2, df)
    ll <- est - t1*se
    ul <- est + t1*se
  } else {
    v2 <- diag(sd^2)%*%(solve(diag(n)))
    se <- sqrt(t(v)%*%v2%*%v)
    df = (se^4)/sum(((v^4)*(sd^4)/(n^2*(n-1))))
    t2 <- qt(1 - alpha/2, df)
    ll <- est - t2*se
    ul <- est + t2*se
  }
  out <- cbind(est, se, ll, ul, df)
  colnames(out) <- c("Estimate", "SE", "LL", "UL", "df")
  rownames(out) <- "Contrast"
  return(out)
}


# meta.lc.prop1 ==============================================
#' Confidence interval for a linear contrast of proportions 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and an adjusted Wald confidence 
#' interval for a linear contrast of proportions from two or more studies.
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     f      	vector of frequency counts
#' @param     n      	vector of sample sizes
#' @param     v       vector of contrast coefficients
#' 
#'
#' @return
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#' @examples
#' f <- c(26, 24, 38)
#' n <- c(60, 60, 60)
#' v <- c(-.5, -.5, 1)
#' meta.lc.prop1(.05, f, n, v)
#'
#' # Should return: 
#' #           Estimate         SE         LL        UL
#' # Contrast 0.2119565 0.07602892 0.06294259 0.3609705
#' 
#' 
#' @references
#' \insertRef{Price2004}{vcmeta}
#'
#'
#' @importFrom stats qnorm
#' @export
meta.lc.prop1 <- function(alpha, f, n, v) {
  z <- qnorm(1 - alpha/2)
  m <- length(v) - length(which(v==0))
  p <- (f + 2/m)/(n + 4/m)
  est <- t(v)%*%p
  se <- sqrt(t(v)%*%diag(p*(1 - p))%*%solve(diag(n + 4/m))%*%v)
  ll <- est - z*se
  ul <- est + z*se
  out <- cbind(est, se, ll, ul)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- "Contrast"
  return(out)
}


#  meta.lc.gen =====================================================
#' Confidence interval for a linear contrast of effect sizes 
#' 
#' 
#' @description
#' Computes the estimate, standard error, and confidence interval for a 
#' linear contrast of any type of effect size from two or more studies.
#'
#'
#' @param     alpha 	alpha level for 1-alpha confidence
#' @param     est   	vector of parameter estimates
#' @param     se    	vector of standard errors
#' @param     v     	vector of contrast coefficients
#' 
#' 
#' @return
#' Returns 1-row matrix with the following columns: 
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#'   
#' @examples
#' est <- c(.55, .59, .44, .48, .26, .19)
#' se <- c(.054, .098, .029, .084, .104, .065)
#' v <- c(.5, .5, -.25, -.25, -.25, -.25)
#' meta.lc.gen(.05, est, se, v)
#'
#' # Should return: 
#' #          Estimate         SE        LL        UL
#' # Contrast   0.2275 0.06755461 0.0950954 0.3599046
#' 
#' @importFrom stats qnorm
#' @export
meta.lc.gen <- function(alpha, est, se, v) {
  z <- qnorm(1 - alpha/2)
  con <- t(v)%*%est 
  se <- sqrt(t(v)%*%(diag(se^2))%*%v)
  ll <- con - z*se
  ul <- con + z*se
  out <- cbind(con, se, ll, ul)
  colnames(out) <- c("Estimate", "SE", "LL", "UL")
  rownames(out) <- "Contrast"
  return(out)
}

