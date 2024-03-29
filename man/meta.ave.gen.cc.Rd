% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_ave.R
\name{meta.ave.gen.cc}
\alias{meta.ave.gen.cc}
\title{Confidence interval for an average effect size using a constant
coefficient model}
\usage{
meta.ave.gen.cc(alpha, est, se, bystudy = TRUE)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{est}{vector of parameter estimates}

\item{se}{vector of standard errors}

\item{bystudy}{logical to also return each study estimate (TRUE) or not}
}
\value{
Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
\itemize{
\item Estimate - estimated effect size
\item SE - standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes the estimate, standard error, and confidence interval for a
weighted average effect from two or more studies using the constant
coefficient (fixed-effect) meta-analysis model.
}
\details{
The weighted average estimate will be biased regardless of the number of
studies or the sample size in each study. The actual confidence interval
coverage probability can be much smaller than the specified confidence
level when the population effect sizes are not identical across studies.

The constant coefficient model should be used with caution, and the varying
coefficient methods in this package are the recommended alternatives. The
varying coefficient methods do not require effect-size homogeneity across
the selected studies. This constant coefficient meta-analysis function is
included in the vcmeta package primarily for classroom  demonstrations to
illustrate the problematic characteristics of the constant coefficient
meta-analysis model.
}
\examples{
est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
meta.ave.gen.cc(.05, est, se, bystudy = TRUE)

# Should return:
#             Estimate         SE          LL        UL
# Average    0.3127916 0.06854394  0.17844794 0.4471352
# Study 1    0.0220000 0.12400000 -0.22103553 0.2650355
# Study 2    0.7510000 0.46400000 -0.15842329 1.6604233
# Study 3    0.4210000 0.10200000  0.22108367 0.6209163
# Study 4    0.2870000 0.59200000 -0.87329868 1.4472987
# Study 5    0.0520000 0.86400000 -1.64140888 1.7454089
# Study 6    0.1460000 0.24100000 -0.32635132 0.6183513
# Study 7    0.5620000 0.25200000  0.06808908 1.0559109
# Study 8    0.9040000 0.31800000  0.28073145 1.5272685


}
\references{
\itemize{
\item \insertRef{Hedges1985}{vcmeta}
\item \insertRef{Borenstein2009}{vcmeta}
}
}
\seealso{
\link[vcmeta]{meta.ave.gen}
}
