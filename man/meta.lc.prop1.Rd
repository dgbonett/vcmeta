% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_comp.R
\name{meta.lc.prop1}
\alias{meta.lc.prop1}
\title{Confidence interval for a linear contrast of proportions}
\usage{
meta.lc.prop1(alpha, f, n, v)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{f}{vector of frequency counts}

\item{n}{vector of sample sizes}

\item{v}{vector of contrast coefficients}
}
\value{
Returns 1-row matrix with the following columns:
\itemize{
\item Estimate -estimated linear contrast
\item SE - standard error
\item LL - lower limit of the adjusted Wald confidence interval
\item UL - upper limit of the adjusted Wald confidence interval
}
}
\description{
Computes the estimate, standard error, and an adjusted Wald confidence
interval for a linear contrast of proportions from two or more studies.
}
\examples{
f <- c(26, 24, 38)
n <- c(60, 60, 60)
v <- c(-.5, -.5, 1)
meta.lc.prop1(.05, f, n, v)

# Should return: 
#           Estimate         SE         LL        UL
# Contrast 0.2119565 0.07602892 0.06294259 0.3609705


}
\references{
\insertRef{Price2004}{vcmeta}
}
