% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_ave.R
\name{meta.ave.slope}
\alias{meta.ave.slope}
\title{Confidence interval for an average slope coefficient}
\usage{
meta.ave.slope(alpha, n, cor, sdy, sdx, bystudy = TRUE)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{n}{vector of sample sizes}

\item{cor}{vector of estimated correlations}

\item{sdy}{vector of estimated SDs of y}

\item{sdx}{vector of estimated SDs of x}

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
\item df - degrees of freedom
}
}
\description{
Computes the estimate, standard error, and confidence interval for an
average slope coefficient in a simple linear regression model from two
or more studies. A Satterthwaite adjustment to the degrees of freedom
is used to improve the accuracy of the confidence interval.
}
\examples{
n <- c(45, 85, 50, 60)
cor <- c(.24, .35, .16, .20)
sdy <- c(12.2, 14.1, 11.7, 15.9)
sdx <- c(1.34, 1.87, 2.02, 2.37)
meta.ave.slope(.05, n, cor, sdy, sdx, bystudy = TRUE)

# Should return:
#          Estimate        SE         LL       UL       df
# Average 1.7731542 0.4755417  0.8335021 2.712806 149.4777
# Study 1 2.1850746 1.3084468 -0.4536599 4.823809  43.0000
# Study 2 2.6390374 0.7262491  1.1945573 4.083518  83.0000
# Study 3 0.9267327 0.8146126 -0.7111558 2.564621  48.0000
# Study 4 1.3417722 0.8456799 -0.3510401 3.034584  58.0000


}
