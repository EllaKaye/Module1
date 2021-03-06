% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Module1_code.R
\name{BLB.1d}
\alias{BLB.1d}
\title{Bag of Little Bootstraps in 1d}
\usage{
BLB.1d(data, gamma, FUN, ..., s = 20, r = 100)
}
\arguments{
\item{data}{a vector of the original sample}

\item{gamma}{controls the size of the subsamples - ideally in [0.5, 1]}

\item{FUN}{a function to calculate the estimator of the parameter of interest}

\item{s}{the number of subsamples}

\item{r}{the number of bootstrap replications (resamples) per subsample}
}
\value{
the mean across the s subsamples of the standard errors of the parameter estimates
}
\description{
Implements Bag of Little Bootstraps on a 1-d data vector, with a generic function for the estimator
}
\examples{
X <- rnorm(5000)
BLB.1d(X, mean, gamma=0.5)
}

