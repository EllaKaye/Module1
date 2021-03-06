% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Module1_code.R
\name{BLB.adapt}
\alias{BLB.adapt}
\title{BLB with adaptive selection on r and s}
\usage{
BLB.adapt(data, gamma = 0.7, w_s = 3, w_r = 20, lambda = 10^(-5),
  epsilon = 0.05, alpha = 0.05)
}
\arguments{
\item{data}{a matrix or dataframe}

\item{gamma}{specifies the subsample size by b=n^gamma}

\item{w_s}{window size for the adaptive selection of s}

\item{w_r}{window size for the adaptive selection of r}

\item{lambda}{specifies the L2 penalty in Ridge regression}

\item{alpha}{the level of the confidence which is set to 0.05}

\item{epsilon}{relative error for the convergence criterion in adaptive selection}
}
\value{
a list consisting of
\item{s}{the total number of subsamples from the original dataset}
\item{r}{the number of resamples for each subsample}
\item{final_CI}{the final confidence interval for all the components of the parameter estimate }
\item{CI_widths}{the widths of the final confidence interval for all the components of the parameter estimate }
\item{mean_width}{mean_width the average marginal width of the confidence intervals across dimensions}
\item{final_se}{the final standard error for all the componennts of the parameter estimate}
\item{mean_se}{mean_se the average standard error across the parameter dimensions}
}
\description{
BLB function which returns the average (across dimensions) confidence interal width or the average standard error of the parameter and the values of
r and s which were chosen adaptively within the function
}
\examples{
X <- mvrnorm(100, mu = rep(0,2), Sigma <- diag(2))
epsilon <- rnorm(100, mean=0, sd=sqrt(10))
t_theta <- as.matrix(rep(1, 2))
Y <- X \%*\% t_theta + epsilon
my_data <- cbind(X, Y)
my_result <- BLB.adapt(my_data, gamma=0.7, w_s=3, w_r=20, lambda=10^(-5), epsilon=0.05, alpha=0.05)
}

