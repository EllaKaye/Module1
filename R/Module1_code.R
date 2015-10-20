# libraries
library(ggplot2)
library(boot)

### Bootstrap
# returns bootstrap replicates and their standard error
# data is a vector, matrix or data-frame
# FUN must act on a vector in 1-dim case, or a matrix or data-frame in N-dim case
bootstrap <- function(data, FUN, ..., B=1000) {
  T_boot <- numeric(B)

  # when data is a vector
  if (is.null(nrow(data))) {
    n <- length(data)

    for (b in 1:B) {
      X_star <- sample(data, n, replace=TRUE)
      T_boot[b] <- FUN(X_star, ...)
    }
  }

  # when data is a matrix or dataframe.
  else {
    n <- nrow(data)

    for (b in 1:B) {
      index <- sample(1:n, n, replace=TRUE)
      X_star <- data[index,]
      T_boot[b] <- FUN(X_star, ...)
    }
  }

  se <- sd(T_boot)
  return(list(se=se, T_boot=T_boot))
}

# correlation function on a data-frame or matrix

#' Pearson's correlation coefficient of bivariate data
#'
#' @param bivar_data A data-frame or matrix, the first two columns of which contain the two vectors for which the correlation coefficient is required
#' @return the Pearson correlation coefficient of bivar_data[,1] and bivar_data[,2]
#' @examples
#' A <- data.frame(x=rnorm(10), y=rnorm(10, mean=2))
#' cor_df(A)
cor_df <- function(bivar_data) {
  cor(bivar_data[,1], bivar_data[,2])
}

# test function on law school data
law_school <- data.frame(LSAT=c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594), GPA=c(3.39,3.30,2.81,3.03,3.44,3.07,3.00,3.43,3.36,3.13,3.12,2.74,2.76,2.88,2.96))
cor(law_school$LSAT, law_school$GPA)
bootstrap_law <- bootstrap(law_school, cor_df, B=1000)
bootstrap_law$se
bootstrap_law_df <- data.frame(T_boot = bootstrap_law$T_boot - cor_df(law_school))
qplot(T_boot, data=bootstrap_law_df, geom="histogram", binwidth=0.02)

# compare our function with boot
set.seed(1)
cor_boot <- function(data, i) cor(data[i,1], data[i,2])
cb <- boot(law_school, cor_boot, 1000)
cb

### Bayesian Bootstrap

# takes a dataset (either vector, matrix or dataframe),
# a function which takes data and weights and returns the statistic of interest,
# and the number of bootstrap samples.
# returns a dataframe with the statistic of interest for each bootstrap sample
# and the standard errors of those statistics
BB <- function(data, FUN, ..., B=1000) {

  # find number of original samples
  # if data is a vector
  if (is.null(nrow(data))) n <- length(data)

  # if data is a matrix or dataframe
  else n <- nrow(data)

  # find dim of FUN output and set up storage space (which will allow plotting)
  init_weights <- rep(1/n, n)
  dims <- dim(FUN(data, init_weights, ...))
  T_boot <- matrix(0, B, prod(dims))

  # generate bootstrap replicates
  for (b in 1:B) {
    # generate weights from uniform Dirichlet
    u <- runif(n-1)
    g <- diff(c(0, sort(u), 1))

    # calculate the statistic, using above weights
    T_boot[b,] <- FUN(data, g, ...)
  }

  # T_boot_df allows for plotting of the bootstrap replicates using ggplot2
  T_boot_df <- data.frame(T_boot=T_boot)
  return(list(se=sd(T_boot), replicates=T_boot_df))
}

# test the BB function with mean
X <- rnorm(100, 4, 2)
BB_mean <- BB(X, weighted.mean, B=10000)
BB_mean$se
qplot(T_boot, data=BB_mean$replicates, geom="histogram", binwidth=0.05)

## test the BB function with correlation and law school data
# takes bivarate data and weights (g) and returns a weighted correlation
weighted.cor <- function(data, g) {
  y <- data[,1]
  z <- data[,2]
  wc_num <- sum(g*y*z)-(sum(g*y))*(sum(g*z))
  wc_den <- (sum(g*y^2)-(sum(g*y))^2)*(sum(g*z^2)-(sum(g*z))^2)
  wc_num/sqrt(wc_den)
}

# test BB using weighted.cor and law school data
set.seed(1)
BB_law <- BB(law_school, weighted.cor, B=1000)
BB_law$se
qplot(T_boot, data=BB_law$replicates, geom="histogram", binwidth=0.02)

### Bag of Little Bootstraps
# data is a vector
# gamma is a number ideally between [0.5, 1] which controls size of subsample
# s is number of nodes
# r is number of replicates per node
# returns the mean of the standard errors of the bootstrap replicates.
BLB.1d <- function(data, FUN, ..., gamma, s=20, r=100) {
  n <- length(data)
  b <- round(n^gamma)
  xis <- numeric(s)

  subsamp_mat <- matrix(0, s, b)
  subsamp <- numeric(b)
  resample_mat <- matrix(0, r, n)

  for (i in 1:s) {
    subsamp <- sample(data, b)
    resample_mat <- matrix(sample(subsamp, r*n, replace = TRUE), r, n)
    theta <- apply(resample_mat, 1, FUN, ...)
    xis[i] <- sd(theta)
  }

  return(mean(xis))
}

# Test BLB.1d
X <- rnorm(5000)
BLB.1d(X, mean, gamma=0.5)
