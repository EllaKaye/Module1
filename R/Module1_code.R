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

# takes bivarate data and weights (g) and returns a weighted correlation
weighted.cor <- function(data, g) {
  y <- data[,1]
  z <- data[,2]
  wc_num <- sum(g*y*z)-(sum(g*y))*(sum(g*z))
  wc_den <- (sum(g*y^2)-(sum(g*y))^2)*(sum(g*z^2)-(sum(g*z))^2)
  wc_num/sqrt(wc_den)
}

### Bag of Little Bootstraps
# data is a vector
# gamma is a number ideally between [0.5, 1] which controls size of subsample
# s is number of nodes
# r is number of replicates per node
# returns the mean of the standard errors of the bootstrap replicates.
BLB.1d <- function(data, gamma, FUN, ..., s=20, r=100) {
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

#' BLB function that computes the (1-alpha)\% CI and the standard errors of the parameter estimate components
#' @title BLB for multivariate data
#' @param data a matrix or dataframe
#' @param s number of subsamples
#' @param r number of resamples per subsample
#' @param gamma specifies the size of the subsample by b=n^gamma
#' @param lambda is the penalty in L2 norm for the ridge regression
#' @param alpha the level of the confidence for CI, respecified to 0.05
#' @return a list consisting of
#' \item{output_CI}{average confidence intervals for the parameter estimate components}
#' \item{output_se}{average standard errors for the parameter estimates components}
#' @examples
#' n <- 100
#' d <- 5
#' X <- mvrnorm(n ,mu = rep(0 ,d ) ,Sigma <- diag( d ) )
#' epsilon <- rnorm(n ,mean = 0 ,sd = sqrt(10 ) )
#' t_theta <- as.matrix( rep(1 ,d ) )
#' Y <- X%*%t_theta + epsilon
#' data <- cbind(X ,Y )
#' my_result <- BLB.multi(data, 0,7, 15 ,100 ,10^(-5) )
BLB.multi <- function(data, gamma=0.7, s=20, r=100, lambda=10^-5, alpha=0.05) {
  #CI = confidence interval
  #se = standard error

  n <- nrow(data)
  d <- ncol(data)-1
  b <- round(n^gamma)
  samp_ind <- 0
  subsamp <- matrix(0 ,b ,d+1 )
  resamp_ind <- 0
  resamp <- array(0 ,dim = c(n ,d+1 ,r ))
  re_theta <- matrix(0 ,d ,r )
  theta <- matrix(0 ,d ,s )
  sd_theta <- matrix(0 ,d ,s )

  #to store the values for CI's we create an array
  # and for se a matrix
    xis <- array(0 ,dim=c(d ,2 ,s ) )
    output_CI <- matrix(0 ,d ,2 )
    output_se <- matrix(0 ,d ,1 )

  # subsample from the original data
  for ( i in 1:s ) {
    samp_ind <- sample(1:n ,b ,replace=FALSE )
    subsamp <- data[samp_ind , ]

    #resample from the subsample and compute the parameter estimate
    for( j in 1:r ) {
      resamp_ind <- sample(1:b ,n ,replace=TRUE )
      resamp[,,j] <- subsamp[resamp_ind , ]
      Y <- resamp[ ,d+1 ,j ]
      X <- resamp[ ,1:d ,j ]
      re_theta[ ,j ] <- as.matrix(lm.ridge(Y~X,lambda = lambda )$coef )
    }

    #compute the se of the parameter estimates for each subsample which is used in both CI and se
    sd_theta[ ,i ] <- apply(re_theta ,1 ,sd )

    #if the task is CI compute it for each subsample
    theta[,i] <- rowMeans(re_theta)
    xis[,,i] <- cbind(theta[ ,i ],theta[ ,i ]) + cbind(qnorm(alpha/2) * sd_theta[ ,i ] ,qnorm(1-alpha/2) * sd_theta[ ,i ])
  }

  # average CI across subsamples
  output_CI <- apply(xis ,c(1 ,2 ) ,mean )
  CI_width <- output_CI[,2]-output_CI[,1]

  # average se across subsample
  output_se <- apply(sd_theta ,1 ,mean )

  return(list(CI=output_CI, CI_width=CI_width, se=output_se))
}

#' BLB function which returns the average (across dimensions) confidence interal width or the average standard error of the parameter and the values of
#' r and s which were chosen adaptively within the function
#' @title BLB with adaptive selection on r and s
#' @param data a matrix or dataframe
#' @param gamma specifies the subsample size by b=n^gamma
#' @param w_s window size for the adaptive selection of s
#' @param w_r window size for the adaptive selection of r
#' @param lambda specifies the L2 penalty in Ridge regression
#' @param err relative error for the convergence criterion in adaptive selection
#' @param alpha the level of the confidence which is set to 0.05
#' @return a list consisting of
#' \item{s}{the total number of subsamples from the original dataset}
#' \item{r}{the number of resamples for each subsample}
#' \item{mean_width}{mean_width the average marginal width of the confidence intervals across dimensions}
#' \item{mean_se}{mean_se the average standard error across the parameter dimensions}
#' @examples X <- mvrnorm(100 ,mu = rep(0,2) ,Sigma <- diag(2) )
#' epsilon <- rnorm(100,mean=0,sd=sqrt(10))
#' t_theta <- as.matrix( rep(1 ,2 ) )
#' Y <- X%*%t_theta + epsilon
#' my_data <- cbind(X ,Y )
#' my_result <- BLB.adapt(my_data, gamma=0.7, w_s=3, w_r=20, lambda=10^(-5), epsilon=0.05, alpha=0.05)
BLB.adapt <- function(data, gamma=0.7, w_s=3, w_r=20, lambda=10^(-5), epsilon=0.05, alpha=0.05) {
  # initialise storage
  if (!is.matrix(data)) data <- as.matrix(data)
  n <- nrow( data )
  d <- ncol( data ) - 1
  r_max <- 500  #maximum number of resamples of each subsampe (to be found)
  s_max <- 50 #maximum number of subsamples (to be found)
  b <- round( n^gamma )
  samp_ind <- 0
  subsamp <- matrix(0 ,b ,d+1 )
  resamp_ind <- 0 # indictor
  resamp <- array(0 ,dim = c(n ,d+1 ,r_max ) ) # all the resamples from subsample
  re_theta <- matrix(0,d,r_max) # theta on resamples
  re_theta_sd <- matrix(0 ,d ,r_max ) # matrix that stores the sd of the thetas in each resample
  differ_r <- numeric( d )
  differ_s <- numeric( d)
  theta <- matrix(0 ,d ,s_max )
  sd_theta <- width <- final_sd_theta <- matrix(0 ,d ,s_max ) # final sd of each component
  xis <- array(0 ,dim = c(d ,2 ,s_max ) ) # final confidence intervals
  diff_vec_r <- numeric( w_r ) # vector for checking convergence
  diff_vec_s <- numeric(w_s)
  out_s <- numeric( w_s ) # vectors used in the convergence testing
  out_r <- numeric( w_r ) # vectors used in the convergence testing
  vec_of_r <- numeric( s_max ) # vector to store number of resamples for each subsamples
  final_standard_dev <- numeric( d )
  mean_width<-mean_se<-0
  STATUS_S <- FALSE
  STATUS_R <- FALSE
  s<-0 # initialise number of subsamples

  #adaptive selection of s
  while( STATUS_S == FALSE ) {
    s <- s+1

    #subsample from the original dataset
    samp_ind <- sample(1:n ,b ,replace = FALSE )
    subsamp <- data[samp_ind ,]
    r <- 0
    STATUS_R <- FALSE

    # adaptive selection of r
    # for that subsample resample data of size n until the parameter estimate converges
    while( STATUS_R == FALSE ) {
      r <- r+1
      resamp_ind <- sample(1:b ,n ,replace = TRUE )
      resamp[,,r] <- subsamp[resamp_ind,]
      Y <- resamp[ ,d+1 ,r ]
      X <- resamp[ ,1:d ,r ]
      re_theta[ ,r ] <- as.matrix(lm.ridge( Y~X ,lambda = lambda)$coef )

      # store values of z which are the se's of the components of the parameter estimate
      #if (r==1) {
      #re_theta_sd[,1]<-c(rep(NA,d))
      #}

      if (r > 1)  {re_theta_sd[ ,(r-1) ] <- apply(re_theta[,1:r],1,sd)}


      #beyond the window value test the convegence condition
      if (r > w_r) {
        out_r<-numeric(w_r)

        #for the past w values we test the value of the relative error
        for ( k in 1:w_r ) {
          differ_r <- re_theta_sd[ ,(r-1) ] - re_theta_sd[ ,( r-1-k )]
          diff_vec_r[ k ] <- sum( abs( differ_r ) / abs( re_theta_sd[ ,(r-1) ] ) ) / d

          # to chec whether the condition for all the w values is satisfied
          #store 1 if yes and 0 if not
          if ( (diff_vec_r[ k ] < epsilon )| (diff_vec_r[ k ] == epsilon ))
            out_r[k]<-1
        }

        #terminate the loop if all the values of the relative error are less than epsilon

        if ( sum( out_r ) == w_r )

          STATUS_R <-TRUE
      }

    }


    # vector that stores the values of r selected adaptively for each subsample
    vec_of_r[ s ] <- r

    # store the parameter estimates for each subsample
    theta[ ,s ] <- rowMeans( re_theta[ ,1:r ] )
    sd_theta[ ,s ] <- apply( re_theta[ ,1:r ] ,1 ,sd )
    xis[ , ,s ] <- cbind(theta[ ,s ] ,theta[ ,s ] ) + cbind( -qnorm(1-alpha/2) * sd_theta[ ,s ], qnorm(1-alpha/2) * sd_theta[ ,s ])
    #width[ ,s ] <- 2 * 1.959964 * sd_theta[ ,s ]

    # store values of z which are the se's of the components of parameter estimate for that subsample
    if ( s == 1 ) {
      final_sd_theta[,1]<-sd_theta[,1]
    }

    if ( s > 1 ) {
      final_sd_theta[,s]<-apply(sd_theta[,1:s],1,mean)
    }

    #beyond the window value check the convergence condition
    if (s > w_s) {
      out <- numeric( w_s )
      for ( l in 1:w_s )  {
        differ_s <- final_sd_theta[ ,s ] - final_sd_theta[ ,(s-l) ]
        diff_vec_s[ l ] <- sum(abs(differ_s)/abs(final_sd_theta[ ,s ]))/d
        if ( (diff_vec_s[ l ] < epsilon ) | (diff_vec_s[ l ] == epsilon ))
          out_s[ l ] <- 1
      }

      #terminate the loop if all the values of the relative error are less than epsilon
      if ( sum( out_s ) == w_s)
        STATUS_S <-TRUE

    }
  }

  #store the final (average) se of the parameter estimates
  final_standard_dev <- final_sd_theta[ ,s ]
  #compute the CI for the parameter estimates across dimensions
  final_xi <- apply(xis[ , ,1:s ] ,c(1,2) ,mean )
  #compute the width of the C's
  widths <- final_xi[ ,2 ] - final_xi[ ,1 ]
  #return the final number of subsamples, resamples or each subsample and average CI width

  return( list(s = s, r = vec_of_r[ 1:s ] , mean_width = mean(widths), mean_se=mean(final_standard_dev) ) )

}


