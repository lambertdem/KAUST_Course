# NEW 2022
# GPD fitting functions & functions for transforming to Laplace scale based on semi-parametric cdf fit
#===============================================================================
#===============================================================================
# negative GPD log-likelihood
# ARGS:
# pars = (0.1*xi, sigma(1 + xi)); we use the orthogonal parameterization; factor of 0.1 is for better numerical stability
# x: data vector 
# threshold: scalar threshold above which 
#============================
negllGPD <- function(pars, x, threshold, xi.tol=0.01){
  x <- x[x > threshold] - threshold  # threshold excesses 
  xi  <- 0.1*pars[1] 
  sig <- pars[2] / (1 + xi)    
  n <- length(x)
  if(sig < 0){
    ll <- -1e10
  } else {
  if(abs(xi) <= xi.tol){ # use Gumbel likelihood when xi is close to zero
    ll <- -n*log(sig) - (1/sig)*sum(x)
  } else {
    ll <- max(-1e10, -n*log(sig) - (1/xi + 1)*sum(log(pmax(1+xi*x/sig, 0)))) # need max as second term may eval to -Inf
   }
  }
  -ll
}
#=============================
# function to find GPD MLEs; xi & sigma in classic parameterization (and rate estimate)  
my.fit.gpd <- function(pars, x, threshold, xi.tol=0.01,...){
  fit <- optim(par=pars, fn=negllGPD, x=x, threshold=threshold, xi.tol=xi.tol,...)
  fit.mle <- fit$par
  xi  <- 0.1*fit.mle[1]
  sig <- fit.mle[2] / (1 + xi)
  mles <- c(xi, sig)
  names(mles) <- c("shape", "scale")
  rate <- sum(x > threshold) / length(x) # proportion of data that exceed threshold
  list(mles=mles, rate=rate)
}
#=======================================================================================
#=======================================================================================
#=======================================================================================
#=======================================================================================
#=======================================================================================
# function for transforming data to Laplace scale 

# ARGS:
# x: numeric vector that we wish to transform to Laplace margins
# thresh: threshold used, above which GPD model is used for estimation

ToLaplace <- function(x, threshold, pars=c(-1,2),...){
  # first step is to get PIT values for each data point using semi-parametric model; 
  # store in vector called cdf
  ranks <- rank(x) # get ranks
  ecdf  <- ranks / (length(x) + 0.5) # empirical cdf
  exc   <- which(x > threshold) 
  non.exc <- which(x <= threshold)
  fit1 <- my.fit.gpd(pars=pars, x=x, threshold=threshold,...) 
  sig <- fit1[["mles"]]["scale"]
  xi  <- fit1[["mles"]]["shape"]
  exec.rate <- fit1[["rate"]] 
  cdf <- NULL
  cdf[non.exc] <- ecdf[non.exc]
  cdf[exc] <- 1 - exec.rate*(pmax((1 + xi*(x[exc] - threshold)/sig), 0))^(-1/xi)
  lower.half <- which(cdf < 0.5)
  upper.half <- which(cdf >= 0.5)
  lap <- NULL # now transform to Laplace scale
  lap[lower.half] <- log(2*cdf[lower.half])
  lap[upper.half] <- -log(2*(1 - cdf[upper.half]))
  # return transformed x on Laplace scale and GPD parameter estimates & exceedance rate
  return(list(lap=lap, sig=sig, xi=xi,rate=exec.rate)) 
}
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
# now write a function that takes as input a dataframe with (at least) columns TX and year and returns 
# TX on Laplace scale
# here data is a dataframe as opposed to previous function where it is a numeric vector

ToLaplace_df <- function(data, thresh){
  TX.lap <- ToLaplace(x=data$TX, thresh=thresh) # returns a list
  data.lap <- TX.lap[[1]] # get just Laplace data
  data$TX  <- data.lap
  list(data=data, sig=TX.lap[["sig"]], xi=TX.lap[["xi"]], rate=TX.lap[["rate"]])
}
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
# function to back transform to original scale from Laplace scale
#==========================
# ARGS: 
# x.lap: numeric vector on Laplace scale
# y.star: threshold that was used for fitting GPD model (= 29.7 in analysis)
# rate: estimated rate/probability of exceedaing y.star 
# sig: estimated GPD scale parameter
# xi:  estimated GPD shape parameter
# x.orig: numeric vector of observations on original scale (will be TX column of a dataframe in our application)
BackTransform <- function(x.lap, y.star, rate, sig, xi, x.orig){
  single.transform <- function(x){ # first function to transform a single numeric
    if(x <= 0){
      p <- exp(x) / 2
    } else {
      p <- 1 - exp(-x) / 2
    }
    if(p <= (1 - rate)){
      q <- unname(quantile(x.orig, probs=p)) 
    } else {
      q <- y.star + (sig/xi)*(-1 + ((1 - p)/rate)^(-xi))
    }
    q
  } # end of subfunction
  # now apply subfunction to all elements of x.lap
  sapply(1:length(x.lap), function(i) single.transform(x=x.lap[i]))
}
#===============================================================================
#===============================================================================
