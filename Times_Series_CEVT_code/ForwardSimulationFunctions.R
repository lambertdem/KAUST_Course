# NEW 2022
# function for forward simulation of clusters simulating forward from a threshold exceedance)

# first we write a function to obtain z (residual vector) given alpha & beta parameters 
# given data in the form as returned by the function filterData from MarkovGaussZ.R script
#================
# ARGS:
# filt.data: a list, filtered data of the form produced by the function filterData from MarkovGaussZ.R script
# lengthsList: a numeric vector which gives the length of each elements of lengthsList
# alphas: vector of alpha values of length k
# beta: numeric of length 1; the fitted value of beta
# Value:
# a list of residual z vectors each of length k
get.z <- function(filt.data, lengthsList, alphas, beta, k){
  if(length(alphas)!=k){
    stop("alphas vector is not of length k")
  }else{
  get.ind <- which(lengthsList >= k+1)    # get clusters in which there is a at least k subsequent values after exceedance i.e. vectors must be of length k+1
  z.list  <- lapply(get.ind, function(i) (filt.data[[i]][2:(k+1)] - alphas[1:k]*filt.data[[i]][1]) / filt.data[[i]][1]^beta)
  }
  z.list
}
#===============================================================================
#===============================================================================
#===============================================================================
# function for simulation forward from a threshold exceedance (on Laplace scale)
#================
# ARGS:
# u.thresh: single numeric, the threshold above which we simulate an exceedance
# alphas: vector of lagged alpha parameters; length is the number of values after exceedance we wish to simulate
# n.sim = number of simulated clusters

forwardSim <- function(u.thresh, alphas, beta, z.list, n.sim){
  single.sim <- function(){     # function for generating a single cluster
  x.init <- u.thresh + rexp(1)  # initial exceedance 
  n <- length(z.list)
  which.z <- sample(1:n, size=1)
  forward.vals <- alphas*x.init + (x.init^beta)*z.list[[which.z]] 
  c(x.init, forward.vals) # simulated cluster
  }
  lapply(1:n.sim, function(i) single.sim())
}
#===============================================================================
#===============================================================================
#===============================================================================

