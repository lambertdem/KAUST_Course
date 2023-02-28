# Functions for fitting TSCE model with correlation function structure in alpha parameters
# and the semi-parametric approach (using working assumption of Gaussian residuals)
library(tidyverse)
#===============================================================================
#===============================================================================
#===============================================================================
# function to take orleans data (on Laplace scale) and return a list of numeric vectors;
# first element of each vector corresponds to a threshold exceedance and is followed by k subsequent values
# data should be in format of orleans_data_lap.RData and have columns named TX and year
# (norming argument can be removed; only classic is implemented)

filterData <- function(data, k=20, u=-log(2*(1-0.95))){
  N <- nrow(data %>% filter(TX > u)) # get number of exceedances
  out.data <- vector("list", N)
  M <- 0  # counting variable
  for(i in 1946:2012){
    year.i.data <- data  %>% filter(year==i)
    year.i.TX   <- year.i.data$TX
    exc <- which(year.i.TX > u) # threshold exceedance times in year i 
    for(j in 1:length(exc)){
      if(length(exc)!=0){
      M <- M + 1
      out.data[[M]] <-  c(year.i.TX[exc[j]], 
                          year.i.TX[(min(exc[j]+1, nrow(year.i.data))):(min(exc[j]+k, nrow(year.i.data)))])
      } 
    }  
  }
  out.data
}
#===============================================================================
#===============================================================================
#===============================================================================
# ARGS:
# pars: numeric vector of length 1 (transformed parameters alpha and beta)
# data: a list of numeric vectors; first element of each vector is an exceedance; other elements are subsequent values 
# (data will be of format produces by filterData function above)
Markov1 <- function(pars, 
                    data,
                    k = 20,
                    lengthsList,
                    norming=c("classic", "alt")){
  norming <- match.arg(norming)
  alpha <- (2/pi)*atan(pars[1])
  beta <- 1 / (1 + exp(-pars[2]))
  alphas <- alpha^(1:k) # alpha parameters for lags 1 to k 
  # now write function to get profile likelihood estimates of (mu_i, sigma_i), i=1...k
  prof.Mu.Sig <- function(no.Lags){
    list.indices <- which(lengthsList >= no.Lags+1) # + 1 is because initial exceedance is included in list
    lagged.data <- sapply(list.indices, function(i) data[[i]][no.Lags+1]) 
    data.exc <- sapply(list.indices, function(i) data[[i]][1])
    if(norming=="classic"){
    z <- (lagged.data - alphas[no.Lags]*data.exc) / (data.exc^beta) 
    } else {
    z <- (lagged.data - alphas[no.Lags]*data.exc) / (1 + (alphas[no.Lags]*data.exc)^beta)  
    }
    mu.z <- mean(z); sig.z <- sd(z)
    c(mu.z, sig.z)
  } 
  # now sapply prof.Mu.Sig function to get all estimates of mu_i & sigma_i
  mus.sigs.z <- lapply(1:k, function(i) prof.Mu.Sig(i)) 
  mus.prof   <- sapply(1:k, function(i) mus.sigs.z[[i]][1]) 
  sigs.prof  <- sapply(1:k, function(i) mus.sigs.z[[i]][2]) 
  # now get log-lik contributions for each lag from 1 to k
  ll_cont <- rep(0, k) # to store likelihood contributions from different lags
  for(i in 1:k){
    list.indices <- which(lengthsList >= i+1) 
    y <- sapply(list.indices, function(j) data[[j]][1])
    y.lag <- sapply(list.indices, function(j) data[[j]][i+1]) 
    if(norming=="classic"){
    mu.y.lag <- alphas[i]*y + (y^beta)*mus.prof[i]
    sd.y.lag <- (y^beta)*sigs.prof[i]
    }else{
    mu.y.lag <- alphas[i]*y + (1 + (alphas[i]*y)^beta)*mus.prof[i]
    sd.y.lag <- (1 + (alphas[i]*y)^beta)*sigs.prof[i]
    }
    ll_cont[i] <- sum(dnorm(x=y.lag, mean=mu.y.lag, sd=sd.y.lag, log=TRUE))
  }
  -sum(ll_cont)
}
#===============================================================================
#===============================================================================
#===============================================================================
# ARGS:
# pars: numeric vector of length 3 (transformed parameters pi1, pi2 and beta)
# data: a list of numeric vectors; first element of each vector is an exceedance; other elements are subsequent values 
# (norming argument can be removed; only classic is implemented)
Markov2 <- function(pars, 
                    data,
                    k = 20,
                    lengthsList,
                    norming=c("classic", "alt")){
  norming <- match.arg(norming)
  # pi1 and pi2 correspond to partial autocorrelations in autocorrelation function parametrization
  pi1 <- (2/pi)*atan(pars[1])
  pi2 <- (2/pi)*atan(pars[2])
  beta <- 1 / (1 + exp(-pars[3]))
  # theta1 and theta2 correspond to AR coefficients
  theta1 <- pi1*(1 - pi2)
  theta2 <- pi2 
  alphas <- rep(0, k)
  alphas[1] <- theta1 / (1 - theta2)
  alphas[2] <- theta1*alphas[1] + theta2*1  # alpha0 = 1 
  for(i in 3:k){
    alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2]
  }
  # now write function to get profile likelihood estimates of (mu_i, sigma_i), i=1...k
  prof.Mu.Sig <- function(no.Lags){
    list.indices <- which(lengthsList >= no.Lags+1) # + 1 is because initial exceedance is included in list
    lagged.data <- sapply(list.indices, function(i) data[[i]][no.Lags+1]) 
    data.exc <- sapply(list.indices, function(i) data[[i]][1])
    if(norming=="classic"){
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (data.exc^beta) 
    } else {
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (1 + (alphas[no.Lags]*data.exc)^beta)  
    }
    mu.z <- mean(z); sig.z <- sd(z)
    c(mu.z, sig.z)
  } 
  # now sapply prof.Mu.Sig function to get all estimates of mu_i & sigma_i
  mus.sigs.z <- lapply(1:k, function(i) prof.Mu.Sig(i)) 
  mus.prof   <- sapply(1:k, function(i) mus.sigs.z[[i]][1]) 
  sigs.prof  <- sapply(1:k, function(i) mus.sigs.z[[i]][2]) 
  # now get log-lik contributions for each lag from 1 to k
  ll_cont <- rep(0, k) # to store likelihood contributions from different lags
  for(i in 1:k){
    list.indices <- which(lengthsList >= i+1) 
    y <- sapply(list.indices, function(j) data[[j]][1])
    y.lag <- sapply(list.indices, function(j) data[[j]][i+1]) 
    if(norming=="classic"){
      mu.y.lag <- alphas[i]*y + (y^beta)*mus.prof[i]
      sd.y.lag <- (y^beta)*sigs.prof[i]
    }else{
      mu.y.lag <- alphas[i]*y + (1 + (alphas[i]*y)^beta)*mus.prof[i]
      sd.y.lag <- (1 + (alphas[i]*y)^beta)*sigs.prof[i]
    }
    ll_cont[i] <- sum(dnorm(x=y.lag, mean=mu.y.lag, sd=sd.y.lag, log=TRUE))
  }
  -sum(ll_cont)
}
#===============================================================================
#===============================================================================
#===============================================================================
# now function for Markov order 3 model
#===============================================================================
# ARGS:
# pars: numeric vector of length 4 (transformed parameters pi1, pi2,pi3 and beta)
# data: a list of numeric vectors; first element of each vector is an exceedance; other elements are subsequent values 
# (norming argument can be removed; only classic is implemented)
Markov3 <- function(pars, 
                    data,
                    k = 20,
                    lengthsList,
                    norming=c("classic", "alt")){
  norming <- match.arg(norming)
  # pi1...pi3 correspond to partial autocorrelations
  pi1 <- (2/pi)*atan(pars[1])
  pi2 <- (2/pi)*atan(pars[2])
  pi3 <- (2/pi)*atan(pars[3])
  beta <-  1 / (1 + exp(-pars[4])) 
  # theta1, theta2 & theta3 correspond to AR coefficients
  theta1 <- pi1 - pi1*pi2 - pi2*pi3 
  theta2 <- pi2 - pi1*pi3 + pi1*pi2*pi3 
  theta3 <- pi3 
  alphas <- rep(0, k) # to store sequence of lagged alpha values
  alphas[1] <- (theta1 + theta2*theta3) / (1 - theta2 - theta1*theta3 - theta3^2) 
  alphas[2] <- theta2 + (theta1 + theta3)*alphas[1]  
  alphas[3] <- theta1*alphas[2] + theta2*alphas[1] + theta3*1  # alpha0 = 1
  for(i in 4:k){
    alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2] + theta3*alphas[i-3]
  }
  # now write function to get profile likelihood estimates of (mu_i, sigma_i), i=1...k
  prof.Mu.Sig <- function(no.Lags){
    list.indices <- which(lengthsList >= no.Lags+1) # + 1 is because initial exceedance is included in list
    lagged.data <- sapply(list.indices, function(i) data[[i]][no.Lags+1]) 
    data.exc <- sapply(list.indices, function(i) data[[i]][1])
    if(norming=="classic"){
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (data.exc^beta) 
    } else {
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (1 + (alphas[no.Lags]*data.exc)^beta)  
    }
    mu.z <- mean(z); sig.z <- sd(z)
    c(mu.z, sig.z)
  }  
  # now sapply prof.Mu.Sig function to get all estimates of mu_i & sigma_i
  mus.sigs.z <- lapply(1:k, function(i) prof.Mu.Sig(i)) 
  mus.prof   <- sapply(1:k, function(i) mus.sigs.z[[i]][1]) 
  sigs.prof  <- sapply(1:k, function(i) mus.sigs.z[[i]][2]) 
  # now get log-lik contributions for each lag from 1 to k
  ll_cont <- rep(0, k) # to store likelihood contributions from different lags
  for(i in 1:k){
    list.indices <- which(lengthsList >= i+1) 
    y <- sapply(list.indices, function(j) data[[j]][1])
    y.lag <- sapply(list.indices, function(j) data[[j]][i+1]) 
    if(norming=="classic"){
      mu.y.lag <- alphas[i]*y + (y^beta)*mus.prof[i]
      sd.y.lag <- (y^beta)*sigs.prof[i]
    }else{
      mu.y.lag <- alphas[i]*y + (1 + (alphas[i]*y)^beta)*mus.prof[i]
      sd.y.lag <- (1 + (alphas[i]*y)^beta)*sigs.prof[i]
    }
    ll_cont[i] <- sum(dnorm(x=y.lag, mean=mu.y.lag, sd=sd.y.lag, log=TRUE))
  }
  -sum(ll_cont)
}
#===============================================================================
#===============================================================================
# function for fitting Markov model of specified order and returning fitted parameters

fit.Markov <- function(pars, 
                       data,
                       k = 20,
                       lengthsList,
                       norming=c("classic", "alt"),
                       orderMarkov=c("one","two","three"),...){
  
  orderMarkov <- match.arg(orderMarkov)
  if(orderMarkov=="one"){
    fit <- optim(par=pars,fn=Markov1, data=data, k=k, lengthsList=lengthsList, norming=norming)
    raw.pars <- fit$par
    alpha <- (2/pi)*atan(raw.pars[1])
    beta <- 1 / (1 + exp(-raw.pars[2]))
    alphas <- alpha^(1:k) 
    out.list <- list(alpha=alpha, beta=beta, alphas=alphas)
  } else if(orderMarkov=="two"){
    fit <- optim(par=pars,fn=Markov2, data=data, k=k, lengthsList=lengthsList, norming=norming)
    raw.pars <- fit$par
    pi1 <- (2/pi)*atan(raw.pars[1]) # pi1 & pi2 are partial autocorrelations 
    pi2 <- (2/pi)*atan(raw.pars[2])
    beta <- 1 / (1 + exp(-raw.pars[3]))
    # theta1 and theta2 correspond to AR coefficients
    theta1 <- pi1*(1 - pi2)
    theta2 <- pi2 
    alphas <- rep(0, k)
    alphas[1] <- theta1 / (1 - theta2)
    alphas[2] <- theta1*alphas[1] + theta2*1  # alpha0 = 1 
    for(i in 3:k){
      alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2]
    }
    out.list <- list(pi1=pi1, pi2=pi2, theta1=theta1,theta2=theta2, beta=beta, alphas=alphas)
  } else{
    fit <- optim(par=pars,fn=Markov3, data=data, k=k, lengthsList=lengthsList, norming=norming)
    raw.pars <- fit$par
    # pi1...pi3 correspond to partial autocorrelations
    pi1 <- (2/pi)*atan(raw.pars[1])
    pi2 <- (2/pi)*atan(raw.pars[2])
    pi3 <- (2/pi)*atan(raw.pars[3])
    beta <-  1 / (1 + exp(-raw.pars[4])) 
    # theta1, theta2 & theta3 correspond to AR coefficients
    theta1 <- pi1 - pi1*pi2 - pi2*pi3 
    theta2 <- pi2 - pi1*pi3 + pi1*pi2*pi3 
    theta3 <- pi3 
    alphas <- rep(0, k)
    alphas[1] <- (theta1 + theta2*theta3) / (1 - theta2 - theta1*theta3 - theta3^2) 
    alphas[2] <- theta2 + (theta1 + theta3)*alphas[1]  
    alphas[3] <- theta1*alphas[2] + theta2*alphas[1] + theta3*1  # alpha0 = 1
    for(i in 4:k){
      alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2] + theta3*alphas[i-3]
    }
    out.list <- list(pi1=pi1, pi2=pi2,pi3=pi3, 
                     theta1=theta1,theta2=theta2,theta3=theta3, 
                     beta=beta, alphas=alphas)
  }
  out.list
}
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
# wrapper function to simplify markov correlation fits; 
# takes raw data on Laplace scale, filters it and produces 
# lengths list as required by fit markov function
fit.Markov.wrapper <- function(pars,
                               u=-log(2*(1-0.9)),
                               raw.data.lap,
                               k = 20,
                               norming=c("classic", "alt"),
                               orderMarkov=c("one","two","three"),...){
  filt.data1 <- filterData(data=raw.data.lap, k=20, u=u)
  lengthsList <- sapply(filt.data1, function(i) length(i))
  # now call fit.Markov
  fit.Markov(pars = pars, 
             data=filt.data1,
             lengthsList = lengthsList,
             k=k,
             norming=norming,
             orderMarkov=orderMarkov,
  )
}
#============





