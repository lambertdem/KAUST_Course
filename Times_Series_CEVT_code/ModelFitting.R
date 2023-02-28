# NEW 2022
# Script for fitting correlation function models to lagged alpha parameters; 
# We consider correlation functions corresponding to stationary Markov processes of 
# orders 1, 2 & 3. 
#===============================================================================
rm(list=ls())
setwd("C:/Users/lambe/Documents/Edinburgh/PhD/CodeSharing_IP_LDM/Time_Series_Conditional_Extremes_code/") 
library(tidyverse)
source("MarkovFitsGaussZ.R")
source("ForwardSimulationFunctions.R")
source("GPDfunctions.R")



# -------------------------------------------------------------------
# Data loading and manipulation
# -------------------------------------------------------------------

load("orleans_data.RData")

# Convert orleans_data$TX to Laplace margins with semi-parametric model
u <- quantile(orleans_data$TX,0.90) # u above which rank transform becomes parametric GPD
orleans_data_lap <- ToLaplace_df(data = orleans_data, thresh = u)


# Sanity check for the Back transformation to data scale
back <- BackTransform(x.lap = orleans_data_lap$data$TX, 
                      y.star = u, 
                      rate = orleans_data_lap$rate, 
                      sig = orleans_data_lap$sig, 
                      xi = orleans_data_lap$xi, 
                      x.orig = orleans_data$TX)
plot(back,orleans_data$TX)

# Thresholds of exceedances above which to fit the models
thresh.lap <- quantile(orleans_data_lap$data$TX, probs=seq(0.85,0.95, length=11))

# -------------------------------------------------------------------
# Model fitting
# -------------------------------------------------------------------

# Number of lags ahead to fit
L <- 20

# fit order 1 model at a range of threshold (for stability plots)
markov1_origdata_all_thresholds <- lapply(thresh.lap, function(i) fit.Markov.wrapper(pars=c(1,0), 
                                                                                     u=i, 
                                                                                     raw.data.lap = orleans_data_lap$data,
                                                                                     k=L,
                                                                                     orderMarkov = "one"))

# save(markov1_origdata_all_thresholds, file="markov1_origdata_all_thresholds.RData")



# now fit to original data set on Laplace scale
markov2_origdata_all_thresholds <- lapply(thresh.lap, function(i) fit.Markov.wrapper(pars=c(1,1,0), 
                                                                      u=i, 
                                                                      raw.data.lap = orleans_data_lap$data, 
                                                                      k=L,
                                                                      orderMarkov = "two"))

# save(markov2_origdata_all_thresholds, file="markov2_origdata_all_thresholds.RData")

# now order 3
markov3_origdata_all_thresholds <- lapply(thresh.lap, function(i) fit.Markov.wrapper(pars=c(1.5,0.75,0.35,0), 
                                                                      u=i,
                                                                      raw.data.lap = orleans_data_lap$data,
                                                                      k=L, 
                                                                      orderMarkov = "three"))

# save(markov3_origdata_all_thresholds, file="markov3_origdata_all_thresholds.RData")

# The above needs to be repeated for each bootstrap sample to produce parameter statbility plots...

# -------------------------------------------------------------------
# Performing forward simulation
# -------------------------------------------------------------------

# Select model to forward-simulate from
fitted <- markov1_origdata_all_thresholds

# Filter data from the Orleans data that exceed the 0.95 quantile on Laplace margins
filt.data <- filterData(orleans_data_lap$data)

# Get cluster lengths
lengthsList <- lapply(filt.data,length)

# Get the fitted alphas and betas using exceedances of 0.95 quantile on Laplace margins
alphas <- fitted$`95%`$alphas
beta <- fitted$`95%`$beta

# Get a list of residuals from the CEVT model fit
z.list <- get.z(filt.data, lengthsList, alphas, beta, L)

# Number of forward simulations to run given X.0 > u
n.sim <- 1000

# Select a threshold above which to simulate X.0 with standard exp.
quantile(orleans_data$TX,0.96)
X.0 <- quantile(orleans_data_lap$data$TX,0.96)

# Run the forward simulation based on X.0, fitted alphas and beta, residuals
clusters <- forwardSim(X.0, alphas, beta, z.list, n.sim)

# Store the forward simulation in observations scale
frwrd.sim <- matrix(NA,nrow=n.sim,ncol=L+1)
for(i in 1:n.sim){
  frwrd.sim[i,] <- BackTransform(x.lap = clusters[[i]], 
                                 y.star = u, 
                                 rate = orleans_data_lap$rate, 
                                 sig = orleans_data_lap$sig, 
                                 xi = orleans_data_lap$xi, 
                                 x.orig = orleans_data$TX)
}

# Plot forward simulations on Laplace scale
par(mfrow=c(1,1),mar=c(4.1,5.1,0.5,1.1),mgp=c(2.6,0.8,0))
plot(clusters[[1]],type="l",ylim=c(min(unlist(clusters)),max(unlist(clusters))))
for(i in 2:n.sim){
  lines(clusters[[i]])
}

# Plot forward simulations on observation scale
par(mfrow=c(1,1),mar=c(4.1,5.1,0.5,1.1),mgp=c(2.6,0.8,0))
plot(frwrd.sim[i,],type="l",ylim=c(min(frwrd.sim),max(frwrd.sim)))
for(i in 2:n.sim){
  lines(frwrd.sim[i,])
}

frwrd.sim

# -------------------------------------------------------------------
# Estimation of p(v, d, r)
# -------------------------------------------------------------------

# Function to compute p(u, d, j)
pmf.clust.length.scalar <- function(frwrd.sim,u,d,j){
  res <-apply(frwrd.sim[,1:d],1,function(x){
              as.numeric(sum(x>u)==j)}
              )
  mean(res)
}

# Vectorize function pmf.clust.length.scalar over d,j, and u
pmf.clust.length <- Vectorize(pmf.clust.length.scalar,vectorize.args = c("d","j","u"))

ds <- c(5,10,20)
us <- c(32.1,35)
js <- c(1:L)
g <- as.data.frame(expand.grid(j=js,d=ds,u=us)) %>% dplyr::filter(j<=d)

# Compute p(u, d, j) for each combination in grid g
g$prob <- pmf.clust.length(frwrd.sim,u=g$u,d=g$d,j=g$j)

# Plot p(u, d, j=1:d) for different u's and d's
par(mfrow=c(length(ds),length(us)),mar=c(4.1,5.1,2,1.1),mgp=c(2.6,0.8,0))
for(r in 1:length(ds)){
  for(c in 1:length(us)){
    g.rc <- (g%>% filter(u==us[c]&d==ds[r]))$prob
    if(r==1&c==1){
      plot(g.rc,main=paste(rep("V",2*c),collapse=""),
           ylab=paste(rep("Y",2*r),collapse=""),xlab="",
           cex.lab=2,cex.axis=1.6,cex.main=2)
    }
    else if(r==length(ds)&c==1){
      plot(g.rc,ylab=paste(rep("Y",2*r),collapse=""),xlab="XXX",
           cex.lab=2,cex.axis=1.6,cex.main=2)
    }
    else if(r==1){
      plot(g.rc,main=paste(rep("V",2*c),collapse=""),xlab="",ylab="",
           cex.lab=2,cex.axis=1.6,cex.main=2)
    }
    else if(c==1){
      plot(g.rc,ylab=paste(rep("Y",2*r),collapse=""),
           xlab="",
           cex.lab=2,cex.axis=1.6,cex.main=2)
    }
    
    else if(r==length(ds)){
      plot(g.rc,xlab="XXX",ylab="",cex.lab=2,cex.axis=1.6,cex.main=2)
    }
    else{
      plot(g.rc,xlab="",ylab="",cex.lab=2,cex.axis=1.6,cex.main=2)
    }
  }
}
