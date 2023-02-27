rm(list=ls())

ChernickAR <- function(k,n){
  p <- c(0:(k-1))/k
  eps <- sample(p,n,replace=T)
  
  X <- rep(NA,n); X[1] <- eps[1]
  for(i in 2:n){
    X[i] <- X[i-1]/k+eps[i]
  }
  X
}

blockbootstrap <- function(X,blocksize){
  n <- length(X)
  boot.samp <- rep(NA,n)
  ind <- c(1:n)
  ind.strt <- sample(ind,n/blocksize,replace=T)
  for(i in 1:(n/blocksize)){
    block.ind <- c(ind.strt[i]:(ind.strt[i]+blocksize-1)) %% n
    block.ind[block.ind==0] <- n
    
    ind.new.samp <- c(((i-1)*blocksize+1): (i*blocksize))
    boot.samp[ind.new.samp] <- X[block.ind]
  }
  boot.samp
}

theta.bootstrap <- function(X,n.samp,us,block.size){
  CIs <- matrix(NA,ncol=6,nrow=length(us))
  colnames(CIs) <- c("lower_r1","upper_r1","lower_r2","upper_r2","lower_fs","upper_fs")
  for(i in seq_along(us)){
    print(paste("Bootsrap for u =",us[i]))
    thetas <- matrix(NA,ncol=3,nrow=n.samp)
    for(j in 1:n.samp){
      X.boot <- blockbootstrap(X,block.size)
      thetas[j,1] <- exi(X.boot,u=us[i],r=1)
      thetas[j,2] <- exi(X.boot,u=us[i],r=2)
      thetas[j,3] <- exi(X.boot,u=us[i],r=0)
    }
    CIs[i,c(1,2)] <- quantile(thetas[,1],c(0.025,0.975))
    CIs[i,c(3,4)] <- quantile(thetas[,2],c(0.025,0.975))
    CIs[i,c(5,6)] <- quantile(thetas[,3],c(0.025,0.975))
  }
  CIs
}

library(evd)

set.seed(333)
ks <- c(2:15)
thetas <- matrix(NA,ncol=3,nrow=length(ks))
for(i in 1:length(ks)){
  X <- ChernickAR(ks[i],1000)
  theta <- exi(X,u=0.95,r=1)
  CIs <- theta.bootstrap(X,n.samp=250,us=c(0.95),block.size=50)
  thetas[i,] <- c(2*theta-CIs[,2],theta,2*theta-CIs[,1])
}

par(mfrow=c(1,1),mar=c(4.1,6.1,0.5,2.1))
plot(ks,thetas[,2],type="l",ylim=c(0.2,1),
     xlab="XXX",ylab="YYY",cex.main=2,
     cex.lab=2,cex.axis=1.5,lwd=3)
lines(ks,thetas[,1],lty="dotted")
lines(ks,thetas[,3],lty="dotted")
lines(ks,(ks-1)/ks,col="red")

