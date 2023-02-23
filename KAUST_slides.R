rm(list=ls())
set.seed(444)
# ---------------------------- Slide 16 -----------------------------
library(evd)
# Sample 100 iid standard Frechet RVs
par(mfrow=c(1,1),mar=c(4.1,6.1,0.5,2.1))

y <- rgev(100,0,1,1)

# Generate the process X_t
max.ind <- c()
x <- c(y[1])
for(i in 3:length(y)){
  x <- c(x,max(c(y[i],y[i-1])))
  if(max(c(y[i],y[i-1]))==y[i-1]){
    max.ind <- c(max.ind,i-1)
  }
}

par(mfrow=c(2,1),mar=c(4.1,6.1,0.5,2.1))
plot(y,cex.lab=4,xlab="",cex.axis=1.3)
plot(x,cex.lab=4,xlab="t",cex.axis=1.3)
points(max.ind,x[max.ind],col="red")

# ---------------------------- Slide 29 -----------------------------
par(mfrow=c(1,4),mgp=c(3,0.8,0),mar=c(5,5,2.5,0.6))
alpha <- c(1/3,1/2,1/6)
ns <- c(10,100,1000,10000)

# For samples sizes n=10,100,..., sample n standard Frechet and
#   and produce the MovingMax(alpha) process
for(k in seq_along(ns)){
  n <- ns[k]+2
  
  y <- rgev(n,0,1,1)
  x <- c()
  for(i in 3:length(y)){
    x <- c(x,sum(y[(i-2):i]*alpha))
  }
  #expression(X[t]/n)
  if(k==1){
    plot(c(1:length(x))/(length(x)+1),x/length(x),ylab="YYYYY",main=paste0("n=",length(x)),
         xlab="XXXXX",cex.lab=2,cex.axis=1.4,cex.main=2)
  }
  else{
    plot(c(1:length(x))/(length(x)+1),x/length(x),ylab="",main=paste0("n=",length(x)),
         xlab="XXXXX",cex.lab=2,cex.axis=1.4,cex.main=2)
  }
}

# ---------------------------- Slide 39 -----------------------------

library(ismev)

# Obtain the rainfall data set
data("rain")


dts <- rep(c(1914:1961),each=365)
dts <- c(dts,rep(1962,length(rain)-length(dts)))

par(mfrow=c(1,1),mar=c(4.1,6.1,0.5,2.1))
plot(rain,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,xaxt="n")
axis(1,at=seq(1,(1962-1914)*365,by=365*4),labels=seq(1914,1961,by=4))

# ---------------------------- Slide 40 -----------------------------

library(evd)

thetas1 <- c()
thetas2 <- c()
thetas3 <- c()
us <- c(4:35)
for(u in us){
  thetas1 <- c(thetas1,exi(rain,u=u,r=1))
  thetas2 <- c(thetas2,exi(rain,u=u,r=2))
  thetas3 <- c(thetas3,exi(rain,u=u,r=0))
}

par(mfrow=c(1,1),mar=c(4.1,6.1,0.5,2.1))
plot(us,thetas1,type="l",ylim=c(0.3,1.05),
     xlab="u",ylab="theta",cex.main=2,
     cex.lab=2,cex.axis=1.5,lwd=2)
lines(us,thetas2,col="blue",lwd=2)
lines(us,thetas3,col="red",lwd=2)
legend(25,0.6,legend=c("m=1","m=2"),col=c("black","blue"),
       lty=c(1,1),lwd = c(2,2))
abline(h=1,lty="dotted")

# Function to obtain moving block bootstrap from data set X
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

# ACF of rain data to select block size
acf(rain)

# Pick 47 since length(rain)/47 is integer
length(rain)/47

CIs <- theta.bootstrap(rain,n.samp,us,47)

par(mfrow=c(1,2),mar=c(4.1,6.1,2,0.2))
plot(us,thetas1,type="l",ylim=c(0.3,1.05),
     xlab="u",ylab="theta",main="WWWW",
     cex.main=2,cex.lab=2,cex.axis=1.5,lwd=2)
polygon(x=c(us,rev(us)),y=c(CIs[,1],rev(CIs[,2])),
        col="lightgrey",border="lightgrey")
polygon(x=c(us,rev(us)),y=c(CIs[,3],rev(CIs[,4])),
        col="lightgrey",border="lightgrey")
lines(us,thetas1,lwd=2)
lines(us,thetas2,lwd=2)
abline(h=1,lty="dotted")


plot(us,thetas3,type="l",ylim=c(0.3,1.05),
     xlab="u",ylab="theta",main="YYYYYYYYYYYYY",
     cex.main=2,cex.lab=2,cex.axis=1.5,lwd=2)
polygon(x=c(us,rev(us)),y=c(CIs[,5],rev(CIs[,6])),col="lightgrey",border="lightgrey")
lines(us,thetas3,lwd=2)
abline(h=1,lty="dotted")

# ---------------------------- Slide 41 -----------------------------

# Set current wd to KAUST_Course
# setwd(".../KAUST_Course/")
setwd("C:/Users/lambe/Documents/Edinburgh/PhD/KAUST/KAUST_Course/")
SP500 <- read.table("Data/SP500_daily.csv", header=TRUE, sep=",")

par(mfrow=c(1,1),mgp=c(3,0.8,0),mar=c(5,5,0.5,0.5))
plot(as.Date(SP500$Date),SP500$Close,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

rets <- c()
for(i in 1:(length(SP500$Close)-1)){
  rets <- c(rets,(SP500$Close[i]-SP500$Close[i+1])/SP500$Close[i])
}
plot(as.Date(SP500$Date)[1:(length(SP500$Date)-1)],rets,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

# ---------------------------- Slide 42 -----------------------------

# Get the squared returns
rets.2 <- rets^2

plot(as.Date(SP500$Date)[1:(length(SP500$Date)-1)],rets.2,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

# ---------------------------- Slide 43 -----------------------------

us <- seq(0.0001,0.0016,by=0.00005)

# Compute the extremal index for a range of u's (us)
fexi <- function(u,returns,r) exi(returns,u=u,r=r)
thetas.sp500 <- unlist(lapply(us,fexi,returns = rets.2,r=10))

plot(us,thetas.sp500,ylim=c(0,1),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")
abline(v=0.0004,lty="dotted")

# Display the mean of thetas above 0.0003
mean.theta <- mean(thetas.sp500[5:length(thetas.sp500)])
lines(cbind(c(0.0004,0.002),rep(mean.theta,2)),col="red")

# ---------------------------- Slide 44 -----------------------------

# Moving average length
h <- 20

# Store moving averages, moving sd, residuals
means <- c()
sds <- c()
res <- c()
for(i in (1+h/2):(length(rets)-h/2)){
  m <- mean(rets[(i-h/2):(i+h/2)])
  means <- c(means,m)
  s <- sd(rets[(i-h/2):(i+h/2)])
  sds <- c(sds,s)
  res <- c(res,(rets[i]-m)/s)
}

# Compute standardised returns
res <- (rets[(h/2+1):(length(rets)-h/2)]-means)/sds

# Plot the squared standardised returns 
subseq.dts <- as.Date(SP500$Date)[(1+h/2):(length(SP500$Date)-(1+h/2))]
plot(subseq.dts,res^2,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

# ---------------------------- Slide 45 -----------------------------

# Compute extremal indexes for m=1,2 and various u's (us)
thetas.res1 <- c()
thetas.res2 <- c()
us <- seq(0,6,by=0.01)
for(u in us){
  thetas.res1 <- c(thetas.res1,exi(res^2,u=u,r=1))
  thetas.res2 <- c(thetas.res2,exi(res^2,u=u,r=2))
}

# Plot extremal indexes computed for m=1,2 and various u's (us)
par(mfrow=c(1,1),mar=c(4.1,6.1,2,0.2))
plot(us,thetas.res1,ylim=c(0,1),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l",lwd=2)
lines(us,thetas.res2,col="blue",lwd=2)
legend(4,0.4,legend=c("m=1","m=2"),col=c("black","blue"),
       lty=c(1,1),lwd = c(3,3))


# ---------------------------- Slide 47 -----------------------------

# Obtain list of clusters defined for u=0.0003
u <- 0.0004
cs <- clusters(rets.2,u=u,r=10)
length(cs)

# Store the maxima from each cluster
maxima <- c()

# Store index of cluster maxima within all squared returns
ind.max <- c()

for(i in 1:length(cs)){
  # Get maxima of cluster i
  maxima <- c(maxima,max(cs[[i]]))
  
  # Get index of all values in cluster i
  ind.clust.i <- as.numeric(names(cs[[i]]))
  
  # Select index corresponding to maxima within cluster i
  ind <- ind.clust.i[which.max(cs[[i]])]
  
  ind.max <- c(ind.max,ind)
}

# Get index of values exceeding u (i.e. those in clusters)
ind.c <- which(rets.2>u)

dts <- as.Date(SP500$Date)[1:(length(SP500$Date)-1)]
plot(dts,rets.2,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")
points(dts[ind.c],rets.2[ind.c],col="blue",pch=20)
points(dts[ind.max],rets.2[ind.max],col="red",pch=20)

# ---------------------------- Slide 49 -----------------------------

library(ismev)
fit <- gpd.fit(maxima,threshold=u)
diag <- gpd.diag(fit)

sig.u <- fit$mle[1]
xi <- fit$mle[2]

n.exc.u <- length(rets.2[rets.2>u])

sig <- sig.u*n.exc.u^xi
mu <- u + (sig-sig.u)/xi

sig
mu

# ---------------------------- Slide 50 -----------------------------

ps <- c(1:length(cs))/(length(cs)+1)

# PP plot of GPD fit on maxima of clusters of squared standardised returns
par(mfrow=c(1,2),mgp=c(3,0.8,0),mar=c(5,5,2.5,0.5))
plot(ps,sort(pgpd(maxima,u,fit$mle[1],fit$mle[2])),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,cex.main=2,main="ZZZZZ")
abline(a=0,b=1)

# QQ plot of GPD fit on maxima of clusters of squared standardised returns
qs <- qgpd(ps,u,fit$mle[1],fit$mle[2])
plot(qs,sort(maxima),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,cex.main=2,main="WWWWW")
abline(a=0,b=1)
