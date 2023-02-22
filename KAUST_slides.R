rm(list=ls())
set.seed(444)
# ---------------------------- Slide 16 -----------------------------

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
qs <- c(4:35)
for(q in qs){
  thetas1 <- c(thetas1,exi(rain,u=q,r=1))
  thetas2 <- c(thetas2,exi(rain,u=q,r=2))
}


par(mfrow=c(1,1),mar=c(4.1,6.1,0.5,2.1))
plot(qs,thetas1,type="l",ylim=c(0.3,1.05),
     xlab="u",ylab="theta",cex.main=2,
     cex.lab=2,cex.axis=1.5)
lines(qs,thetas2,col="blue")
legend(25,0.6,legend=c("m=1","m=2"),col=c("black","blue"),
       lty=c(1,1),lwd = c(2,2))
abline(h=1,lty="dotted")

# ---------------------------- Slide 41 -----------------------------

# Set current wd to KAUST_Course
# setwd(".../KAUST_Course/")
myd <- read.table("Data/SP500-10d.csv", header=TRUE, sep=",")

par(mfrow=c(1,1),mgp=c(3,0.8,0),mar=c(5,5,0.5,0.5))
plot(as.Date(myd$Date),myd$Close,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

rets <- c()
for(i in 1:(length(myd$Close)-1)){
  rets <- c(rets,(myd$Close[i]-myd$Close[i+1])/myd$Close[i])
}
plot(as.Date(myd$Date)[1:(length(myd$Date)-1)],rets,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

# ---------------------------- Slide 42 -----------------------------

# Get the squared returns
rets.2 <- rets^2

plot(as.Date(myd$Date)[1:(length(myd$Date)-1)],rets.2,ylab="YYYYY",xlab="XXXXX",
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
subseq.dts <- as.Date(myd$Date)[(1+h/2):(length(myd$Date)-(1+h/2))]
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

# Store index of cluster maxima within all squred returns
ind.max <- c()

for(i in 1:length(cs)){
  # Get maxima of cluster i
  maxima <- c(maxima,max(cs[[i]]))
  
  # Get index of all values in cluster i
  ind.clust.i <- as.numeric(names(cs[[34]]))
  
  # Select index corresponding to maxima within cluster i
  ind <- ind.clust.i[which.max(cs[[i]])]
  
  ind.max <- c(ind.max,ind)
}

# Get index of values exceeding u (i.e. those in clusters)
ind.c <- which(rets.2>u)

dts <- as.Date(myd$Date)[1:(length(myd$Date)-1)]
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

