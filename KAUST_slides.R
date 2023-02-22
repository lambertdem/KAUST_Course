rm(list=ls())
set.seed(444)
# ---------------------------- Slide 16 -----------------------------

# Sample 100 iid standard Frechet RVs
y <- rfrechet(100,shape=1)

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
  
  y <- rfrechet(n,shape=1)
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


dts <- c(1914:1961)
dts <- rep(dts,each=365)
dts <- c(dts,rep(1962,11))
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

myd <- read.table(
  "C:/Users/lambe/Documents/McGill/Masters/Courses/MATH80622A/A4/SP500-10d.csv", 
  header=TRUE, sep=",")

par(mfrow=c(1,1),mgp=c(3,0.8,0),mar=c(5,5,0.5,0.5))
plot(as.Date(myd$Date),myd$Close,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

rets <- c()
for(i in 1:(length(myd$Close)-1)){
  rets <- c(rets,(myd$Close[i]-myd$Close[i+1])/myd$Close[i])
}
rets

# ---------------------------- Slide 42 -----------------------------

rets.2 <- rets^2
plot(as.Date(myd$Date)[1:(length(myd$Date)-1)],rets.2,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

# ---------------------------- Slide 43 -----------------------------

thetas.sp500 <- c()
qs <- seq(0.0001,0.0016,by=0.00005)
for(q in qs){
  thetas.sp500 <- c(thetas.sp500,exi(rets.2,u=q,r=10))
}
plot(qs,thetas.sp500,ylim=c(0,1),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")
abline(v=0.0003,lty="dotted")
mean.theta <- mean(thetas.sp500[5:length(thetas.sp500)])
lines(cbind(c(0.0003,0.002),rep(mean.theta,2)),col="red")

# ---------------------------- Slide 44 -----------------------------

h <- 20
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
plot(means)
plot(sds)
subseq.dts <- as.Date(myd$Date)[(1+h/2):(length(myd$Date)-(1+h/2))]
plot(subseq.dts,res^2,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")

# ---------------------------- Slide 45 -----------------------------

thetas.res1 <- c()
thetas.res2 <- c()
qs <- seq(0,6,by=0.01)
for(q in qs){
  thetas.res1 <- c(thetas.res1,exi(res^2,u=q,r=1))
  thetas.res2 <- c(thetas.res2,exi(res^2,u=q,r=2))
}
plot(qs,thetas.res1,ylim=c(0,1),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l",lwd=2)
lines(qs,thetas.res2,col="blue",lwd=2)
legend(4,0.4,legend=c("m=1","m=2"),col=c("black","blue"),
       lty=c(1,1),lwd = c(3,3))

# ---------------------------- Slide 47 -----------------------------

cs <- clusters(rets.2,u=0.0004,r=10)
maxima <- c()
ind.max <- c()
for(i in 1:length(cs)){
  maxima <- c(maxima,max(cs[[i]]))
  
  ind <- as.numeric(names(cs[[i]]))[which.max(cs[[i]])]
  ind.max <- c(ind.max,ind)
}
maxima
ind.max

quantile(rets.2,0.95)
length(cs)/length(rets.2[rets.2>0.0004])

dts <- as.Date(myd$Date)[1:(length(myd$Date)-1)]
plot(dts,rets.2,ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,type="l")
ind.c <- which(rets.2>0.0004)
points(dts[ind.c],rets.2[ind.c],col="blue",pch=20)
points(dts[ind.max],rets.2[ind.max],col="red",pch=20)

as.numeric(names(cs[[4]]))

exi(rets.2,u=quantile(rets.2,0.95),r=10)
library(ismev)
u <- 0.0004
fit <- gpd.fit(maxima,threshold=u)
fit$mle
fit$se
diag <- gpd.diag(fit)


sig <- fit$mle[1]*140^fit$mle[2]

mu <- u +

ps <- c(1:35)/36

par(mfrow=c(1,2),mgp=c(3,0.8,0),mar=c(5,5,2.5,0.5))
plot(ps,sort(pgpd(maxima,u,fit$mle[1],fit$mle[2])),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,cex.main=2,main="ZZZZZ")
abline(a=0,b=1)

qs <- qgpd(ps,u,pars[1],pars[2])
plot(qs,sort(maxima),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,cex.main=2,main="WWWWW")
abline(a=0,b=1)
