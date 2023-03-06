rm(list=ls())
library(evd)
setwd("C:/Users/lambe/Documents/Edinburgh/PhD/CodeSharing_IP_LDM/KAUST/R_exercise_solutions/")
MTL.weather <- read.csv("MTL_1941_2020.csv")

# Subsample only summer data
summer.weather <- MTL.weather[MTL.weather$Month %in% c(6,7,8),]

# Select only date and maximum daily temperature columns
cols <- which(colnames(summer.weather)%in%c("Date.Time","Max.Temp...C."))
summer.temp <- summer.weather[,cols]
colnames(summer.temp) <- c("Date","max.daily.temp")
summer.temp <- summer.temp[which(!is.na(summer.temp$max.daily.temp)),]

# Plot of Montreal summer temperature data
par(mfrow=c(1,1),mar=c(4.1,6.1,0.5,2.1))
plot(as.Date(summer.temp$Date),summer.temp$max.daily.temp,
     xlab="XXX",ylab="YYY",cex.main=2,
     cex.lab=2,cex.axis=1.5)

X <- summer.temp$max.daily.temp

# Sequence of quantiles of the data 
qs <- seq(0.95,0.995,by=0.001)
us <- quantile(summer.temp$max.daily.temp,qs)

# Sequence of run length
rs <- c(1:15)

# Matrix of thetas for specified threshold an run length
thetas <- matrix(NA,nrow=length(us),ncol=length(rs))
for(i in seq_along(us)){
  for(j in seq_along(rs)){
    thetas[i,j] <- exi(X,u=us[i],r=rs[j])
  }
}

# Plot of theta estimates at different threshold for each run length
plot(us,thetas[,1],type="l",ylim=c(0,1))
for(i in 2:dim(thetas)[2]){
  lines(us,thetas[,i])
}
abline(v=us[25],lty="dotted")

# Identifying the clusters
cs <- clusters(X,u=us[25],r=15)

# Store the maxima from each cluster
maxima <- c()

for(i in 1:length(cs)){
  # Get maxima of cluster i
  maxima <- c(maxima,max(cs[[i]]))
}

# Fit a GPD to cluster maxima
library(ismev)
u <- us[25]
fit <- gpd.fit(maxima,threshold=u)
diag <- gpd.diag(fit)

sig.u <- fit$mle[1]
xi <- fit$mle[2]

n.exc.u <- length(X[X>u])

sig <- sig.u*n.exc.u^xi
mu <- u + (sig-sig.u)/xi

sig
mu

# PP plot of GPD fit on maxima of clusters of squared standardised returns
ps <- c(1:length(cs))/(length(cs)+1)
par(mfrow=c(1,2),mgp=c(3,0.8,0),mar=c(5,5,2.5,0.5))
plot(ps,sort(pgpd(maxima,u,fit$mle[1],fit$mle[2])),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,cex.main=2,main="ZZZZZ")
abline(a=0,b=1)

# QQ plot of GPD fit on maxima of clusters of squared standardised returns
qs <- qgpd(ps,u,fit$mle[1],fit$mle[2])
plot(qs,sort(maxima),ylab="YYYYY",xlab="XXXXX",
     cex.lab=2,cex.axis=1.5,cex.main=2,main="WWWWW")
abline(a=0,b=1)
