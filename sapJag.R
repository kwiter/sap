rm(list = ls()) 
setwd('C:/Users/mck14/Dropbox/warming/ENVDATA/processedData/DukeEnv')
setwd('/home/mckwit/Dropbox/warming/ENVDATA/processedData/DukeEnv')
sap = read.csv('sap.csv',header=F)

par(cex=.7,pch=20,bty = 'n',bg = 'white',col = nColor()[1],
    col.axis = "#696969", col.lab = "#696969", col.main = "#696969")
palette(nColor(10,distinct = T))
palette(c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951"))
thecol=c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
thecol=c("#00A0B036","#6A4A3C36","#CC333F36","#EB684136","#EDC95136")

##working paths
if(length(grep('mck14',getwd()))==1){
  hdir <- 'C:/Users/mck14/Dropbox/'
  setwd(hdir)
}
if(length(grep('mckwit',getwd()))==1){
  hdir <- '/home/mckwit/Dropbox/' #home
  setwd(hdir)
}
source('mckFUNCTIONS/mckFUNCTIONS.r')
require('chron')


load(file = paste(hdir,'warming/AssimIndiv.RData',sep=''))
pheno = read.csv(file=paste(hdir,'DissertationEssentials/Biomass/BiomassPhen.csv',sep=''),stringsAsFactors = F)
species = sort(species)
load(file='warming/ENVDATA/processedData/sapENV.RData') #Q,vpd,SP,AT,RH,SM

days = readRDS(file='warming/ENVDATA/processedData/days.rds')
#dayDec = a.n(colnames(SF)) - 365 - 365
#whr = which(dayDec - floor(dayDec) == 0 )  #day starts]
#whr24 = c(1,which(diff(days) == 1)+1)
times = a.n(sap[1,])
times = times[-1]
plots = sap[,1]
plots = plots[-1]
sap = sap[-1,]
sap = t(sap)
colnames(sap) = plots
sap = sap[-1,]
rownames(sap) = times
write.csv(sap,file=paste0(hdir,'sap.csv'))
sap = read.csv(file=paste0(hdir,'sap.csv'),row.names=1)

for(i in 1:ncol(sap)){
  tmp = a.n(sap[,i])
  sap[,i] = round(rowMeans(cbind(c(tmp[-1],0) , c(tmp[-c(1,2)],0,0) , tmp , c(0,tmp[-len(tmp)]) , c(0,0,tmp[-c(len(tmp),len(tmp)-1)])),na.rm=F),5)
  print(i)
}

whrSamp = seq(1,nrow(sap),by=4)
sap = sap[whrSamp,]
Q = Q[whrSamp,]
RH = RH[whrSamp,]
vpd = vpd[whrSamp,]
SM = SM[whrSamp,]
AT = AT[whrSamp,]
SP = SP[whrSamp,]

sap[sap==0] = NA
colnames(sap)[115] = "g12.17.1.qual"
write.csv(sap,file=paste0(hdir,'hourSap.csv'))
sap = read.csv(paste0(hdir,'hourSap.csv'),row.names=1)


plots = TtoC(colnames(sap),sep="\\.")[,1]
heat = a.n(treatment(plots,combineAandC=T) == 'H')
gap = a.n(grepl('g',colnames(sap),ignore.case=T))  #1 is gap
chamb = a.n(treatment(plots) != 'control') #1 is chamber
specs = facTOdummy(TtoC(colnames(sap),sep="\\.")[,4])
whrEnv = match(plots,STDplots(gsub('Q',"",colnames(SM)[grep('Q',colnames(SM))]))) + 8
n = nrow(sap)
xmat = y = numeric()
for(i in 1:ncol(sap)){
  y = c(y,sap[,i])
  xmat = rbind(xmat,
               cbind(
                    rep(1,n),
                    rep(heat[i],n),
                    rep(gap[i],n),
                    rep(chamb[i],n),
                    matrix(rep(specs[i,],n),ncol=3,byrow=T),
                    AT[,whrEnv[i]],
                    Q[,whrEnv[i]],
                    vpd[,whrEnv[i]],
                    SM[,whrEnv[i]]
                    )
               )
  print(i)
}
colnames(xmat) = c('int','hot','gap','cham','litu','qual','quru','at','q','vpd','sm')
whrS = which(!is.na(y))
nZ = len(whrS)
nI = ncol(sap)
sts = seq(1,len(y),by = n)

library('truncnorm')

N <- 1000
n = 100
b <- c(.1, -.1)

sig <- 1
sig2 <- 5

p <- length(b)
x <- cbind(1, matrix(runif((p-1)*N),N, p-1))
#dy <- rbind(0, x[-N,]%*%b) + rtruncnorm(N, 0, Inf, 0, sig)
dy <- rbind(0, x[-N,]%*%b) + rnorm(N, 0, sig)
y <- cumsum(dy)
whrS= sort(sample(1:1000,n,replace=T))
whrS = c(whrS,whrS)
z <- rnorm(length(whrS), y[whrS], sig2)
plot(y)

example1.bug <-
  'model {

    for (i in 1:n) {
      z[whr[i]] ~ dnorm(y[whr[i]], tau2)
    }

    tau2 <- pow(sigma2, -2)
    sigma2 ~ dunif(0, 100)

  y[1] ~ dnorm( 0, tau)
  for (i in 2:(N)) {
    y[i] ~ dnorm(y[i-1] + x[i-1,]%*%beta , tau)
  }

  y[1] ~ dnorm( y[2] - x[1,]%*%beta , tau)
  for (i in 2:(N-1)) {
    y[i] ~ dnorm((y[i-1] + x[i-1,]%*%beta + y[i+1] - x[i,]%*%beta)/2, tau)
  }
  y[N] ~ dnorm(y[N-1] + x[N-1,]%*%beta, tau)

  for (i in 1:p) {
    beta[i] ~ dnorm(0, .0001)
  }

  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)

  }'

#y[1] ~ dnorm( y[2] - x[1,]%*%beta , tau)
#for (i in 2:(N-1)) {
#  y[i] ~ dnorm((y[i-1] + x[i-1,]%*%beta + y[i+1] - x[i,]%*%beta)/2, tau)
#}
#y[N] ~ dnorm(y[N-1] + x[N-1,]%*%beta, tau)

#y[1] ~ dnorm( 0, tau)
#for (i in 2:(N)) {
#  y[i] ~ dnorm(y[i-1] + x[i-1,]%*%beta , tau)
#}


library('rjags')

jags <- jags.model(textConnection(example1.bug),
                   data = list('x' = x,
                               'z' = z,
                               'N' = N,
                               'n' = n,
                               'whr' = whrS,
                               'p'= ncol(x)),
                   n.chains = 4,
                   n.adapt = 100)

example1.bug <-
  'model {
    
  for (i in 1:n) {
      z[whr[i]] ~ dnorm(y[whr[i]], tau2)
    }
    tau2 <- pow(sigma2, -2)
    sigma2 ~ dunif(0, 100)

  for(j in 1:nI) {
    y[sts[j]] ~ dnorm( 0, tau)
    for (i in (1+sts[j]):(sts[j] + N - 1)) {
      y[i] ~ dnorm(y[i-1] + x[i-1,]%*%beta , tau)
    }
  }

  for (i in 1:p) {
    beta[i] ~ dnorm(0, .0001)
  }

  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)

  }'
library(rjags)
jags <- jags.model(textConnection(example1.bug),
                   data = list('x' = xmat,
                               'z' = y,
                               'N' = n,
                               'nI' = nI,
                               'n' = nZ,
                               'whr' = whrS,
                               'p'= ncol(xmat),
                               'sts' = sts),
                            
                   n.chains = 2,
                   n.adapt = 100)

update(jags, 1000)

tmp <- jags.samples(jags,
                    c('y','beta', 'sigma', 'sigma2'),
                    1000)

print(tmp)

plot(y,type='l')
points(whrS,z,col=2)
lines(apply(tmp$y,1,mean),type='l',col=4)
lines(apply(tmp$y,1,quantile,.975),col=3,lty=3)
lines(apply(tmp$y,1,quantile,.0275),col=5,lty=3)

plot(tmp$beta)
plot(tmp$sigma)
plot(tmp$sigma2)

plot(apply(tmp$beta[1,,],1,mean))
plot(apply(tmp$beta[2,,],1,mean))
tmp$beta

plot(tmp$beta)

plot(tmp$sigma)



















## multilevel hockey stick model
## 
hockeyStick <- '
model{
  for (i in 1:n){
    y[i] ~ dnorm(y.hat[i], tau.y[stid[i]])
    y.hat[i] <- beta0[stid[i]]+(beta1[stid[i]]-delta[stid[i]]*step(x[i]-phi[stid[i]]))*(x[i]-phi[stid[i]])
  }
  for (j in 1:n.st){
    beta0[j] ~ dnorm(mu0, tau[1])
    beta1[j] ~ dnorm(mu1, tau[2])
    delta[j] ~ dnorm(mu2, tau[3])
    phi[j] ~ dnorm(mu3,tau[4])T(L[j],U[j])
    tau.y[j] <- pow(sigma.y[j], -2)
    sigma.y[j] ~ dunif(0, 10)
  } 
  mu0 ~ dnorm(0, 0.0001)
  mu1 ~ dnorm(0, 0.0001)
  mu2 ~ dnorm(0, 0.0001)
  mu3 ~ dnorm(0, 0.0001)
  
  for (k in 1:4){
    tau[k] <- pow(sigma[k], -2)    
    sigma[k] ~ dunif(0, 100)
  }
}
'

x = seq(1,30,by=1) - 15.5
n = len(x)
nst = 4
stid = rep(1:4, each = n)
b0 = 10; b1 = .5; de = 10; ph = 15.5; sig = 2
y = rep(0,nst*n)
y[1:n]              = c(b0 + b1*x[1:15],b0 + (b1+de)*x[16:30]) + 10
y[(n+1):(2*n)]      = c(b0 + b1*x[1:15],b0 + (b1+de)*x[16:30]) - 10
y[((2*n)+1):(3*n)]  = c(b0 + b1*x[1:15],b0 + (b1+de)*x[16:30]) + 20
y[((3*n)+1):(4*n)]  = c(b0 + b1*x[1:15],b0 + (b1+de)*x[16:30]) - 20
y = rnorm(len(y),y,sig ) 
x = rep(x,4)
L <- tapply(x, stid, min)
U <- tapply(x, stid, max)
plot(x,y)

#bugs.dat <- list(n=n, nst=n.st, y=y, x=x, stid=stid,              
#                 minYr=min(L), maxYr=max(U))

dat <- list(
  n = len(y),
  n.st = nst,
  y = y,
  stid = stid,
  x =  x,
  L = L,
  U = U
)

parameters <- c("y.hat","beta0","beta1","delta", "phi", "mu0", 
                "mu1","mu2","mu3","sigma", "sigma.y")

library(rjags)
jags <- jags.model(textConnection(hockeyStick),
                   data =dat,
                   n.chains = 2,
                   n.adapt = 100
                   )

update(jags, 1000)

tmp <- jags.samples(jags,
                    parameters,
                    1000)
summary(tmp)
tmp$phi

plot(x,y,type='p')
points(x,apply(tmp$y,1,mean),col=4,pch=20)
for(i in 1:len(x)){
  lines(c(x[i],x[i]),c(apply(tmp$y,1,quantile,.975)[i],apply(tmp$y,1,quantile,.025)[i]),col=3,lty=1)
}





## multilevel hockey stick model
## 
hockeyStick <- '
model{
  for (i in 1:n){
    y[i] ~ dnorm(y.hat[i], tau.y[stid[i]])
    y.hat[i] <- beta0[stid[i]]+(beta1[stid[i]]-delta[stid[i]]*step(x[i]-phi[stid[i]]))*(x[i]-phi[stid[i]])
  }
  for (j in 1:n.st){
    beta0[j] ~ dnorm(mu0, tau[1])
    beta1[j] ~ dnorm(mu1, tau[2])
    delta[j] ~ dnorm(mu2, tau[3])
    phi[j] ~ dnorm(mu3,tau[4])T(L[j],U[j])
    tau.y[j] <- pow(sigma.y[j], -2)
    sigma.y[j] ~ dunif(0, 10)
  } 
  mu0 ~ dnorm(0, 0.0001)
  mu1 ~ dnorm(0, 0.0001)
  mu2 ~ dnorm(0, 0.0001)
  mu3 ~ dnorm(0, 0.0001)
  
  for (k in 1:4){
    tau[k] <- pow(sigma[k], -2)    
    sigma[k] ~ dunif(0, 100)
  }
}
'

quantReg <- 
'
model{
 for(i in 1:n){
   mu[i] <- alpha + beta*x[i]
   w[i] ~ dexp(tau)
   me[i] <- (1-2*p)/(p*(1-p))*w[i] + mu[i]
   pe[i] <- (p*(1-p)*tau)/(2*w[i])
   y[i] ~ dnorm(me[i],pe[i])
 }

 #priors for regression
 alpha ~ dnorm(0,1E-6)
 beta ~ dnorm(0,1E-6)

 lsigma ~ dunif(-5,15)
 sigma <- exp(lsigma/2)
 tau <- pow(sigma,-2)
}
'

dat <- list(
  n = len(y),
  y = y,
  x =  x,
  p = .95
)

parameters <- c("alpha","beta","mu")

library(rjags)
jags <- jags.model(textConnection(quantReg),
                   data =dat,
                   n.chains = 2,
                   n.adapt = 100
)

update(jags, 1000)

tmp <- jags.samples(jags,
                    parameters,
                    1000)
summary(tmp)
tmp$phi

plot(x,y,type='p')
points(x,apply(tmp$mu,1,mean),col=4,pch=20)
for(i in 1:len(x)){
  lines(c(x[i],x[i]),c(apply(tmp$mu,1,quantile,.975)[i],apply(tmp$mu,1,quantile,.025)[i]),col=3,lty=1)
}




x = seq(1,30,by=1) -15.5
n = len(x)
nst = 4
stid = rep(1:4, each = n)
b0 = 10; b1 = .5; de = 10; ph = 15.5; sig = 10
y = rep(0,nst*n)
y[1:n]              = c(b0 + b1*x[1:15],b0 + (b1+de)*x[16:30]) + 10
y[(n+1):(2*n)]      = c(b0 + b1*x[1:15],b0 + (b1+de)*x[16:30]) - 10
x = seq(1,30,by=1) -20.5
y[((2*n)+1):(3*n)]  = c(b0 + (20*b1)*x[1:20],b0 + ((b1*20)-20)*x[21:30]) + 20
y[((3*n)+1):(4*n)]  = c(b0 + (20*b1)*x[1:20],b0 + ((b1*20)-20)*x[21:30]) - 20
y = rnorm(len(y),y,sig ) 
x = c(rep(seq(1,30,by=1) -15.5,2),rep(seq(1,30,by=1) -20.5,2))
plot(x,y)
x = rep(seq(1,30,by=1),4)
L <- tapply(x, stid, min)
U <- tapply(x, stid, max)


## multilevel hockey stick model
## 
stickQuant <- '
model{
  for (i in 1:n){
    x.hat[i] <- x[i] - phi[stid[i]]
    w[i] ~ dexp(tau.y[stid[i]])
    me[i] <- (1-2*p)/(p*(1-p))*w[i] + y.hat[i]
    pe[i] <- (p*(1-p)*tau.y[stid[i]])/(2*w[i])
    y[i] ~ dnorm(me[i],pe[i])
    y.hat[i] <- beta0[stid[i]]+(beta1[stid[i]]-delta[stid[i]]*step(x.hat[i]-phi[stid[i]]))*(x.hat[i]-phi[stid[i]])
  }
  for (j in 1:n.st){
    beta0[j] ~ dnorm(mu0, tau[1])
    beta1[j] ~ dnorm(mu1, tau[2])
    delta[j] ~ dnorm(mu2, tau[3])
    phi[j] ~ dnorm(mu3,tau[4])T(L[j],U[j])
    tau.y[j] <- pow(sigma.y[j], -2)
    sigma.y[j] ~ dunif(0, 10)
  } 
  mu0 ~ dnorm(0, 0.0001)
  mu1 ~ dnorm(0, 0.0001)
  mu2 ~ dnorm(0, 0.0001)
  mu3 ~ dnorm(0, 0.0001)

  for (k in 1:4){
    tau[k] <- pow(sigma[k], -2)    
    sigma[k] ~ dunif(0, 100)
  }
}
'


dat <- list(
  n = len(y),
  n.st = nst,
  y = y,
  stid = stid,
  x =  rep(seq(1,30,by=1),4) ,
  L = L,
  U = U,
  p = .25
)

parameters <- c("y.hat","beta0","beta1","delta", "phi", "mu0", 
                "mu1","mu2","mu3","sigma", "sigma.y")

library(rjags)
jags <- jags.model(textConnection(stickQuant),
                   data =dat,
                   n.chains = 4,
                   n.adapt = 1000
)

update(jags, 3000)

tmp <- jags.samples(jags,
                    parameters,
                    3000)
summary(tmp)
tmp$phi
x = rep(seq(1,30,by=1),4) 
plot(x,y,type='p')
points(x,apply(tmp$y,1,mean),col=4,pch=20)
hi = apply(tmp$y,1,quantile,.975)
lo = apply(tmp$y,1,quantile,.025)       
for(i in 1:len(x)){
  lines(c(x[i],x[i]),c(hi[i],lo[i]),col=3,lty=1)
}

