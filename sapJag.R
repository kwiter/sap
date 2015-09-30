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