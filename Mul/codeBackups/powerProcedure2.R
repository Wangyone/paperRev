rm(list = ls())
######################### Type I Error or Power #########################
setwd('X:/Desktop/JointCUSUM/Mul')
library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
sourceCpp('statLMul.cpp')
sourceCpp('eucLMul2.cpp')
source('dataGenerationSupply.R')
## boundary function 
Qk<-function(s0,gamma0){
  u=(1+s0)*(s0/(1+s0))^gamma0
}
############################  Model para
d <- 5 # dimension of data
DGP <- 1
lag <- 0
kernelTunta <- 3  
## 1 for gauss  && 2 for energy && 3 for quadexp
####################  data length 
trainingSize <- 200
monitoringPeriod <- 1
## para for monitomg statistics
gamma<-0
deltagauss <- 1*d
nsim<- 1000


yt <- dataGeneration(monitoringPeriod,trainingSize,model = DGP,dataDimension = d)
# MTSplot(yt)


E <- eucL2(yt[1:trainingSize,],lag)
medianE <- (median(E[E > 0]))
rm(E)

alphaE <- 1  ##energy
if (kernelTunta == 1||kernelTunta == 2) {
  delta0 <- deltagauss  ## gauss\energy
} else if (kernelTunta == 3) {
  delta0 <-medianE   ## quadexp
}

if (kernelTunta == 1) {
  kernal <- "gauss" ##  gauss
} else if (kernelTunta == 2) {
  kernal <- "energy"  ## energy
}else if (kernelTunta == 3) {
  kernal <- "quadexp"  
}

kernal
delta0 

###### 统计量的权重 ###########
Tm <- trainingSize-lag
TypeI <- 0.05
N<-(1+monitoringPeriod)*trainingSize
k<-1:(N-trainingSize)
## 权重函数
vaha<-(Tm+k)^2/(Tm*Qk(k/Tm,gamma)^2)



#######  Bootstrap ################## 

KT <- max(5, sqrt(log10(trainingSize)))  ## select m^
sqrtLog <- sqrt(log10(trainingSize)/trainingSize)

## flat-top lag-window
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}


Tn <- matrix(0, nrow=nsim, ncol=1) 
## Calculate the statistics for each generated data
Tn.star <- matrix(0, nrow=nsim, ncol=1) 
## Calculate statistics for each Bootsrap data
rej0 <- rep(0,1) ## reject rate

startTime <- Sys.time()
for (s in 1:nsim) {
  yt <- dataGeneration(monitoringPeriod,trainingSize,model = DGP,dataDimension = d)
  Block <- sapply(1:d, function(i) {
    ut <- yt[1:trainingSize, i]
    R <- as.vector(acf(ut, lag.max=2*KT, plot=FALSE)[[1]])
    tmp <- which(abs(R[1:KT]) < 2*sqrt(log10(trainingSize)/trainingSize))
    M <- if(length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
    ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1])
    Ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1]*abs(-M:M))
    max(1, round((2*Ghat^2/(2*ghat^2))^(1/3)*trainingSize^(1/3)))
  })
  #Block
  blocksize <- round( mean(Block) )
  # blocksize
  #生成Bootstrap样本
  y.star <- matrix(nrow = N, ncol = d)
  current_row <- 1
  while(current_row <= N) {
    idx.start <- sample.int(trainingSize, 1)
    idx.len <- rgeom(1, 1/blocksize) + 1
    idx <- (idx.start + 0:(idx.len-1) - 1) %% trainingSize + 1
    needed <- min(idx.len, N - current_row + 1)
    y.star[current_row:(current_row+needed-1), ] <- yt[idx[1:needed], ]
    current_row <- current_row + needed
  }
  #计算原始数据x和Bootstrap样本x.star的各个参数a对应的统计量，并且存储在Tn和Tn.star矩阵中
  Tn[s] <- max(statL(yt,trainingSize,m=lag,alphaE,vaha,kernal,delta0))  ## gauss quadexp energy
  Tn.star[s] <- max(statL(y.star,trainingSize,m=lag,alphaE,vaha,kernal,delta0))
  #每50次循环输出一次拒绝原假设的比例
  if (s %% 50 == 0) {
    critvals <-  quantile(Tn.star[1:s,], 1-TypeI)
    rej0<- sum(critvals < Tn[1:s])
    
    cat('\n')
    rej <- round(100*rej0/s,1)
    print(rej)
  }
  elapsed = difftime(Sys.time(),startTime,units='mins')
  remaining = round((nsim-s)*elapsed/s)
  cat('\r', floor(s/nsim*100),'% complete (',remaining,' min remaining)    ',sep='')
  flush.console()
}

cat('\n')
kernal
delta0
trainingSize
monitoringPeriod

cat('\n')
# rej0
rej
