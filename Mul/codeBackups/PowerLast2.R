rm(list = ls())
######################### Type1 error #########################
setwd('D:/Desktop/JointEcf/Mul')
library(ggplot2)
library(corpcor)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
sourceCpp('statLMul.cpp')
sourceCpp('eucLMul.cpp')
source('gendataMulLast.R')
## weighting function 
Qk<-function(s0,gamma0){
  u=(1+s0)*(s0/(1+s0))^gamma0
}
d <- 2
DGP <- 206
m0 <- 0
tunta <- 1  ## 1 for gauss  && 2 for energy && 3 for quadexp
#################### para
T <- 200
L <- 1


deltagauss <- 1
gamma<-0


yt <- genDataMul(L,T,model = DGP)
MTSplot(yt)

var(yt)

E <- eucL(yt[1:T,],T,m0) # eucL(yt[1:T,],T,m0)
## 1 for gauss  && 2 for energy && 3 for quadexp
alphaE <- 1  ##energy

if (tunta == 1||tunta == 2) {
  delta0 <- deltagauss  ## gauss\energy
} else if (tunta == 3) {
  delta0 <- (median(E[E > 0]))  ## quadexp
}


if (tunta == 1) {
  kernal <- "gauss" ##  gauss
} else if (tunta == 2) {
  kernal <- "energy"  ## energy
}else if (tunta == 3) {
  kernal <- "quadexp"  ## quadexp
}

kernal
delta0 

###### 统计量的权重 ###########
Tm <- T-m0
TypeI <- 0.05
N<-(1+L)*T
k<-1:(N-T)
## 权重函数

vaha<-(Tm+k)^2/(Tm*Qk(k/Tm,gamma)^2)



#######  Bootstrap ################## 

KT <- max(5, sqrt(log10(T)))  ## select m^
## flat-top lag-window
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}

nsim<-1000
Tn <- matrix(0, nrow=nsim, ncol=1) ## Calculate the statistics for each generated data
Tn.star <- matrix(0, nrow=nsim, ncol=1) ## Calculate statistics for each Bootsrap data
rej0 <- rep(0,1) ## reject rate

startTime <- Sys.time()
for (s in 1:nsim) {
  
  yt <- genDataMul(L,T,model = DGP)
  Block <- numeric(d)
  for (i in 1:d) {
    ut <- yt[1:T,i]
    R <- as.vector(acf(ut, lag.max=2*KT, plot=FALSE)[[1]])
    #使用acf函数计算时间序列数据x的自相关系数。x[1:T]表示将x截取为长度为T的子序列，
    # lag.max=2*KT表示计算自相关系数的最大滞后阶数
    tmp <- which(abs(R[1:KT]) < 2*sqrt(log10(T)/T)) ##满足条件的m^+k,k=1,...,KT，即寻找m^
    if (length(tmp) > 0) {
      M <- 2*(tmp[1] - 1)
    } else {
      M <- 2*(KT - 1)
    }
    ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1])
    Ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1]*abs(-M:M))
    D.SB <- 2*ghat^2
    blocksize <- max(1, round((2*Ghat^2/D.SB)^(1/3)*T^(1/3)))
  Block[i] <- blocksize
  }
  #Block
  blocksize <- round( mean(Block) )
  # blocksize
  #生成Bootstrap样本
  y.star <- matrix(rep(NA,1*2),nrow=1,ncol=d)

  while (nrow(y.star) < (N+T/2)) {
    idx.start <- sample.int(T, 1, replace=TRUE) #生成{1,...,T}的离散均匀分布值，作为随机块的起点
    idx.len <- rgeom(1, 1/blocksize) + 1       # 
    idx <- idx.start:(idx.start + idx.len)
    idx <- (idx - 1) %% T + 1
    idx
    y.star <- rbind(y.star, yt[(idx),])
  }
  y.star <- y.star[-1,]
  y.star <- y.star[1:N,]
  
  #MTSplot(yt)
  #计算原始数据x和Bootstrap样本x.star的各个参数a对应的统计量，并且存储在Tn和Tn.star矩阵中
  Tn[s] <- max(statL(yt,T,m=m0,alphaE,vaha,kernal,delta0))  ## gauss quadexp energy
  Tn.star[s] <- max(statL(y.star,T,m=m0,alphaE,vaha,kernal,delta0))
  
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
#cat('\n')
rej0
rej

kernal
delta0



