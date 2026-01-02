rm(list  =  ls())
source("D:/Desktop/JointCUSUM/Mul/dataBlock.R", echo=TRUE)
# parameter
d <- 3
AllData <- 200
# 相关系数
phi <- seq(0,0.9,by=0.1)
#######  Bootstrap ################## 
## flat-top lag-window
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}
## block
blockYt <- function(Xt){
  d <- NCOL(Xt)
  trainingSize <- NROW(Xt)
  KT <- max(5, sqrt(log10(trainingSize)))  ## select m^
  Block <- numeric(d)
  for (i in 1:d) {
    ut <- Xt[,i]
    R <- as.vector(acf(ut, lag.max=2*KT, plot=FALSE)[[1]])
    #使用acf函数计算时间序列的自相关系数
    # lag.max=2*KT表示计算自相关系数的最大滞后阶数
    tmp <- which(abs(R[1:KT]) < 2*sqrt(log10(trainingSize)/trainingSize)) 
    ##满足条件的m^+k,k=1,...,KT，即寻找m^
    if (length(tmp) > 0) {
      M <- 2*(tmp[1] - 1)
    } else {
      M <- 2*(KT - 1)
    }
    ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1])
    Ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1]*abs(-M:M))
    D.SB <- 2*ghat^2
    blocksize <- max(1, round((2*Ghat^2/D.SB)^(1/3)*trainingSize^(1/3)))
    Block[i] <- blocksize
  }
  #Block
  blocksize <- ceiling( mean(Block) )
  return(blocksize)
}

## 储存结果
nsim <- 10
BlockMatrix <- matrix(NA,nrow=nsim,ncol=NROW(phi))



for(i in 1:nsim){
  for (j  in 1:NROW(phi)) {
    beta_y <- phi[j]
    yt <- dataBlock(AllData,beta_y,d)
    BlockMatrix[i,j] <- blockYt(yt)
  }
}



