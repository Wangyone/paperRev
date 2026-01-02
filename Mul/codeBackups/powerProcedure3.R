rm(list = ls())
library(progress)
library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)

############ Type I Error or Power ################
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

# 加载外部文件
Rcpp::sourceCpp('statLMul.cpp')
Rcpp::sourceCpp("eucLMulNew.cpp")
source('dataGenerationNew.R')

## boundary function 
Qk <- function(s0,gamma0){
  return(1 + s0)*(s0/(1 + s0))^gamma0
}

############################  Model para
d <- 2 # dimension of data
lag <- 1
kernelTunta <- 3  
## 1 for gauss  && 2 for energy && 3 for quadexp

####################  data length 
trainingSize <- 200
monitoringPeriod <- 1

## para for monitomg statistics
gamma <- 0
deltagauss <- 1*d
nsim <- 1000

DGP <- 1

# 预生成一次数据以初始化参数
yt <- dataGeneration(monitoringPeriod, trainingSize,
                     cptModel = DGP, dataDimension = d)

# 计算 delta0
E <- eucL(yt[1:trainingSize,], lag)
medianE <- (median(E[E > 0]))
rm(E)

alphaE <- 1  ## energy
if (kernelTunta == 1 || kernelTunta == 2) {
  delta0 <- deltagauss  ## gauss\energy
} else if (kernelTunta == 3) {
  delta0 <- medianE     ## quadexp
}

if (kernelTunta == 1) {
  kernal <- "gauss"   ## gauss
} else if (kernelTunta == 2) {
  kernal <- "energy"  ## energy
} else if (kernelTunta == 3) {
  kernal <- "quadexp"  
}

# 打印参数确认
print(kernal)
print(delta0) 

###### 统计量的权重 ###########
Tm <- trainingSize - lag
TypeI <- 0.05
N <- (1 + monitoringPeriod) * trainingSize
k <- 1:(N - trainingSize)
## 权重函数
vaha <- (Tm + k)^2 / (Tm * Qk(k/Tm, gamma)^2)


#######  Bootstrap 参数 ################## 
KT <- max(5, sqrt(log10(trainingSize)))  ## select m^
sqrtLog <- sqrt(log10(trainingSize)/trainingSize)

## flat-top lag-window
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}

# 编译关键函数以加速
statL_compiled <- cmpfun(statL)
lambda_compiled <- cmpfun(lambda)

# 预计算常量
sqrt_log_ts <- 2 * sqrt(log10(trainingSize)/trainingSize)
pb <- progress_bar$new(total = nsim, format = "[:bar] :percent (:eta remaining)")

# 初始化结果容器
Tn <- matrix(0, nrow=nsim, ncol=1) 
Tn.star <- matrix(0, nrow=nsim, ncol=1) 
rej0 <- rep(0,1) 

# --- 主循环 ---
for (s in 1:nsim) {
  # 1. 生成数据
  yt <- dataGeneration(monitoringPeriod, trainingSize, cptModel = DGP, dataDimension = d)
  
  # 2. 自动选择 Block Size (Politis & White/Romano method)
  Block <- sapply(1:d, function(i) {
    ut <- yt[1:trainingSize, i]
    R <- drop(acf(ut, lag.max=2*KT, plot=FALSE, demean=FALSE)$acf)
    tmp <- which.max(abs(R[1:KT]) >= sqrt_log_ts)
    M <- if(tmp > 1) 2*(tmp - 1) else 2*(KT - 1)
    M_vals <- -M:M
    R_vals <- R[abs(M_vals) + 1]
    lambda_vals <- lambda_compiled(M_vals/M)
    ghat <- sum(lambda_vals * R_vals)
    Ghat <- sum(lambda_vals * R_vals * abs(M_vals))
    max(1, round((Ghat^2/ghat^2)^(1/3) * trainingSize^(1/3)))
  })
  
  blocksize <- round(mean(Block))
  # 避免 blocksize 为 0 或过小的极端情况
  if(blocksize < 1) blocksize <- 1 
  
  # 3. Stationary Bootstrap 数据重采样 (包含修复逻辑)
  
  # (3.1) 生成块长 (Lengths)
  # 初始估计需要的块数 (乘以2作为余量)
  num_blocks_est <- ceiling(N / blocksize) * 2
  lengths <- rgeom(num_blocks_est, 1/blocksize) + 1
  
  # [修复 BUG 的关键步骤]：检查总长度是否足够，不够则循环追加
  while (sum(lengths) < N) {
    # 计算还缺多少长度，并估算需要补充的块数 (+5为了保险)
    needed <- N - sum(lengths)
    add_blocks <- ceiling(needed / blocksize) + 5
    lengths <- c(lengths, rgeom(add_blocks, 1/blocksize) + 1)
  }
  
  # (3.2) 根据最终确定的块数量，生成起始点 (Starts)
  # 注意：Starts 的数量必须和 Lengths 的数量完全一致
  total_blocks_final <- length(lengths)
  starts <- sample.int(trainingSize, total_blocks_final, replace=TRUE)
  
  # (3.3) 构建索引列表
  idx_list <- lapply(1:total_blocks_final, function(j) {
    # 环形索引逻辑 (Circular Bootstrap)
    (starts[j] + 0:(lengths[j]-1) - 1) %% trainingSize + 1
  })
  
  # (3.4) 拼接并截取前 N 个样本
  # 因为上面保证了 sum(lengths) >= N，这里绝对不会越界
  y.star <- do.call(rbind, lapply(idx_list, function(idx) {
    yt[idx, , drop=FALSE]
  }))[1:N, ]
  
  # 4. 计算统计量
  Tn[s] <- max(statL_compiled(yt, trainingSize, m=lag, alphaE, vaha, kernal, delta0))
  Tn.star[s] <- max(statL_compiled(y.star, trainingSize, m=lag, alphaE, vaha, kernal, delta0))
  
  # 5. 定期输出结果
  if (s %% 100 == 0) {
    critvals <- quantile(Tn.star[1:s], 1-TypeI, na.rm=TRUE)
    rej <- round(100 * mean(Tn[1:s] > critvals, na.rm=TRUE), 1)
    cat('\nRejection rate:', rej, '%\n')
  }
  
  pb$tick()
  # 手动GC防止内存泄漏
  if (s %% 100 == 0) gc()
}

cat('\nFinal Results:\n')
print(paste("Kernel:", kernal))
print(paste("Delta0:", delta0))
print(paste("Training Size:", trainingSize))
print(paste("Monitoring Period:", monitoringPeriod))

critvals_final <- quantile(Tn.star, 1-TypeI, na.rm=TRUE)
rej_final <- mean(Tn > critvals_final, na.rm=TRUE)
print(paste("Final Rejection Rate:", round(100 * rej_final, 2), "%"))