rm(list = ls())
### 多参数比较：Type I Error or Power ###

# --- 1. 环境设置与依赖加载 ---
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(compiler) 
library(progress) 

# 加载外部文件
sourceCpp('statLMul.cpp')
sourceCpp('eucLMulNew.cpp')
source("Gendata4PaperNew.R") # 确保这是刚才修复过的版本

# --- 2. 辅助函数定义与编译 ---

## Boundary function 
Qk <- function(s0, gamma0){
  (1 + s0) * (s0 / (1 + s0))^gamma0
}

## Flat-top lag-window
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1))
}

# 编译 R 函数以提高速度
lambda_compiled <- cmpfun(lambda)
Qk_compiled <- cmpfun(Qk)

# --- 3. 基础参数设置 ---

DGP <- 1         # Model ID
lag <- 0         # m0
kernelTunta <- 1 # 1=gauss, 2=energy, 3=quadexp

trainingSize <- 100
monitoringPeriod <- 1
gamma <- 0
nsim <- 1000     # 模拟次数
TypeI <- 0.05

# --- 4. 关键修复：预先确定维度 d ---
# 在设置 deltagauss 之前，必须先知道 d 是多少
# 生成一个临时样本来探测维度
cat("正在初始化参数与探测维度...\n")
temp_yt <- gendataPaper(monitoringPeriod, trainingSize, model = DGP)
d <- ncol(temp_yt)
cat(sprintf("检测到数据维度 d = %d\n", d))

# 定义 delta 参数组 (依赖于 d)
deltagauss <- c(0.5, 1, 2, 4, 6, 8, 10) * d
simulationResult <- numeric(length(deltagauss))
names(simulationResult) <- paste0("Delta_", deltagauss)

# 预计算常量
Tm <- trainingSize - lag
N <- (1 + monitoringPeriod) * trainingSize
k <- 1:(N - trainingSize)
vaha <- (Tm + k)^2 / (Tm * Qk_compiled(k/Tm, gamma)^2)
alphaE <- 1
KT <- max(5, sqrt(log10(trainingSize)))
log_T_limit <- 2 * sqrt(log10(trainingSize) / trainingSize)

# 确定 Kernel 名称
if (kernelTunta == 1) {
  kernal <- "gauss"
} else if (kernelTunta == 2) {
  kernal <- "energy"
} else if (kernelTunta == 3) {
  kernal <- "quadexp"
}

# 如果是 quadexp，通常 delta 是固定的中位数
if (kernelTunta == 3) {
  E <- eucL2(temp_yt[1:trainingSize, ], trainingSize, lag)
  delta0_fixed <- median(E[E > 0])
}

# --- 5. 模拟主循环 ---

cat("开始模拟...\n")
cat("Kernel:", kernal, "| T:", trainingSize, "| Sim:", nsim, "\n")

# 外层循环：遍历参数 delta
for (var in 1:length(deltagauss)) {
  
  # 确定当前的 delta0
  if (kernelTunta == 1 || kernelTunta == 2) {
    curr_delta0 <- deltagauss[var]
  } else {
    curr_delta0 <- delta0_fixed
    # 如果是 quadexp，外层循环其实只跑一次或只跑不同的 delta0_fixed 倍数，
    # 这里假设如果是 quadexp，我们忽略 deltagauss 的变化，或者您可以根据需要修改逻辑
  }
  
  Tn <- numeric(nsim)
  Tn.star <- numeric(nsim)
  
  # 初始化进度条
  pb <- progress_bar$new(
    format = paste0(" Para Delta=", curr_delta0, " [:bar] :percent eta: :eta"),
    total = nsim, clear = FALSE, width = 60
  )
  
  # 内层循环：Monte Carlo 模拟
  for (s in 1:nsim) {
    
    # 1. 生成数据
    yt <- gendataPaper(monitoringPeriod, trainingSize, model = DGP)
    # 确保 d 在循环内也是正确的 (虽然上面已经获取)
    d <- ncol(yt) 
    
    # 2. 选择 Block Size (向量化优化)
    Block_sizes <- sapply(1:d, function(i) {
      ut <- yt[1:trainingSize, i]
      # 修复: lag.max 必须为整数
      acf_res <- acf(ut, lag.max = as.integer(2*KT), plot = FALSE)$acf
      R <- drop(acf_res) # 降维处理
      
      tmp <- which(abs(R[1:as.integer(KT)]) < log_T_limit)
      M <- if (length(tmp) > 0) 2 * (tmp[1] - 1) else 2 * (KT - 1)
      
      idx_range <- -M:M
      # R下标从1开始，所以 +1
      vals <- R[abs(idx_range) + 1]
      weights <- lambda_compiled(idx_range / M)
      
      ghat <- sum(weights * vals)
      Ghat <- sum(weights * vals * abs(idx_range))
      D.SB <- 2 * ghat^2
      
      if (is.na(D.SB) || D.SB < 1e-10) return(1)
      max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainingSize^(1/3)))
    })
    
    blocksize <- max(1, round(mean(Block_sizes)))
    
    # 3. 生成 Stationary Bootstrap 样本 (极速优化版)
    expected_blocks <- ceiling((N + trainingSize) / blocksize) + 5
    len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
    
    while(sum(len_vec) < N) {
      len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
    }
    
    cumsum_len <- cumsum(len_vec)
    cutoff_idx <- which(cumsum_len >= N)[1]
    len_vec <- len_vec[1:cutoff_idx]
    num_blocks <- length(len_vec)
    
    start_vec <- sample.int(trainingSize, num_blocks, replace = TRUE)
    
    idx_list <- lapply(1:num_blocks, function(j) {
      (start_vec[j] + 0:(len_vec[j]-1) - 1) %% trainingSize + 1
    })
    
    bootstrap_idx <- unlist(idx_list)[1:N]
    
    # 修复: 保持矩阵结构 drop=FALSE
    y.star <- yt[bootstrap_idx, , drop = FALSE]
    
    # 4. 计算统计量
    Tn[s] <- max(statL(yt, trainingSize, m = lag, alphaE, vaha, kernal, curr_delta0))
    Tn.star[s] <- max(statL(y.star, trainingSize, m = lag, alphaE, vaha, kernal, curr_delta0))
    
    pb$tick()
  }
  
  # 计算拒绝率
  critvals <- quantile(Tn.star, 1 - TypeI, na.rm = TRUE)
  rej_rate <- mean(Tn > critvals, na.rm = TRUE)
  
  simulationResult[var] <- round(rej_rate * 100, 2)
  
  cat(sprintf("\n[Result] Delta: %.2f -> Rej Rate: %.2f%%\n", curr_delta0, simulationResult[var]))
}

cat("\n==============================\n")
cat("最终多参数模拟结果:\n")
print(simulationResult)