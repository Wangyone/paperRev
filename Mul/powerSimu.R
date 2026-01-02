rm(list = ls())
######################### Simulation: Type1 I error && Power   
#########################################################

# --- 0. 环境与依赖 ---
setwd("D:/Desktop/paperRev/JointCUSUM/Mul/")

library(compiler)     # 用于字节编译R函数
library(progress)     # 进度条
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)

# 加载 C++ 和 R 源码
sourceCpp('statLMul.cpp')
sourceCpp('eucLMulNew.cpp')
source("Gendata4PaperNew.R")

# --- 1. 参数设置 ---

# Boundary function 
Qk <- function(s0, gamma0){
    (1 + s0) * (s0 / (1 + s0))^gamma0
}

# 核心参数
gamma <- 0

cptModel <- 201     # Model ID
m0 <- 0             # Lag
tunta <- 1          # Kernel: 1=gauss, 2=energy, 3=quadexp
delta0 <- 1         # Kernel parameter
trainSize <- 100    # Training size 
monitoringPeriod <- 1 # Monitoring period factor 

# 核函数名称映射
if (tunta == 1) {
    kernal <- "gauss"
} else if (tunta == 2) {
    kernal <- "energy"
} else if (tunta == 3) {
    kernal <- "quadexp" 
}

# 统计量权重计算
Tm <- trainSize - m0
TypeI <- 0.05
N <- (1 + monitoringPeriod) * trainSize
k_vec <- 1:(N - trainSize)

# 预计算权重向量
vaha <- (Tm + k_vec)^2 / (Tm * Qk(k_vec/Tm, gamma)^2)
alphaE <- 1

# Bootstrap 参数
KT <- max(5, sqrt(log10(trainSize)))  
sqrt_log_T <- 2 * sqrt(log10(trainSize)/trainSize) 

# Flat-top lag-window function
lambda <- function(t) {
    abst <- abs(t)
    (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}

# --- 2. 性能优化：编译关键函数 ---
lambda_compiled <- cmpfun(lambda)
statL_compiled <- cmpfun(statL)

# --- 3. 模拟主循环 ---
nsim <- 1000 # Number of trials

# 初始化存储容器
Tn <- numeric(nsim)
Tn.star <- numeric(nsim)

# 初始化进度条 (简化显示，详细日志通过 cat 输出)
pb <- progress_bar$new(
    format = "  Running [:bar] :percent | ETA: :eta",
    total = nsim, clear = FALSE, width = 60
)

cat("Starting Simulation...\n")
cat(sprintf("Model: %d | Kernel: %s | T: %d | L: %s\n", cptModel, kernal, trainSize, monitoringPeriod))

for (s in 1:nsim) {
    
    # (A) 生成数据
    yt <- gendataPaper(monitoringPeriod = monitoringPeriod, trainSize = trainSize, model = cptModel)
    d <- ncol(yt)
    # (B) 自动选择 Block Size (向量化优化版)
    Block_sizes <- sapply(1:d, function(i) {
        ut <- yt[1:trainSize, i]
        acf_out <- acf(ut, lag.max=2*KT, plot=FALSE, demean=TRUE) 
        R <- drop(acf_out$acf)
        
        tmp <- which(abs(R[1:KT]) < sqrt_log_T)
        M <- if (length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
        
        M_range <- -M:M
        R_vals <- R[abs(M_range) + 1]
        lam_vals <- lambda_compiled(M_range/M)
        
        ghat <- sum(lam_vals * R_vals)
        Ghat <- sum(lam_vals * R_vals * abs(M_range))
        D.SB <- 2 * ghat^2
        
        if(D.SB <= 1e-8) return(1) 
        max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainSize^(1/3)))
    })
    
    blocksize <- max(1, round(mean(Block_sizes)))
    
    # (C) Stationary Bootstrap (向量化索引优化版)
    expected_blocks <- ceiling(N / blocksize) * 2
    len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
    
    while(sum(len_vec) < N) {
        len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
    }
    
    cumsum_len <- cumsum(len_vec)
    cutoff_idx <- which(cumsum_len >= N)[1]
    len_vec <- len_vec[1:cutoff_idx]
    num_blocks <- length(len_vec)
    start_vec <- sample.int(trainSize, num_blocks, replace = TRUE)
    
    idx_list <- lapply(1:num_blocks, function(j) {
        (start_vec[j] + 0:(len_vec[j]-1) - 1) %% trainSize + 1
    })
    
    bootstrap_idx <- unlist(idx_list)[1:N]
    y.star <- yt[bootstrap_idx, , drop=FALSE]
    
    # (D) 计算统计量
    Tn[s] <- max(statL_compiled(yt, trainSize, m=m0, alphaE, vaha, kernal, delta0))
    Tn.star[s] <- max(statL_compiled(y.star, trainSize, m=m0, alphaE, vaha, kernal, delta0))
    
    # (E) 每100次输出结果
    if (s %% 100 == 0) {
        # 计算当前的拒绝率
        curr_crit <- quantile(Tn.star[1:s], 1 - TypeI, na.rm=TRUE)
        curr_rej <- mean(Tn[1:s] > curr_crit, na.rm=TRUE)
        
        # 格式化输出 (使用 \n 换行以避免覆盖进度条)
        cat(sprintf("\n[Iter %4d/%d] Current Rejection Rate: %.2f%%\n", s, nsim, 100 * curr_rej))
    }
    
    pb$tick()
    
    # 内存清理
    if(s %% 200 == 0) gc(verbose = FALSE)
}

# --- 4. 最终结果输出 ---
final_crit <- quantile(Tn.star, 1 - TypeI, na.rm=TRUE)
final_rej <- mean(Tn > final_crit, na.rm=TRUE)

cat("\n==============================\n")
cat("Simulation Finished.\n")
cat(paste("Kernel:", kernal, "\n"))
cat(paste("Delta0:", delta0, "\n"))
cat(paste("Final Rejection Rate:", round(100 * final_rej, 2), "%\n"))
