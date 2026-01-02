rm(list = ls())
setwd("D:/Desktop/paperRev/JointCUSUM/Mul/") # 请设置工作路径

library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(progress) # 进度条包

# 加载源码
sourceCpp('statLMul.cpp')
sourceCpp('eucLMulNew.cpp')
source("Gendata4PaperNew.R")

# --- 1. 参数设置 ---

# 边界函数
Qk <- function(s0, gamma0){
    (1 + s0) * (s0 / (1 + s0))^gamma0
}

# 基础参数
DGP <- 1           # 1 = H0 (Type I Error), 101 = H1 (Power)
lags <- c(0, 1, 2) # 需要组合的 Lag
nlag <- length(lags)
tunta <- 1         # 1=Gauss Kernel
delta0 <- 1
gamma <- 0
trainSize <- 100   # T
monitoringPeriod <- 1 
TypeI <- 0.05      # α
alphaE <- 1

nsim <- 500        # 模拟次数 (外部循环)
M <- 100           # Bootstrap 次数 (内部循环)

# --- 2. 预计算 (Pre-computation) ---

N <- (1 + monitoringPeriod) * trainSize
kernal <- switch(tunta, "1"="gauss", "2"="energy", "3"="quadexp")

# 预计算权重向量 (存入列表，避免重复计算)
weight_list <- list()
for(lag_val in lags) {
    Tm_lag <- trainSize - lag_val
    k_vec <- 1:(N - trainSize)
    w <- (Tm_lag + k_vec)^2 / (Tm_lag * Qk(k_vec/Tm_lag, gamma)^2)
    weight_list[[as.character(lag_val)]] <- w
}

# Bootstrap 相关参数
KT <- max(5, sqrt(log10(trainSize)))  
sqrt_log_T <- 2 * sqrt(log10(trainSize)/trainSize) 

# 编译函数以加速
lambda <- function(t) {
    abst <- abs(t)
    (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}
lambda_compiled <- cmpfun(lambda)
statL_compiled <- cmpfun(statL) 

# Fisher 统计量计算函数
calc_fisher <- function(pvals) {
    # 加上微小量 1e-10 防止 log(0)
    -2 * sum(log(pmax(pvals, 1e-10)))
}

# --- 3. 结果存储容器 ---
rej_individual <- matrix(0, nrow = nsim, ncol = nlag) # 单独 Lag 的拒绝情况
rej_combined   <- rep(0, nsim)                        # 组合检验的拒绝情况
colnames(rej_individual) <- paste0("Lag_", lags)

# --- 4. 主模拟循环 (Sequential) ---

cat("开始模拟 (非并行模式)...\n")
pb <- progress_bar$new(
    format = "[:bar] :percent | 耗时: :elapsed | ETA: :eta",
    total = nsim, clear = FALSE, width = 60
)

for (s in 1:nsim) {
    
    # [Step 1] 生成数据与原始统计量
    # 假设数据维度 d=2 (根据模型调整)
    yt <- gendataPaper(monitoringPeriod, trainSize, DGP)
    d <- ncol(yt)
    
    # 计算原始数据的统计量 T_orig (向量，长度为 nlag)
    T_orig <- numeric(nlag)
    for (j in 1:nlag) {
        lag_val <- lags[j]
        # statL 返回整个序列，取最大值
        seq_val <- statL_compiled(yt, trainSize, m=lag_val, alphaE, 
                                  weight_list[[as.character(lag_val)]], kernal, delta0)
        T_orig[j] <- max(seq_val)
    }
    
    # [Step 2] 确定 Bootstrap Block Size
    # 使用第一个维度计算即可
    ut <- yt[1:trainSize, 1]
    acf_out <- acf(ut, lag.max=2*KT, plot=FALSE, demean=TRUE) 
    R <- drop(acf_out$acf)
    tmp <- which(abs(R[1:KT]) < sqrt_log_T)
    M_lag <- if (length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
    M_range <- -M_lag:M_lag
    ghat <- sum(lambda_compiled(M_range/M_lag) * R[abs(M_range) + 1])
    Ghat <- sum(lambda_compiled(M_range/M_lag) * R[abs(M_range) + 1] * abs(M_range))
    D.SB <- 2 * ghat^2
    blocksize <- if(D.SB <= 1e-8) 1 else max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainSize^(1/3)))
    
    # [Step 3] Bootstrap 循环
    # 存储 M 个 Bootstrap 样本在所有 Lag 下的统计量
    # T_boot 是一个 M x nlag 的矩阵
    T_boot <- matrix(NA, nrow = M, ncol = nlag)
    
    for (b in 1:M) {
        # 生成平稳 Bootstrap 索引
        expected_blocks <- ceiling(N / blocksize) * 2
        len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
        while(sum(len_vec) < N) len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
        # 截断
        len_vec <- len_vec[1:which(cumsum(len_vec) >= N)[1]]
        
        # 拼接索引
        idx_vec <- numeric(0)
        start_points <- sample.int(trainSize, length(len_vec), replace = TRUE)
        for(k in 1:length(len_vec)) {
            idx_vec <- c(idx_vec, (start_points[k] + 0:(len_vec[k]-1) - 1) %% trainSize + 1)
        }
        y_star <- yt[idx_vec[1:N], , drop=FALSE]
        
        # 计算该 Bootstrap 样本的所有 Lag 统计量
        for (j in 1:nlag) {
            lag_val <- lags[j]
            seq_star <- statL_compiled(y_star, trainSize, m=lag_val, alphaE, 
                                       weight_list[[as.character(lag_val)]], kernal, delta0)
            T_boot[b, j] <- max(seq_star)
        }
    }
    
    # [Step 4] 计算组合 P 值 (核心逻辑)
    
    # 为了计算组合统计量的分布，我们需要每一个 Bootstrap 样本的 P 值向量
    # 将原始数据放在第一行，Bootstrap 数据放在后面 -> (M+1) x nlag
    T_total <- rbind(T_orig, T_boot)
    
    # 计算 P 值矩阵: p_mat[i, j] 表示第 i 个样本在第 j 个 Lag 下的 P 值
    # 这里的 P 值是相对于 Null 分布 (T_boot) 的位置
    p_mat <- matrix(0, nrow = M + 1, ncol = nlag)
    
    for (j in 1:nlag) {
        # 当前 Lag 的 Null 分布
        null_dist <- T_boot[, j]
        
        # 对 T_total 中的每一个值，计算它在 null_dist 中的 P 值
        # 向量化计算加速：
        # rank 方法：P = 1 - (Rank - 1)/(M+1) roughly
        # 这里使用严谨定义: (sum(T* >= T) + 1) / (M + 1)
        
        # 技巧：将当前列与 Null 分布一起排序来快速获得 P 值
        # 或者直接用 sapply (因为 M 只有 500，sapply 还是很快的)
        p_mat[, j] <- sapply(T_total[, j], function(t_val) {
            (sum(null_dist >= t_val) + 1) / (M + 1)
        })
    }
    
    # [Step 5] Fisher 组合
    # 对矩阵的每一行计算 Fisher 统计量 -> 得到 W 的分布向量
    W_vec <- apply(p_mat, 1, calc_fisher)
    
    W_orig <- W_vec[1]       # 原始数据的组合统计量
    W_boot <- W_vec[2:(M+1)] # Bootstrap 的组合统计量 (Null 分布)
    
    # [Step 6] 最终决策
    
    # 1. 组合检验决策
    # 计算 W_orig 的 P 值
    p_combined <- (sum(W_boot >= W_orig) + 1) / (M + 1)
    rej_combined[s] <- as.numeric(p_combined < TypeI)
    
    # 2. 单独 Lag 检验决策
    # 直接看 p_mat 的第一行 (原始数据的 P 值)
    p_individual <- p_mat[1, ]
    rej_individual[s, ] <- as.numeric(p_individual < TypeI)
    
    pb$tick()
}

# --- 5. 输出结果 ---

cat("\n\n================ 结果摘要 ================\n")
cat(sprintf("DGP模型: %d | T: %d | nsim: %d | M: %d\n", DGP, trainSize, nsim, M))
cat("------------------------------------------\n")

# 计算拒绝率
rate_combined <- mean(rej_combined) * 100
rate_indiv <- colMeans(rej_individual) * 100

res_df <- data.frame(
    Test_Type = c(colnames(rej_individual), "Combined"),
    Rejection_Rate_Percent = c(rate_indiv, rate_combined)
)

print(res_df)
cat("==========================================\n")