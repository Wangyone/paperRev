rm(list = ls())
library(progress)
library(compiler)
library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)

# 设置工作目录
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

# 加载C++代码和R函数
Rcpp::sourceCpp('statLMul.cpp')
Rcpp::sourceCpp("eucLMulNew.cpp")
source('dataGenerationSupply.R')

########################### 1. 辅助函数与参数设置 ###########################

# 边界权重函数
Qk <- function(s0, gamma0) {
    return(1 + s0) * (s0 / (1 + s0))^gamma0
}

# Fisher 组合统计量计算函数
# W = -2 * sum(ln(p))
calc_fisher_stat <- function(pvals) {
    # 避免 log(0)
    pvals <- pmax(pvals, 1e-10) 
    return(-2 * rowSums(log(pvals)))
}

# 参数配置
d <- 2  # 维度
DGP <- 1  # 数据生成模型
lags <- c(0, 1, 2)  # 组合的 Lags
nlag <- length(lags)
kernelTunta <- 3  # 3 = quadexp

trainingSize <- 200
monitoringPeriod <- 1
gamma <- 0
nsim <- 1000  # Bootstrap 次数 (Warp-Speed)
TypeI <- 0.05
alphaE <- 1 # energy参数

# 核函数配置
if (kernelTunta == 1) kernal <- "gauss" else if (kernelTunta == 2) kernal <- "energy" else kernal <- "quadexp"

# 预计算 Delta0 (基于 Lag=0)
yt_init <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
if (kernelTunta == 3) {
    E_init <- eucL(yt_init[1:trainingSize,], 0)
    delta0 <- median(E_init[E_init > 0])
} else {
    delta0 <- 1 * d 
}

########################### 2. 预计算不同 Lag 的权重 ###########################

N <- (1 + monitoringPeriod) * trainingSize
weight_list <- list()

cat("预计算权重向量...\n")
for (lag_val in lags) {
    Tm <- trainingSize - lag_val
    k_vec <- 1:(N - trainingSize)
    # 权重公式 vaha
    vaha_vec <- (Tm + k_vec)^2 / (Tm * Qk(k_vec/Tm, gamma)^2)
    weight_list[[as.character(lag_val)]] <- vaha_vec
}

########################### 3. 预编译与初始化 ###########################

lambda <- function(t) {
    abst <- abs(t)
    (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1))
}
statL_compiled <- cmpfun(statL)
lambda_compiled <- cmpfun(lambda)

KT <- max(5, sqrt(log10(trainingSize)))
sqrt_log_ts <- 2 * sqrt(log10(trainingSize) / trainingSize)

# 存储矩阵
Tn <- matrix(0, nrow = nsim, ncol = nlag)      # 原始统计量
Tn_star <- matrix(0, nrow = nsim, ncol = nlag) # Bootstrap 统计量
colnames(Tn) <- paste0("lag", lags)
colnames(Tn_star) <- paste0("lag", lags)

pb <- progress_bar$new(total = nsim, format = "[:bar] :percent 耗时: :elapsedfull")

########################### 4. Warp-Speed 模拟循环 ###########################

cat("\n开始 Warp-Speed 模拟 (包含中间过程监测)...\n")

for (s in 1:nsim) {
    
    # -------------------------------------------------------
    # A. 生成数据与计算原始统计量 (Tn)
    # -------------------------------------------------------
    yt <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
    
    for (i in 1:nlag) {
        m <- lags[i]
        w_vec <- weight_list[[as.character(m)]]
        # 计算序列并取最大值作为该 lag 的统计量
        stat_seq <- statL_compiled(yt, trainingSize, m, alphaE, w_vec, kernal, delta0)
        Tn[s, i] <- max(stat_seq, na.rm = TRUE)
    }
    
    # -------------------------------------------------------
    # B. Bootstrap 采样
    # -------------------------------------------------------
    # 使用 Lag 0 或平均 Block size 进行采样
    Block <- sapply(1:d, function(i) {
        ut <- yt[1:trainingSize, i]
        R <- drop(acf(ut, lag.max = 2 * KT, plot = FALSE, demean = FALSE)$acf)
        tmp <- which.max(abs(R[1:KT]) >= sqrt_log_ts)
        M <- if (tmp > 1) 2 * (tmp - 1) else 2 * (KT - 1)
        M_vals <- -M:M
        lambda_vals <- lambda_compiled(M_vals / M)
        ghat <- sum(lambda_vals * R[abs(M_vals) + 1])
        Ghat <- sum(lambda_vals * R[abs(M_vals) + 1] * abs(M_vals))
        max(1, round((Ghat^2 / ghat^2)^(1/3) * trainingSize^(1/3)))
    })
    blocksize <- round(mean(Block))
    
    total_blocks <- ceiling(N / blocksize) * 2
    starts <- sample.int(trainingSize, total_blocks, replace = TRUE)
    lengths <- rgeom(total_blocks, 1 / blocksize) + 1
    idx_list <- lapply(1:total_blocks, function(j) {
        (starts[j] + 0:(lengths[j] - 1) - 1) %% trainingSize + 1
    })
    y_star <- do.call(rbind, lapply(idx_list, function(idx) {
        yt[idx, , drop = FALSE]
    }))[1:N, ]
    
    # -------------------------------------------------------
    # C. 计算 Bootstrap 统计量 (Tn_star)
    # -------------------------------------------------------
    for (i in 1:nlag) {
        m <- lags[i]
        w_vec <- weight_list[[as.character(m)]]
        stat_seq_star <- statL_compiled(y_star, trainingSize, m, alphaE, w_vec, kernal, delta0)
        Tn_star[s, i] <- max(stat_seq_star, na.rm = TRUE)
    }
    
    # -------------------------------------------------------
    # D. 中间过程监测 (每 100 次)
    # -------------------------------------------------------
    if (s %% 100 == 0) {
        # 提取当前的样本池
        curr_Tn <- Tn[1:s, , drop=FALSE]
        curr_Tn_star <- Tn_star[1:s, , drop=FALSE]
        
        # 计算当前的 P 值 (基于当前 s 个 bootstrap 样本)
        # 注意：Warp-Speed 原理是用 Tn_star 的分布来评估 Tn
        # 这里为了监测，我们批量计算前 s 个样本的 P 值
        
        # 1. 计算 P_obs (Tn 在 Tn_star 分布中的位置)
        curr_P_obs <- matrix(0, nrow = s, ncol = nlag)
        for(j in 1:nlag) {
            # 使用 rank 计算位置比循环更快
            # 统计 Tn_star >= Tn[k] 的个数。
            # 技巧: 将 Tn 和 Tn_star 放在一起排序
            # 但为了简单展示，这里用 sapply (仅在监测步运行，不影响整体速度)
            curr_P_obs[, j] <- sapply(curr_Tn[, j], function(x) (sum(curr_Tn_star[, j] >= x) + 1) / (s + 1))
        }
        
        # 2. 计算 P_star (Bootstrap 自身的 P 值，用于确定临界值)
        curr_P_star <- matrix(0, nrow = s, ncol = nlag)
        for(j in 1:nlag) {
            # Bootstrap 样本的 P 值就是其归一化秩
            # 统计量越大，P 值越小
            curr_P_star[, j] <- rank(-curr_Tn_star[, j]) / (s + 1)
        }
        
        # 3. Fisher 组合
        curr_W_obs <- calc_fisher_stat(curr_P_obs)
        curr_W_star <- calc_fisher_stat(curr_P_star)
        
        # 4. 确定阈值与拒绝率
        crit_val <- quantile(curr_W_star, 1 - TypeI)
        curr_rej_rate <- mean(curr_W_obs > crit_val)
        
        cat(sprintf("\n[Monitor] Iter: %d | Combined Rej Rate: %.1f%%\n", s, curr_rej_rate * 100))
    }
    
    pb$tick()
    if (s %% 200 == 0) gc()
}

########################### 5. 最终结果计算 ###########################

cat("\n正在计算最终 P 值与拒绝率...\n")

# 最终使用完整的 nsim 个样本进行计算
final_P_obs <- matrix(0, nrow = nsim, ncol = nlag)
final_P_star <- matrix(0, nrow = nsim, ncol = nlag)

# 向量化计算最终 P 值
for (i in 1:nlag) {
    # Bootstrap 样本自身的 P 值 (用于 Null 分布)
    final_P_star[, i] <- rank(-Tn_star[, i]) / (nsim + 1)
    
    # 原始样本的 P 值 (相对于 Bootstrap 分布)
    # 利用 findInterval 加速比较
    # sort_star <- sort(Tn_star[, i])
    # rank_idx <- findInterval(Tn[, i], sort_star) 
    # count_ge <- nsim - rank_idx
    # 更加精确的方法:
    final_P_obs[, i] <- sapply(Tn[, i], function(x) (sum(Tn_star[, i] >= x) + 1) / (nsim + 1))
}

# Fisher 组合
final_W_obs <- calc_fisher_stat(final_P_obs)
final_W_star <- calc_fisher_stat(final_P_star)

# 阈值与拒绝决策
final_crit_val <- quantile(final_W_star, 1 - TypeI)
is_rejected <- final_W_obs > final_crit_val
final_rej_rate <- mean(is_rejected)

# 单独 Lag 的拒绝率 (用于对比)
indiv_rej_rates <- numeric(nlag)
for(i in 1:nlag) {
    crit <- quantile(Tn_star[, i], 1 - TypeI)
    indiv_rej_rates[i] <- mean(Tn[, i] > crit)
}

cat("\n================ 最终结果 ================\n")
cat(sprintf("Model: DGP%d | Dimension: %d | Kernel: %s\n", DGP, d, kernal))
cat(sprintf("Lags: {%s}\n", paste(lags, collapse=",")))
cat(sprintf("Bootstrap Replications: %d\n", nsim))
cat("------------------------------------------\n")
cat("单 Lag 拒绝率:\n")
for(i in 1:nlag) {
    cat(sprintf("  Lag %d: %.1f%%\n", lags[i], indiv_rej_rates[i] * 100))
}
cat("------------------------------------------\n")
cat(sprintf("组合检验 (Fisher) 拒绝率: %.1f%%\n", final_rej_rate * 100))
cat("==========================================\n")