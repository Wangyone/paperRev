rm(list = ls())
library(progress)
library(compiler)
library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)

setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')
Rcpp::sourceCpp('statLMul.cpp')
Rcpp::sourceCpp("eucLMulNew.cpp")
source("D:/Desktop/paperRev/JointCUSUM/Mul/codeBackups/dataGenerationSupply.R", echo=TRUE)

############### 1. 辅助函数与参数设置 ###############

# 边界权重函数
Qk <- function(s0, gamma0) {
    return(1 + s0) * (s0 / (1 + s0))^gamma0
}

# Fisher 组合统计量计算函数（文章式2.2）
calc_fisher_stat <- function(pvals, weights = rep(1, length(pvals))) {
    # 避免 log(0) - 使用文章方法后p值不会为0
    return(-2 * sum(weights * log(pvals)))
}

# 参数配置
d <- 2  # 维度
DGP <- 1  # 数据生成模型
lags <- c(0, 1, 2)  # 组合的Lags

nlag <- length(lags)
kernelTunta <- 1  

trainingSize <- 200  # 样本大小
monitoringPeriod <- 1
M <- 500
# Bootstrap次数（严格按照文章）

gamma <- 0
nsim <- 500  # 模拟次数
TypeI <- 0.05 # 显著性水平
alphaE <- 1

# 核函数配置
if (kernelTunta == 1) kernal <- "gauss" else if (kernelTunta == 2) kernal <- "energy" else kernal <- "quadexp"

# 预计算kernel_3的参数Delta0 (基于Lag=0)
if (kernelTunta == 3) {
    yt_init <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
    E_init <- eucL(yt_init[1:trainingSize,], 0)
    delta0 <- median(E_init[E_init > 0])
} else {
    delta0 <- 1 * d 
}

############### 2. 预计算权重向量 ##############
N <- (1 + monitoringPeriod) * trainingSize
weight_list <- list()

cat("预计算权重向量...\n")
for (lag_val in lags) {
    Tm <- trainingSize - lag_val
    k_vec <- 1:(N - trainingSize)
    vaha_vec <- (Tm + k_vec)^2 / (Tm * Qk(k_vec/Tm, gamma)^2)
    weight_list[[as.character(lag_val)]] <- vaha_vec
}

############### 3. 预编译与初始化 ##############

lambda <- function(t) {
    abst <- abs(t)
    (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1))
}

statL_compiled <- cmpfun(statL)
lambda_compiled <- cmpfun(lambda)

KT <- max(5, sqrt(log10(trainingSize)))
sqrt_log_ts <- 2 * sqrt(log10(trainingSize) / trainingSize)

# 存储结果
rej_individual <- matrix(0, nrow = nsim, ncol = nlag)  # 单独检验拒绝率
rej_combined <- rep(0, nsim)  # 组合检验拒绝率
colnames(rej_individual) <- paste0("lag", lags)

pb <- progress_bar$new(total = nsim, 
                       format = "[:bar] :percent 耗时: :elapsedfull")

cat("\n开始模拟...\n")

############### 4. 主模拟循环 ##############

for (sim in 1:nsim) {
    
    # -------------------------------------------------------
    # 步骤1: 生成原始数据并计算原始统计量 T_n^[0]
    # -------------------------------------------------------
    yt <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
    
    # 计算原始统计量 T_n^[0]
    T_original <- numeric(nlag)
    for (i in 1:nlag) {
        m <- lags[i]
        w_vec <- weight_list[[as.character(m)]]
        stat_seq <- statL_compiled(yt, trainingSize, m, alphaE, w_vec, kernal, delta0)
        T_original[i] <- max(stat_seq, na.rm = TRUE)
    }
    
    # 确定blocksize（使用原始数据的自相关）
    Block <- sapply(1:d, function(i) {
        ut <- yt[1:trainingSize, i]
        R <- drop(acf(ut, lag.max = 2*KT, plot = FALSE, demean = FALSE)$acf)
        tmp <- which.max(abs(R[1:KT]) >= sqrt_log_ts)
        M_sel <- if(tmp > 1) 2*(tmp - 1) else 2*(KT - 1)
        M_vals <- -M_sel:M_sel
        R_vals <- R[abs(M_vals) + 1]
        lambda_vals <- lambda_compiled(M_vals/M_sel)
        ghat <- sum(lambda_vals * R_vals)
        Ghat <- sum(lambda_vals * R_vals * abs(M_vals))
        max(1, round((Ghat^2/ghat^2)^(1/3) * trainingSize^(1/3)))
    })
    blocksize <- round(mean(Block))
    # -------------------------------------------------------
    # 步骤2: 生成M个Bootstrap重复样本 T_n^[1], ..., T_n^[M]
    # -------------------------------------------------------
    T_boot <- matrix(0, nrow = M, ncol = nlag)
    
    for (b in 1:M) {
        # 生成Bootstrap样本
        total_blocks <- ceiling(N / blocksize) * 2
        starts <- sample.int(trainingSize, total_blocks, replace = TRUE)
        lengths <- rgeom(total_blocks, 1/blocksize) + 1
        idx_list <- lapply(1:total_blocks, function(j) {
            (starts[j] + 0:(lengths[j] - 1) - 1) %% trainingSize + 1
        })
        y_star <- do.call(rbind, lapply(idx_list, function(idx) {
            yt[idx, , drop = FALSE]
        }))[1:N, ]
        
        # 计算Bootstrap统计量
        for (i in 1:nlag) {
            m <- lags[i]
            w_vec <- weight_list[[as.character(m)]]
            stat_seq_star <- statL_compiled(y_star, trainingSize, m, alphaE, w_vec, kernal, delta0)
            T_boot[b, i] <- max(stat_seq_star, na.rm = TRUE)
        }
    }
    
    # -------------------------------------------------------
    # 步骤3: 计算所有样本（原始+M个Bootstrap）的p值
    # -------------------------------------------------------
    # 合并原始和Bootstrap统计量
    T_all <- rbind(T_original, T_boot)  # (M+1) x nlag
    
    # 计算p值（文章式2.3）
    p_matrix <- matrix(0, nrow = M + 1, ncol = nlag)
    
    for (j in 1:nlag) {
        for (i in 1:(M + 1)) {
        # 计算 p_{n,M}(T_{n,j}^{[i]}) =
        # (0.5 + Σ_{k=1}^M 1(T_{n,j}^{[k]} ≥ T_{n,j}^{[i]})) / (M+1)
            # 注意：这里i=1对应原始统计量，
            # i=2,...,M+1对应Bootstrap统计量
            count <- sum(T_boot[, j] >= T_all[i, j])
            p_matrix[i, j] <- (0.5 + count) / (M + 1)
        }
    }
    
    # -------------------------------------------------------
    # 步骤4: 计算组合统计量 W_{n,M}^{[i]}
    # -------------------------------------------------------
    W_all <- numeric(M + 1)
    
    for (i in 1:(M + 1)) {
        # 使用Fisher组合（文章式2.2）
        W_all[i] <- calc_fisher_stat(p_matrix[i, ])
    }
    
    # -------------------------------------------------------
    # 步骤5: 计算全局统计量的p值并做决策
    # -------------------------------------------------------
    
    # 原始统计量的组合统计量 W_{n,M}^{[0]}（对应i=1）
    W_original <- W_all[1]
    
    # Bootstrap样本的组合统计量（对应i=2,...,M+1）
    W_boot <- W_all[2:(M + 1)]
    
    # 计算组合检验的p值（文章式2.5）
    p_combined <- sum(W_boot >= W_original) / M
    
    # 组合检验决策
    rej_combined[sim] <- as.numeric(p_combined <= TypeI)
    
    # 计算单独检验的p值并做决策
    for (j in 1:nlag) {
        # 原始统计量在Bootstrap分布中的p值
        p_j <- (0.5 + sum(T_boot[, j] >= T_original[j])) / (M + 1)
        rej_individual[sim, j] <- as.numeric(p_j <= TypeI)
    }
    
    # 进度更新
    pb$tick()
    if (sim %% 100 == 0) {
        cat(sprintf("\n[进度] sim=%d | 组合拒绝率: %.2f%% | 单独拒绝率: %s",
                    sim, 
                    mean(rej_combined[1:sim]) * 100,
                    paste(sprintf("%.1f%%", colMeans(rej_individual[1:sim, , drop=FALSE]) * 100), collapse=", ")))
        gc()
    }
}

############### 5. 最终结果 ###############

cat("\n\n====================== 最终结果 ======================\n")

cat("参数配置:\n")
cat(sprintf("  数据维度: d=%d, DGP模型: %d\n", d, DGP))
cat(sprintf("  滞后阶数: {%s}\n", paste(lags, collapse=",")))
cat(sprintf("  核函数: %s, Bootstrap次数: M=%d\n", kernal, M))
cat(sprintf("  样本量: n=%d, 模拟次数: nsim=%d\n", trainingSize, nsim))
cat(sprintf("  显著性水平: α=%.3f\n", TypeI))
cat("-----------------------------------------------\n")

# 计算最终拒绝率
final_rej_individual <- colMeans(rej_individual)
final_rej_combined <- mean(rej_combined)

cat("拒绝率:\n")
for (i in 1:nlag) {
    cat(sprintf("  Lag %d (单独检验): %.4f (%.1f%%)\n", 
                lags[i], final_rej_individual[i], 
                final_rej_individual[i] * 100))
}
cat(sprintf("  组合检验: %.4f (%.1f%%)\n", 
            final_rej_combined, final_rej_combined * 100))
cat("-----------------------------------------------\n")

