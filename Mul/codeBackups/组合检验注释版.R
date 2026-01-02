rm(list = ls())
# 加载必要的库
library(progress)       #用于显示进度条
library(compiler)       #用于编译函数加速（JIT）
library(ggplot2)        #绘图（本代码主要用于计算）
library(mvtnorm)        #多元正态分布
library(Rcpp)           #R与C++接口
library(RcppArmadillo)  #线性代数库
library(MTS)            #多元时间序列包

# 设置工作目录并加载自定义的C++和R函数
# statL: 计算CUSUM型统计量的核心函数
# eucL: 计算欧氏距离矩阵
# dataGeneration: 生成模拟数据
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')
Rcpp::sourceCpp('statLMul.cpp')
Rcpp::sourceCpp("eucLMulNew.cpp")
source('dataGenerationNew.R')

############### 1. 辅助函数与参数设置 ###############

# [边界权重函数]
# 用于调整监测过程中不同时间点的权重，用于防止边缘效应
Qk <- function(s0, gamma0) {
    return(1 + s0) * (s0 / (1 + s0))^gamma0
}


# [Fisher 组合统计量]
# 核心公式: X^2 = -2 * sum(log(p_i))
# 作用: 将多个相关的检验结果（p值）合并为一个统计量
calc_fisher_stat <- function(pvals, weights = rep(1, length(pvals))) {
    # 这里的逻辑是：如果p值很小，log(p)是大的负数，-2*log(p)就是大的正数。
    # 统计量越大，表示越显著。
    return(-2 * sum(weights * log(pvals)))
}

# --- 基础参数配置 ---
d <- 2              # 数据维度
DGP <- 1            # 数据生成过程ID 
lags <- c(0, 1, 2)  # 我们要联合检验的滞后阶数集合
nlag <- length(lags)# 滞后阶数的数量

kernelTunta <- 1    # 核函数选择标记
trainingSize <- 200 # 训练集长度 (T)
monitoringPeriod <- 1 # 监测期长度倍数 (L)，总长度 N = T * (1+L)
M <- 50            # Bootstrap 重抽样次数 (决定了p值的精度)

gamma <- 0          # 权重函数参数
nsim <- 500         # 蒙特卡洛模拟的总次数
TypeI <- 0.05       # 显著性水平 alpha
alphaE <- 1         # Energy Kernel 的指数参数

# 确定核函数名称
if (kernelTunta == 1) kernal <- "gauss" else if (kernelTunta == 2) kernal <- "energy" else kernal <- "quadexp"

# [自适应参数 Delta0]
# 如果使用 quadexp 核，参数通常设为距离的中位数
if (kernelTunta == 3) {
    yt_init <- dataGeneration(monitoringPeriod, trainingSize, cptModel = DGP, dataDimension = d)
    E_init <- eucL(yt_init[1:trainingSize,], 0)
    delta0 <- median(E_init[E_init > 0])
} else {
    delta0 <- 1 * d # Gauss/Energy 核的默认参数
}

############### 2. 预计算权重向量 (性能优化) ##############
# 目的: vaha 权重只与样本长度和 gamma 有关，与具体数据无关。
# 移到循环外计算可以极大提升速度。
N <- (1 + monitoringPeriod) * trainingSize
weight_list <- list()

cat("预计算权重向量...\n")
for (lag_val in lags) {
    Tm <- trainingSize - lag_val
    k_vec <- 1:(N - trainingSize)
    # 计算监测统计量所需的权重序列
    vaha_vec <- (Tm + k_vec)^2 / (Tm * Qk(k_vec/Tm, gamma)^2)
    weight_list[[as.character(lag_val)]] <- vaha_vec
}

############### 3. 预编译与初始化 ##############

# 平顶核函数 (Flat-top kernel)，用于谱密度估计
lambda <- function(t) {
    abst <- abs(t)
    (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1))
}

# 使用 cmpfun 将 R 函数编译为字节码，提高执行速度
statL_compiled <- cmpfun(statL)
lambda_compiled <- cmpfun(lambda)

# 自动选择 Block Size 所需的常数
KT <- max(5, sqrt(log10(trainingSize)))
sqrt_log_ts <- 2 * sqrt(log10(trainingSize) / trainingSize)

# 初始化结果存储矩阵
rej_individual <- matrix(0, nrow = nsim, ncol = nlag)  
# 存储Lag 0, 1, 2 单独检验是否拒绝
rej_combined <- rep(0, nsim)                           
# 存储组合检验是否拒绝
colnames(rej_individual) <- paste0("lag", lags)

# 初始化进度条
pb <- progress_bar$new(total = nsim, 
                       format = "[:bar] :percent 耗时: :elapsedfull")

cat("\n开始模拟...\n")

############### 4. 主模拟循环 ##############

for (sim in 1:nsim) {
    
    # -------------------------------------------------------
    # 步骤1: 生成原始数据并计算原始统计量 T_n^[0]
    # -------------------------------------------------------
    yt <- dataGeneration(monitoringPeriod, trainingSize, cptModel = DGP, dataDimension = d)
    
    # 计算原始数据在 Lag 0, 1, 2 下的统计量最大值
    T_original <- numeric(nlag)
    for (i in 1:nlag) {
        m <- lags[i]
        w_vec <- weight_list[[as.character(m)]]
        # statL_compiled 返回整个监测过程的统计量序列，取 max 作为检验统计量
        stat_seq <- statL_compiled(yt, trainingSize, m, alphaE, w_vec, kernal, delta0)
        T_original[i] <- max(stat_seq, na.rm = TRUE)
    }
    
    # [Politis & White 方法] 自动选择最佳 Block Size
    # Stationary Bootstrap 依赖于 Block Size 来保留时间序列的依赖结构
    Block <- sapply(1:d, function(i) {
        ut <- yt[1:trainingSize, i]
        # 计算自相关系数
        R <- drop(acf(ut, lag.max = 2*KT, plot = FALSE, demean = FALSE)$acf)
        # 寻找显著的自相关截断点
        tmp <- which.max(abs(R[1:KT]) >= sqrt_log_ts)
        M_sel <- if(tmp > 1) 2*(tmp - 1) else 2*(KT - 1)
        # 谱密度估计计算...
        M_vals <- -M_sel:M_sel
        R_vals <- R[abs(M_vals) + 1]
        lambda_vals <- lambda_compiled(M_vals/M_sel)
        ghat <- sum(lambda_vals * R_vals)
        Ghat <- sum(lambda_vals * R_vals * abs(M_vals))
        # 经验公式计算 block size
        max(1, round((Ghat^2/ghat^2)^(1/3) * trainingSize^(1/3)))
    })
    blocksize <- round(mean(Block)) # 对所有维度的 block size 取平均
    
    # -------------------------------------------------------
    # 步骤2: 生成 M 个 Bootstrap 重复样本 T_n^[1], ..., T_n^[M]
    # -------------------------------------------------------
    T_boot <- matrix(0, nrow = M, ncol = nlag)
    

    # 这是一个平稳 Bootstrap 过程：随机选取起始点，块长度服从几何分布
    for (b in 1:M) {
        # 1. 确定块的数量
        total_blocks <- ceiling(N / blocksize) * 2
        # 2. 随机采样起始点
        starts <- sample.int(trainingSize, total_blocks, replace = TRUE)
        # 3. 随机生成块长度 (几何分布)
        lengths <- rgeom(total_blocks, 1/blocksize) + 1
        # 4. 构建索引列表 (使用取模运算实现循环连接)
        idx_list <- lapply(1:total_blocks, function(j) {
            (starts[j] + 0:(lengths[j] - 1) - 1) %% trainingSize + 1
        })
        # 5. 拼装数据 y_star
        y_star <- do.call(rbind, lapply(idx_list, function(idx) {
            yt[idx, , drop = FALSE]
        }))[1:N, ]
        
        # 计算该 Bootstrap 样本在所有 Lag 下的统计量
        for (i in 1:nlag) {
            m <- lags[i]
            w_vec <- weight_list[[as.character(m)]]
            stat_seq_star <- statL_compiled(y_star, trainingSize, m, alphaE, w_vec, kernal, delta0)
            T_boot[b, i] <- max(stat_seq_star, na.rm = TRUE)
        }
    }
    
    # -------------------------------------------------------
    # 步骤3: 计算所有样本（原始+M个Bootstrap）的经验 p 值
    # -------------------------------------------------------
    # 这一步非常关键：为了得到组合统计量 W 的分布，我们需要先把 Bootstrap 统计量也转换成 P 值
    
    # 合并：第一行是原始数据，后面 M 行是 Bootstrap 数据
    T_all <- rbind(T_original, T_boot)  # (M+1) x nlag
    
    p_matrix <- matrix(0, nrow = M + 1, ncol = nlag)
    
    for (j in 1:nlag) {
        # 对于第 j 个 Lag
        for (i in 1:(M + 1)) {
            # 计算第 i 个样本的统计量在“Bootstrap 分布”中的位置（P值）
            # P值定义：比当前统计量更大或相等的比例
            # 使用 T_boot 作为参考分布 (reference distribution)
            count <- sum(T_boot[, j] >= T_all[i, j])
            # (0.5 + count) / (M + 1) 是防止 p 值为 0 或 1 的校正公式
            p_matrix[i, j] <- (0.5 + count) / (M + 1)
        }
    }
    
    # -------------------------------------------------------
    # 步骤4: 计算组合统计量 W_{n,M}^{[i]}
    # -------------------------------------------------------
    W_all <- numeric(M + 1)
    
    for (i in 1:(M + 1)) {
        # 将每个样本在不同 Lag 下的 p 值组合起来
        # W = -2 * ln(p_lag0) - 2 * ln(p_lag1) - 2 * ln(p_lag2)
        W_all[i] <- calc_fisher_stat(p_matrix[i, ])
    }
    
    # -------------------------------------------------------
    # 步骤5: 计算全局统计量的 p 值并做决策
    # -------------------------------------------------------
    
    # 原始数据的组合统计量
    W_original <- W_all[1]
    
    # Bootstrap 样本的组合统计量分布
    W_boot <- W_all[2:(M + 1)]
    
    # 计算最终的组合 P 值：
    # 看 W_original 在 W_boot 分布中是否属于“极端大”的值
    # 注意：Fisher 统计量越大越显著，所以计算右尾概率
    p_combined <- sum(W_boot >= W_original) / M
    
    # 决策：如果 p <= 0.05，拒绝 H0
    rej_combined[sim] <- as.numeric(p_combined <= TypeI)
    
    # --- 同时计算单独 Lag 的拒绝情况 (用于对比) ---
    for (j in 1:nlag) {
        # 单独检验直接比较 T_original 和 T_boot 即可
        p_j <- (0.5 + sum(T_boot[, j] >= T_original[j])) / (M + 1)
        rej_individual[sim, j] <- as.numeric(p_j <= TypeI)
    }
    
    # 进度更新
    pb$tick()
    # 每100次打印一次中间结果
    if (sim %% 100 == 0) {
        cat(sprintf("\n[进度] sim=%d | 组合拒绝率: %.2f%% | 单独拒绝率: %s",
                    sim, 
                    mean(rej_combined[1:sim]) * 100,
                    paste(sprintf("%.1f%%", colMeans(rej_individual[1:sim, , drop=FALSE]) * 100), collapse=", ")))
        gc() # 垃圾回收，释放内存
    }
}

############### 5. 最终结果输出 ###############

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

