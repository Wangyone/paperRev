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

########################### 辅助函数 ###########################

# 边界函数
Qk <- function(s0, gamma0) {
    return(1 + s0) * (s0 / (1 + s0))^gamma0
}

########################### 参数设置 ###########################

d <- 2  # 数据维度
DGP <- 1  # 数据生成过程
lags <- c(0, 1, 2)  # 要组合的lag值
nlag <- length(lags)  # lag的数量
kernelTunta <- 3  # 1 for gauss, 2 for energy, 3 for quadexp

# 数据长度参数
trainingSize <- 100
monitoringPeriod <- 1

# 统计量参数
gamma <- 0
deltagauss <- 1 * d

# 模拟参数
bootRepet <- 500  # Bootstrap次数

# 核函数选择
if (kernelTunta == 1) {
    kernal <- "gauss"
} else if (kernelTunta == 2) {
    kernal <- "energy"
} else if (kernelTunta == 3) {
    kernal <- "quadexp"
}

# 计算delta0参数（使用lag=0）
if (kernelTunta == 1 || kernelTunta == 2) {
    delta0 <- deltagauss  # gauss/energy
} else if (kernelTunta == 3) {
    yt <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
    E <- eucL(yt[1:trainingSize, ], 0)  # 使用lag=0计算
    medianE <- median(E[E > 0])
    delta0 <- medianE
}

# 其他参数
alphaE <- 1  # energy参数
N <- (1 + monitoringPeriod) * trainingSize
TypeI <- 0.05  # 显著性水平

########################### 计算各lag的权重 ###########################

# 计算每个lag的权重序列
weight_list <- list()
for (lag in lags) {
    Tm <- trainingSize - lag
    k <- 1:(N - trainingSize)
    # 权重函数
    vaha <- (Tm + k)^2 / (Tm * Qk(k/Tm, gamma)^2)
    weight_list[[as.character(lag)]] <- vaha
}

########################### 组合函数 ###########################

# Fisher组合函数
fisher_combine <- function(pvals) {
    # 避免p值为0或1导致数值问题
    pvals <- pmin(pmax(pvals, 1e-10), 1 - 1e-10)
    return(-2 * sum(log(pvals)))
}

# 计算p值函数（平滑版本）
calc_pvalue <- function(stat_original, stat_bootstrap, k) {
    if (k == 1) {
        return(0.5)  # 第一次Bootstrap，设为0.5
    } else {
        # 平滑的经验p值（公式2.3）
        return((0.5 + sum(stat_bootstrap[1:k] >= stat_original)) / (k + 1))
    }
}

########################### 预编译函数 ###########################

# flat-top lag-window函数
lambda <- function(t) {
    abst <- abs(t)
    (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1))
}

lambda_compiled <- cmpfun(lambda)

# 编译关键函数
statL_compiled <- cmpfun(statL)

# 预计算
KT <- max(5, sqrt(log10(trainingSize)))
sqrt_log_ts <- 2 * sqrt(log10(trainingSize) / trainingSize)

########################### wrap-speed Bootstrap组合检验 ###########################

# 初始化结果存储
Tn <- matrix(0, nrow = bootRepet, ncol = nlag)  # 原始统计量
Tn_star <- matrix(0, nrow = bootRepet, ncol = nlag)  # Bootstrap统计量
rej_individual <- matrix(0, nrow = bootRepet, ncol = nlag)  # 各单独检验的拒绝标志
rej_combined <- rep(0, bootRepet)  # 组合检验的拒绝标志

colnames(Tn) <- paste0("lag", lags)
colnames(Tn_star) <- paste0("lag", lags)
colnames(rej_individual) <- paste0("lag", lags)

# 创建进度条
pb <- progress_bar$new(
    format = "[:bar] :percent 已用: :elapsed 剩余: :eta",
    total = bootRepet,
    width = 60
)

# 主模拟循环（wrap-speed方法）
for (s in 1:bootRepet) {
    
    # 1. 生成原始数据
    yt <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
    
    # 2. 计算原始统计量（所有lag）
    for (i in 1:nlag) {
        lag_val <- lags[i]
        weight_vec <- weight_list[[as.character(lag_val)]]
        stat_seq <- statL_compiled(yt, trainingSize, m = lag_val, alphaE, weight_vec, kernal, delta0)
        Tn[s, i] <- max(stat_seq, na.rm = TRUE)
    }
    
    # 3. 确定blocksize（使用第一个变量的自相关）
    Block <- sapply(1:d, function(i) {
        ut <- yt[1:trainingSize, i]
        R <- drop(acf(ut, lag.max = 2 * KT, plot = FALSE, demean = FALSE)$acf)
        tmp <- which.max(abs(R[1:KT]) >= sqrt_log_ts)
        M <- if (tmp > 1) 2 * (tmp - 1) else 2 * (KT - 1)
        M_vals <- -M:M
        R_vals <- R[abs(M_vals) + 1]
        lambda_vals <- lambda_compiled(M_vals / M)
        ghat <- sum(lambda_vals * R_vals)
        Ghat <- sum(lambda_vals * R_vals * abs(M_vals))
        max(1, round((Ghat^2 / ghat^2)^(1/3) * trainingSize^(1/3)))
    })
    
    blocksize <- round(mean(Block))
    
    # 4. 生成第s次Bootstrap样本
    total_blocks <- ceiling(N / blocksize) * 2
    starts <- sample.int(trainingSize, total_blocks, replace = TRUE)
    lengths <- rgeom(total_blocks, 1 / blocksize) + 1
    idx_list <- lapply(1:total_blocks, function(j) {
        (starts[j] + 0:(lengths[j] - 1) - 1) %% trainingSize + 1
    })
    y_star <- do.call(rbind, lapply(idx_list, function(idx) {
        yt[idx, , drop = FALSE]
    }))[1:N, ]
    
    # 5. 计算Bootstrap统计量
    for (i in 1:nlag) {
        lag_val <- lags[i]
        weight_vec <- weight_list[[as.character(lag_val)]]
        stat_seq <- statL_compiled(y_star, trainingSize, m = lag_val, alphaE, weight_vec, kernal, delta0)
        Tn_star[s, i] <- max(stat_seq, na.rm = TRUE)
    }
    
    # 6. wrap-speed方法：更新各检验的拒绝决策
    if (s > 1) {
        # 各单独检验的拒绝决策
        for (i in 1:nlag) {
            # 计算临界值（1-TypeI分位数）
            crit_val <- quantile(Tn_star[1:s, i], 1 - TypeI, na.rm = TRUE)
            rej_individual[s, i] <- as.numeric(Tn[s, i] > crit_val)
        }
        
        # 组合检验的拒绝决策
        # 6.1 计算各检验的p值
        pvals <- numeric(nlag)
        for (i in 1:nlag) {
            pvals[i] <- calc_pvalue(Tn[s, i], Tn_star[1:s, i], s)
        }
        
        # 6.2 计算组合统计量
        comb_stat <- fisher_combine(pvals)
        
        # 6.3 计算Bootstrap样本的组合统计量
        comb_stats_bootstrap <- numeric(s)
        for (b in 1:s) {
            # 计算每个Bootstrap样本的p值
            pvals_b <- numeric(nlag)
            for (i in 1:nlag) {
                pvals_b[i] <- calc_pvalue(Tn_star[b, i], Tn_star[1:s, i], s)
            }
            comb_stats_bootstrap[b] <- fisher_combine(pvals_b)
        }
        
        # 6.4 计算组合检验的p值
        pval_combined <- calc_pvalue(comb_stat, comb_stats_bootstrap, s)
        
        # 6.5 记录拒绝决策
        rej_combined[s] <- as.numeric(pval_combined < TypeI)
    }
    
    # 7. 输出进度
    pb$tick()
    
    # 定期输出中间结果
    if (s %% 100 == 0) {
        # 计算当前拒绝率
        if (s > 1) {
            current_rej_individual <- colMeans(rej_individual[1:s, , drop = FALSE])
            current_rej_combined <- mean(rej_combined[1:s])
            
            cat(sprintf("\n=== 完成 %d/%d 次模拟 ===\n", s, bootRepet))
            cat("当前拒绝率（水平=", TypeI, "）：\n")
            for (i in 1:nlag) {
                cat(sprintf("  %-10s: %.3f\n", colnames(rej_individual)[i], current_rej_individual[i]))
            }
            cat(sprintf("  %-10s: %.3f\n", "combined", current_rej_combined))
        }
    }
    
    # 定期垃圾回收
    if (s %% 100 == 0) gc()
}

########################### 结果分析 ###########################

cat("\n\n=================== 最终结果 ===================\n")
cat("参数设置：\n")
cat(sprintf("  数据维度 d = %d\n", d))
cat(sprintf("  训练样本大小 T = %d\n", trainingSize))
cat(sprintf("  监测周期 L = %d\n", monitoringPeriod))
cat(sprintf("  核函数: %s\n", kernal))
cat(sprintf("  delta0 = %.4f\n", delta0))
cat(sprintf("  组合的lags: %s\n", paste(lags, collapse = ", ")))

# 计算最终拒绝率（只使用最后100次以确保稳定）

start_idx <- 1
final_rej_individual <- colMeans(rej_individual[start_idx:bootRepet, , drop = FALSE])
final_rej_combined <- mean(rej_combined[start_idx:bootRepet])

cat("\n拒绝率（水平=", TypeI, "，基于最后", (bootRepet - start_idx + 1), "次模拟）：\n")
for (i in 1:nlag) {
    cat(sprintf("  %-10s: %.3f (%.1f%%)\n", 
                colnames(rej_individual)[i], 
                final_rej_individual[i],
                final_rej_individual[i] * 100))
}
cat(sprintf("  %-10s: %.3f (%.1f%%)\n", 
            "combined", 
            final_rej_combined,
            final_rej_combined * 100))



########################### 结果可视化 ###########################

# # 1. 拒绝率柱状图
# rej_df <- data.frame(
#     Test = c(colnames(rej_individual), "combined"),
#     RejectionRate = c(final_rej_individual * 100, final_rej_combined * 100),
#     SE = c(se_individual * 100, se_combined * 100)
# )
# 
# p1 <- ggplot(rej_df, aes(x = Test, y = RejectionRate, fill = Test)) +
#     geom_bar(stat = "identity", alpha = 0.8) +
#     geom_errorbar(aes(ymin = RejectionRate - 1.96 * SE, ymax = RejectionRate + 1.96 * SE),
#                   width = 0.2) +
#     geom_hline(yintercept = TypeI * 100, linetype = "dashed", color = "red", size = 1) +
#     labs(title = paste0("拒绝率比较 (水平 = ", TypeI * 100, "%)"),
#          subtitle = paste0("数据生成过程: DGP", DGP, ", d=", d),
#          x = "检验方法", y = "拒绝率 (%)") +
#     theme_minimal() +
#     theme(legend.position = "none",
#           axis.text.x = element_text(angle = 45, hjust = 1),
#           plot.title = element_text(hjust = 0.5),
#           plot.subtitle = element_text(hjust = 0.5))
# 
# print(p1)
# 
# # 2. 拒绝率随时间变化图（wrap-speed方法特性）
# if (bootRepet >= 100) {
#     # 计算累积拒绝率
#     cum_rej_individual <- matrix(0, nrow = bootRepet, ncol = nlag)
#     cum_rej_combined <- numeric(bootRepet)
#     
#     for (s in 1:bootRepet) {
#         if (s == 1) {
#             cum_rej_individual[1, ] <- rej_individual[1, ]
#             cum_rej_combined[1] <- rej_combined[1]
#         } else {
#             cum_rej_individual[s, ] <- colMeans(rej_individual[1:s, , drop = FALSE])
#             cum_rej_combined[s] <- mean(rej_combined[1:s])
#         }
#     }
#     
#     # 创建数据框
#     cum_rej_df <- data.frame(
#         Iteration = rep(1:bootRepet, nlag + 1),
#         RejectionRate = c(as.vector(cum_rej_individual), cum_rej_combined),
#         Test = c(rep(colnames(rej_individual), each = bootRepet), rep("combined", bootRepet))
#     )
#     
#     p2 <- ggplot(cum_rej_df, aes(x = Iteration, y = RejectionRate, color = Test)) +
#         geom_line(size = 0.8) +
#         geom_hline(yintercept = TypeI, linetype = "dashed", color = "red", size = 1) +
#         labs(title = "拒绝率随Bootstrap次数变化",
#              subtitle = "wrap-speed方法：每次Bootstrap后更新临界值",
#              x = "Bootstrap次数", y = "累积拒绝率") +
#         theme_minimal() +
#         theme(legend.position = "bottom",
#               plot.title = element_text(hjust = 0.5),
#               plot.subtitle = element_text(hjust = 0.5))
#     
#     print(p2)
# }

