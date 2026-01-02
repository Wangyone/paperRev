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
source('dataGenerationSupply.R')

############### 1. 辅助函数与参数设置 ###############

# 边界权重函数
Qk <- function(s0, gamma0) {
    return(1 + s0) * (s0 / (1 + s0))^gamma0
}

# Fisher 组合统计量计算函数
calc_fisher_stat <- function(pvals) {
    # 避免 log(0)
    pvals <- pmax(pvals, 1e-10) 
    return(-2 * rowSums(log(pvals)))
}

# 参数配置
d <- 2  # 维度
DGP <- 3  # 数据生成模型
lags <- c(0, 1, 2)  # 组合的 Lags

nlag <- length(lags)
kernelTunta <- 1  

trainingSize <- 100
monitoringPeriod <- 1

gamma <- 0
nsim <- 1000  # Bootstrap 次数 (Warp-Speed)
TypeI <- 0.05 # 显著性水平
alphaE <- 1

# 核函数配置
if (kernelTunta == 1) kernal <- "gauss" else if (kernelTunta == 2) kernal <- "energy" else kernal <- "quadexp"

# 预计算kernel_3的参数Delta0 (基于 Lag=0)
yt_init <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
if (kernelTunta == 3) {
    E_init <- eucL(yt_init[1:trainingSize,], 0)
    delta0 <- median(E_init[E_init > 0])
} else {
    delta0 <- 1 * d 
}

############### 2. 预计算不同 Lag下统计量 的权重###############
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

# 存储矩阵
Tn <- matrix(0, nrow = nsim, ncol = nlag)      # 原始统计量
Tn_star <- matrix(0, nrow = nsim, ncol = nlag) # Bootstrap 统计量
colnames(Tn) <- paste0("lag", lags)
colnames(Tn_star) <- paste0("lag", lags)

pb <- progress_bar$new(total = nsim, 
                       format = "[:bar] :percent 耗时: :elapsedfull")

############### 4. 生成所有Bootstrap样本 ###############

cat("\n生成所有 Bootstrap 样本...\n")


# 生成所有原始数据和Bootstrap数据
for (s in 1:nsim) {
    # -------------------------------------------------------
    # A. 生成原始数据与计算原始统计量 (Tn)
    # -------------------------------------------------------
    yt <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
    
    # 确定blocksize（使用第一组数据的自相关）
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
    # cat(sprintf("使用的 blocksize: %d\n", blocksize))
    
    # 计算不同lag下的原始统计量
    for (i in 1:nlag) {
        m <- lags[i]
        w_vec <- weight_list[[as.character(m)]]
        stat_seq <- statL_compiled(yt, trainingSize, m, alphaE, w_vec, kernal, delta0)
        Tn[s, i] <- max(stat_seq, na.rm = TRUE)
    }
    
    # -------------------------------------------------------
    # B. 生成Bootstrap样本并计算统计量 (Tn_star)
    # -------------------------------------------------------
    total_blocks <- ceiling(N / blocksize) * 2
    starts <- sample.int(trainingSize, total_blocks, replace = TRUE)
    lengths <- rgeom(total_blocks, 1 / blocksize) + 1
    idx_list <- lapply(1:total_blocks, function(j) {
        (starts[j] + 0:(lengths[j] - 1) - 1) %% trainingSize + 1
    })
    y_star <- do.call(rbind, lapply(idx_list, function(idx) {
        yt[idx, , drop = FALSE]
    }))[1:N, ]
    
    for (i in 1:nlag) {
        m <- lags[i]
        w_vec <- weight_list[[as.character(m)]]
        stat_seq_star <- statL_compiled(y_star, trainingSize, m, alphaE, w_vec, kernal, delta0)
        Tn_star[s, i] <- max(stat_seq_star, na.rm = TRUE)
    }
    
    pb$tick()
    if (s %% 200 == 0) gc()
}

############### 5. Wrap-Speed 方法计算拒绝率 ###############

cat("\n使用 Wrap-Speed 方法计算拒绝率...\n")

# 初始化存储
rejection_history_individual <- matrix(0, nrow = nsim, ncol = nlag)
rejection_history_combined <- rep(0, nsim)
colnames(rejection_history_individual) <- paste0("lag", lags)

# Warp-Speed: 每次Bootstrap后更新临界值和拒绝决策
for (s in 1:nsim) {
    # -------------------------------------------------------
    # A. 计算各lag的p值（基于前s次Bootstrap）
    # -------------------------------------------------------
    p_obs <- numeric(nlag)
    p_star_matrix <- matrix(0, nrow = s, ncol = nlag)
    
    for (i in 1:nlag) {
        # 1. 原始统计量的p值（相对于s次Bootstrap）
        p_obs[i] <- (sum(Tn_star[1:s, i] >= Tn[s, i]) + 1) / (s + 1)
        
        # 2. Bootstrap样本的p值（用于构建组合统计量的分布）
        for (b in 1:s) {
            # 每个Bootstrap样本在s个样本中的p值
            p_star_matrix[b, i] <- (sum(Tn_star[1:s, i] >= Tn_star[b, i]) + 1) / (s + 1)
        }
    }
    
    # -------------------------------------------------------
    # B. Fisher组合
    # -------------------------------------------------------
    W_obs <- -2 * sum(log(p_obs))
    
    # Bootstrap样本的组合统计量
    W_star <- numeric(s)
    for (b in 1:s) {
        W_star[b] <- -2 * sum(log(p_star_matrix[b, ]))
    }
    
    # -------------------------------------------------------
    # C. 确定临界值与拒绝决策
    # -------------------------------------------------------
    # 各单独检验的拒绝决策
    for (i in 1:nlag) {
        crit_val_individual <- quantile(Tn_star[1:s, i], 1 - TypeI, na.rm = TRUE)
        rejection_history_individual[s, i] <- as.numeric(Tn[s, i] > crit_val_individual)
    }
    
    # 组合检验的拒绝决策
    if (s > 1) {
        crit_val_combined <- quantile(W_star, 1 - TypeI, na.rm = TRUE)
        rejection_history_combined[s] <- as.numeric(W_obs > crit_val_combined)
    }
    
    # -------------------------------------------------------
    # D. 中间过程监测
    # -------------------------------------------------------
    if (s %% 100 == 0) {
        # 计算当前累积拒绝率
        if (s > 10) {  # 避免初期样本太少
            curr_rej_individual <- colMeans(rejection_history_individual[1:s, , drop = FALSE])
            curr_rej_combined <- mean(rejection_history_combined[1:s])
            
            cat(sprintf("\n[Monitor] s = %d | Combined Rej Rate: %.2f%%", s, curr_rej_combined * 100))
            cat(sprintf(" | Individual: %s", 
                        paste(sprintf("%.1f%%", curr_rej_individual * 100), collapse = ", ")))
        }
    }
}

############### 6. 最终结果计算与分析 ###############

cat("\n\n====================== 最终结果 ======================\n")

# 使用所有nsim次模拟计算最终拒绝率
final_rej_individual <- colMeans(rejection_history_individual)
final_rej_combined <- mean(rejection_history_combined)


cat("参数配置:\n")
cat(sprintf("  DGP: %d, 维度: d=%d, 核函数: %s\n", DGP, d, kernal))
cat(sprintf("  Lags: {%s}, nsim: %d, 水平: %.3f\n", 
            paste(lags, collapse=","), nsim, TypeI))
cat("-----------------------------------------------\n")
cat("拒绝率 (Wrap-Speed 方法):\n")
for (i in 1:nlag) {
    cat(sprintf("  Lag %d: %.3f (%.1f%%) \n", 
                lags[i], final_rej_individual[i], 
                final_rej_individual[i] * 100))
}
cat(sprintf("  组合检验: %.3f (%.1f%%) \n", 
            final_rej_combined, final_rej_combined * 100))
cat("-----------------------------------------------\n")



############### 7. 可视化 ###############

# 1. 拒绝率随Bootstrap次数变化图
cum_rej_individual <- apply(rejection_history_individual, 2, function(x) {
    cumsum(x) / 1:nsim
})
cum_rej_combined <- cumsum(rejection_history_combined) / 1:nsim

# 创建数据框
plot_df <- data.frame(
    Iteration = rep(1:nsim, nlag + 1),
    RejectionRate = c(as.vector(cum_rej_individual), cum_rej_combined),
    Test = c(rep(colnames(rejection_history_individual), each = nsim), 
             rep("Combined", nsim))
)

# 每50次取一个点以减少图形复杂度
thin_idx <- seq(1, nsim, by = 50)
if (max(thin_idx) < nsim) thin_idx <- c(thin_idx, nsim)

plot_df_thin <- plot_df[plot_df$Iteration %in% thin_idx, ]

p1 <- ggplot(plot_df_thin, aes(x = Iteration, y = RejectionRate, color = Test)) +
    geom_line(size = 0.8) +
    geom_hline(yintercept = TypeI, linetype = "dashed", color = "red", size = 0.8) +
    labs(title = "Wrap-Speed方法：拒绝率随Bootstrap次数变化",
         subtitle = sprintf("DGP=%d, d=%d, Lags={%s}", DGP, d, paste(lags, collapse=",")),
         x = "Bootstrap次数", y = "累积拒绝率") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9)) +
    scale_color_brewer(palette = "Set1")

print(p1)

# 2. 最终拒绝率柱状图
final_df <- data.frame(
    Test = c(colnames(rejection_history_individual), "Combined"),
    RejectionRate = c(final_rej_individual, final_rej_combined) * 100
)

p2 <- ggplot(final_df, aes(x = Test, y = RejectionRate, fill = Test)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = TypeI * 100, linetype = "dashed", 
               color = "red", size = 1) +
    labs(title = "最终拒绝率比较",
         subtitle = sprintf("B=%d, 水平=%.1f%%", nsim, TypeI * 100),
         x = "检验方法", y = "拒绝率 (%)") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set2")

print(p2)

