rm(list = ls())
setwd("D:/Desktop/paperRev/JointCUSUM/Mul/") # 请设置你的工作路径

library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(foreach)
library(doParallel)

# --- 1. 全局配置与参数 ---

# 组合方法选择器
# 1 = Fisher (对数加和, 适合弥散信号)
# 2 = Tippett (Min-P, 适合稀疏信号, 如 Model 306)
# 3 = Stouffer (Z分数加和, 适合正态化权衡)
# 4 = Max-P (Pearson, 通常用于检测一致性, 较少用于变点)
CombiningMethod <- c(1, 2, 3) 

# 基础参数
DGP <- 306         # 模型 ID (如 306 为稀疏信号模型)
lags <- c(0, 1, 2) # 组合的 Lags
nlag <- length(lags)
tunta <- 1         # 核函数: 1=Gauss
delta0 <- 1
gamma <- 0
trainSize <- 100   # 样本大小 T
monitoringPeriod <- 1 
TypeI <- 0.05      # 显著性水平
alphaE <- 1

nsim <- 800        # 模拟次数
M <- 299           # Bootstrap 次数

# 边界函数
Qk <- function(s0, gamma0){
    (1 + s0) * (s0 / (1 + s0))^gamma0
}

# 组合统计量计算核心函数
calc_combined_stats_vec <- function(pvals, methods_indices) {
    # 预处理：防止数值溢出
    pvals <- pmax(pvals, 1e-10)
    pvals <- pmin(pvals, 1 - 1e-10)
    
    results <- numeric(0)
    
    # Method 1: Fisher
    if(1 %in% methods_indices) {
        fisher_stat <- -2 * sum(log(pvals))
        results <- c(results, fisher_stat)
    }
    
    # Method 2: Tippett (1 - min_p)
    # 使用 1-min 是为了保持"值越大越显著"的统一逻辑
    if(2 %in% methods_indices) {
        tippett_stat <- 1 - min(pvals)
        results <- c(results, tippett_stat)
    }
    
    # Method 3: Stouffer
    if(3 %in% methods_indices) {
        stouffer_stat <- sum(qnorm(1 - pvals)) / sqrt(length(pvals))
        results <- c(results, stouffer_stat)
    }
    
    # Method 4: Max-P (Pearson)
    if(4 %in% methods_indices) {
        max_p_stat <- -2 * sum(log(1 - pvals))
        results <- c(results, max_p_stat)
    }
    
    return(results)
}

# --- 2. 预计算 (主进程) ---

N <- (1 + monitoringPeriod) * trainSize
kernal <- switch(tunta, "1"="gauss", "2"="energy", "3"="quadexp")

# 预计算权重向量
weight_list <- list()
for(lag_val in lags) {
    Tm_lag <- trainSize - lag_val
    k_vec <- 1:(N - trainSize)
    w <- (Tm_lag + k_vec)^2 / (Tm_lag * Qk(k_vec/Tm_lag, gamma)^2)
    weight_list[[as.character(lag_val)]] <- w
}

KT <- max(5, sqrt(log10(trainSize)))  
sqrt_log_T <- 2 * sqrt(log10(trainSize)/trainSize) 

# --- 3. 启动并行环境 ---

# 自动检测核心数
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("已启动并行计算，使用核心数: %d\n", num_cores))

# 配置子进程环境 (必须步骤)
clusterEvalQ(cl, {
    library(mvtnorm)
    library(Rcpp)
    library(RcppArmadillo)
    library(MTS)
    library(compiler)
    
    # === 请确保此处路径正确 ===
    setwd("D:/Desktop/paperRev/JointCUSUM/Mul/") 
    
    # 加载源码
    sourceCpp('statLMul.cpp')
    sourceCpp('eucLMulNew.cpp')
    if(file.exists("Gendata4PaperNew.R")) source("Gendata4PaperNew.R")
    
    # 编译辅助函数
    lambda <- function(t) {
        abst <- abs(t)
        (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
    }
    lambda_compiled <- cmpfun(lambda)
    statL_compiled <- cmpfun(statL) 
})

# --- 4. 主模拟循环 ---

cat("开始多方法组合检验模拟...\n")
cat(sprintf("选定的方法 ID: %s\n", paste(CombiningMethod, collapse=",")))
startTime <- Sys.time()

# 确定返回列数: nlag (单独) + 选定的组合方法数
num_comb_methods <- length(CombiningMethod)
total_cols <- nlag + num_comb_methods

final_results <- foreach(s = 1:nsim, .combine = rbind, 
                         .export = c("weight_list", "lags", "nlag", "Qk",
                                     "trainSize", "monitoringPeriod", "DGP", "N", "alphaE", 
                                     "kernal", "delta0", "M", "KT", "sqrt_log_T", 
                                     "calc_combined_stats_vec", "CombiningMethod")) %dopar% {
                                         
                                         # === [Safety Check 1] 数据生成 ===
                                         yt <- tryCatch({
                                             gendataPaper(monitoringPeriod, trainSize, DGP)
                                         }, error = function(e) return(NULL))
                                         
                                         if(is.null(yt) || any(!is.finite(yt))) return(rep(NA, total_cols))
                                         
                                         d <- ncol(yt)
                                         
                                         # 1. 计算原始统计量
                                         T_orig <- numeric(nlag)
                                         for (j in 1:nlag) {
                                             lag_val <- lags[j]
                                             seq_val <- statL_compiled(yt, trainSize, m=lag_val, alphaE, 
                                                                       weight_list[[as.character(lag_val)]], kernal, delta0)
                                             val <- max(seq_val)
                                             T_orig[j] <- if(is.finite(val)) val else -Inf 
                                         }
                                         
                                         # 2. 计算 Block Size (带容错)
                                         blocksize <- 1 
                                         try({
                                             ut <- yt[1:trainSize, 1]
                                             if(var(ut) > 1e-10) {
                                                 acf_out <- acf(ut, lag.max=2*KT, plot=FALSE, demean=TRUE) 
                                                 R <- drop(acf_out$acf)
                                                 if(all(is.finite(R))) {
                                                     tmp <- which(abs(R[1:KT]) < sqrt_log_T)
                                                     M_lag <- if (length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
                                                     M_range <- -M_lag:M_lag
                                                     lam_vals <- lambda_compiled(M_range/M_lag)
                                                     ghat <- sum(lam_vals * R[abs(M_range) + 1])
                                                     Ghat <- sum(lam_vals * R[abs(M_range) + 1] * abs(M_range))
                                                     D.SB <- 2 * ghat^2
                                                     if(!is.na(D.SB) && D.SB > 1e-8) {
                                                         blocksize <- max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainSize^(1/3)))
                                                     }
                                                 }
                                             }
                                         }, silent = TRUE)
                                         if(is.na(blocksize) || blocksize < 1) blocksize <- 5
                                         
                                         # 3. Bootstrap 循环
                                         T_boot <- matrix(NA, nrow = M, ncol = nlag)
                                         num_blocks_est <- ceiling(N / blocksize) + 10
                                         
                                         for (b in 1:M) {
                                             len_vec <- rgeom(num_blocks_est, 1/blocksize) + 1
                                             cumsum_len <- cumsum(len_vec)
                                             cut_point <- which(cumsum_len >= N)[1]
                                             if(is.na(cut_point)) cut_point <- length(len_vec)
                                             valid_lens <- len_vec[1:cut_point]
                                             
                                             start_points <- sample.int(trainSize, length(valid_lens), replace = TRUE)
                                             idx_list <- lapply(seq_along(valid_lens), function(k) {
                                                 (start_points[k] + 0:(valid_lens[k]-1) - 1) %% trainSize + 1
                                             })
                                             idx_vec <- unlist(idx_list)[1:N]
                                             y_star <- yt[idx_vec, , drop=FALSE]
                                             
                                             for (j in 1:nlag) {
                                                 lag_val <- lags[j]
                                                 seq_star <- statL_compiled(y_star, trainSize, m=lag_val, alphaE, 
                                                                            weight_list[[as.character(lag_val)]], kernal, delta0)
                                                 val <- max(seq_star)
                                                 T_boot[b, j] <- if(is.finite(val)) val else -Inf
                                             }
                                         }
                                         
                                         # 4. P 值计算 (Rank Method)
                                         T_total <- rbind(T_orig, T_boot)
                                         p_mat <- matrix(0, nrow = M + 1, ncol = nlag)
                                         
                                         for (j in 1:nlag) {
                                             ranks <- rank(-T_total[, j], ties.method = "max", na.last = TRUE)
                                             p_mat[, j] <- (ranks + 1) / (M + 2)
                                         }
                                         
                                         # 5. 多种组合统计量计算
                                         # 对每一行计算选定方法的统计量 -> 得到 (M+1) x num_methods 的矩阵
                                         W_mat <- t(apply(p_mat, 1, calc_combined_stats_vec, methods_indices=CombiningMethod))
                                         
                                         # 6. 决策
                                         decisions <- numeric(total_cols)
                                         
                                         # 6.1 单独 Lag 决策 (前 nlag 列)
                                         decisions[1:nlag] <- as.numeric(p_mat[1, ] < TypeI)
                                         
                                         # 6.2 组合方法决策 (后 num_comb_methods 列)
                                         # 比较 W_orig (第1行) 和 W_boot (第2:M+1行)
                                         if(ncol(W_mat) > 0) {
                                             for(k in 1:ncol(W_mat)) {
                                                 W_col <- W_mat[, k]
                                                 # 计算组合 P 值
                                                 p_comb <- (sum(W_col[2:(M+1)] >= W_col[1]) + 1) / (M + 1)
                                                 decisions[nlag + k] <- as.numeric(p_comb < TypeI)
                                             }
                                         }
                                         
                                         decisions
                                     }

# --- 5. 结果处理与输出 ---

stopCluster(cl)

# 移除无效行
final_results <- na.omit(final_results)

elapsed <- difftime(Sys.time(), startTime, units="secs")
cat(sprintf("\n模拟完成! 有效样本: %d/%d, 耗时: %.2f 秒\n", nrow(final_results), nsim, as.numeric(elapsed)))

# 构建列名
method_names <- c("Fisher", "Tippett", "Stouffer", "MaxP")
selected_names <- method_names[CombiningMethod]
all_colnames <- c(paste0("Lag_", lags), selected_names)
colnames(final_results) <- all_colnames

# 计算拒绝率 (%)
rejection_rates <- colMeans(final_results) * 100

cat("\n================ 结果摘要 ================\n")
cat(sprintf("DGP模型: %d | nsim: %d | M: %d\n", DGP, nsim, M))
cat("------------------------------------------\n")

res_df <- data.frame(
    Test_Method = names(rejection_rates),
    Rejection_Rate_Percent = rejection_rates
)
# 格式化输出
print(res_df, row.names = FALSE)
cat("==========================================\n")