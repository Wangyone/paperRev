rm(list = ls())
######################### Simulation: Fisher Combination Test (Full Bootstrap) #########################

# --- 1. 环境设置 ---
# 请确保工作目录下包含:
# 1. statLMax.cpp (新的 C++ 文件)
# 2. Gendata4PaperNew.R (数据生成函数)
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(foreach)
library(doParallel)

# 检查依赖
if(!file.exists("Gendata4PaperNew.R")) stop("Gendata4PaperNew.R missing")
if(!file.exists("statLMax.cpp")) stop("statLMax.cpp missing")

# =========================================================================
# 2. 核心组合统计量函数 
# =========================================================================
calc_combined_stats_vec <- function(pvals, methods_indices, weights) {
    if(any(is.na(pvals))) return(rep(NA, length(methods_indices)))
    
    # 预处理防止 log(0)
    pvals <- pmax(pvals, 1e-10)
    pvals <- pmin(pvals, 1 - 1e-10)
    
    results <- numeric(0)
    
    # Method 1: Weighted Fisher
    if(1 %in% methods_indices) {
        fisher_stat <- -2 * sum(weights * log(pvals))
        results <- c(results, fisher_stat)
    }
    return(results)
}

# =========================================================================
# 3. 参数设置
# =========================================================================
m_list      <- c(100, 200)      
L_list      <- c(1, 2, 3)        
model_list  <- 1:8              

# 方法配置
CombiningMethod <- c(1)        # 1 = Fisher
lags            <- c(0, 1, 2)  
nlag            <- length(lags)
lag_weights     <- c(1.0, 1.0, 1.0) 

# 模拟参数
gamma    <- 0
nsim     <- 100       # 外层循环次数 (建议正式运行时设为 1000)
M        <- 200       # 内层 Bootstrap 次数 (计算 P 值所需，建议 200+)
TypeI    <- 0.05
kernal   <- "gauss"
delta0   <- 1.0
alphaE   <- 1

# 辅助函数
Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }

cat("==========================================================\n")
cat("Start Simulation (Full Bootstrap - Fisher) with statL_max\n")
cat(sprintf("Models: N1-N8 | m: %s | L: %s | Sim: %d | Boot(M): %d\n", 
            paste(m_list, collapse=","), paste(L_list, collapse=","), nsim, M))
cat("==========================================================\n\n")

# --- 4. 启动并行 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 环境预热：每个 Worker 加载 C++ 和必要函数
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    library(mvtnorm); library(Rcpp); library(RcppArmadillo); library(MTS); library(compiler)
    
    # 1. 加载新的 C++ 函数
    sourceCpp('statLMax.cpp')
    
    # 2. 加载数据生成
    if(file.exists("Gendata4PaperNew.R")) source("Gendata4PaperNew.R")
    
    # 3. 辅助函数
    lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1)) }
    lambda_compiled <- cmpfun(lambda)
    
    # Qk 需要在 Worker 中可用（或者通过 export 传递，这里定义更安全）
    Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
    Qk_compiled <- cmpfun(Qk)
})

# --- 5. 主循环 ---
startTime <- Sys.time()
final_results <- data.frame()

for (m_val in m_list) {
    for (L_val in L_list) {
        
        current_row <- list(Method = "Fisher", m = m_val, L = L_val)
        cat(sprintf("Computing m=%d, L=%d ... ", m_val, L_val))
        
        # 预计算权重列表 (避免在循环中重复计算)
        N_total <- (1 + L_val) * m_val
        k_vec   <- 1:(N_total - m_val)
        weight_list_per_lag <- list()
        
        # 预先计算好每个 lag 的 vaha 向量
        for(lag_v in lags) {
            Tm_lag <- m_val - lag_v
            # 注意: 这里在主进程计算一次即可，公式需与 C++ 逻辑一致
            # 为了确保精度，这里使用临时函数，Worker 中我们直接用 Qk_compiled
            w <- (Tm_lag + k_vec)^2 / (Tm_lag * Qk(k_vec/Tm_lag, gamma)^2)
            weight_list_per_lag[[as.character(lag_v)]] <- w
        }
        
        KT <- max(5, sqrt(log10(m_val)))
        sqrt_log_T <- 2 * sqrt(log10(m_val)/m_val)
        
        for (mod in model_list) {
            
            # 使用 foreach 并行
            decisions <- foreach(s = 1:nsim, .combine = c, .errorhandling = "remove",
                                 .export = c("lags", "nlag", "lag_weights", "CombiningMethod",
                                             "weight_list_per_lag", "alphaE", "kernal", "delta0",
                                             "M", "KT", "sqrt_log_T", "TypeI", "N_total", 
                                             "calc_combined_stats_vec")) %dopar% {
                                                 
                                                 # [A] 数据生成
                                                 yt <- tryCatch({
                                                     gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod)
                                                 }, error = function(e) return(NULL))
                                                 
                                                 if(is.null(yt)) return(NA)
                                                 if(!is.matrix(yt)) yt <- as.matrix(yt)
                                                 
                                                 # [B] 计算原始统计量 (T_orig) - 使用 statL_max
                                                 T_orig <- numeric(nlag)
                                                 for (j in 1:nlag) {
                                                     lag_v <- lags[j]
                                                     # 直接调用 statL_max，无需 max()
                                                     T_orig[j] <- statL_max(yt, m_val, lag_v, alphaE, 
                                                                            weight_list_per_lag[[as.character(lag_v)]], kernal, delta0)
                                                 }
                                                 
                                                 # [C] 自动 Block Size
                                                 blocksize <- max(1, round(m_val^(1/3))) # Fallback
                                                 try({
                                                     # 仅使用第一维进行 Block Size 估计
                                                     ut <- yt[1:m_val, 1] 
                                                     acf_out <- acf(ut, lag.max=2*KT, plot=FALSE, demean=TRUE)
                                                     R <- as.vector(acf_out$acf)
                                                     
                                                     tmp <- which(abs(R[1:KT]) < sqrt_log_T)
                                                     M_lag <- if (length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
                                                     
                                                     M_range <- -M_lag:M_lag
                                                     valid_idx <- which(abs(M_range) + 1 <= length(R))
                                                     
                                                     if(length(valid_idx) > 0) {
                                                         M_range <- M_range[valid_idx]
                                                         lam_vals <- lambda_compiled(M_range/M_lag)
                                                         ghat <- sum(lam_vals * R[abs(M_range) + 1])
                                                         Ghat <- sum(lam_vals * R[abs(M_range) + 1] * abs(M_range))
                                                         D.SB <- 2 * ghat^2
                                                         if(!is.na(D.SB) && D.SB > 1e-8) {
                                                             blocksize <- max(1, round((2 * Ghat^2 / D.SB)^(1/3) * m_val^(1/3)))
                                                         }
                                                     }
                                                 }, silent=TRUE)
                                                 
                                                 # [D] Full Bootstrap 循环 (耗时部分)
                                                 T_boot <- matrix(NA, nrow = M, ncol = nlag)
                                                 
                                                 # 预估需要的块数
                                                 num_blocks_est <- ceiling(N_total / blocksize) + 10
                                                 
                                                 for (b in 1:M) {
                                                     # 生成 SB 索引
                                                     len_vec <- rgeom(num_blocks_est, 1/blocksize) + 1
                                                     cumsum_len <- cumsum(len_vec)
                                                     cutoff <- which(cumsum_len >= N_total)[1]
                                                     # 保护措施
                                                     if(is.na(cutoff)) {
                                                         while(sum(len_vec) < N_total) len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
                                                         cumsum_len <- cumsum(len_vec)
                                                         cutoff <- which(cumsum_len >= N_total)[1]
                                                     }
                                                     len_vec <- len_vec[1:cutoff]
                                                     
                                                     start_points <- sample.int(m_val, length(len_vec), replace = TRUE)
                                                     
                                                     # 构建索引向量
                                                     idx_list <- lapply(seq_along(len_vec), function(k) {
                                                         (start_points[k] + 0:(len_vec[k]-1) - 1) %% m_val + 1
                                                     })
                                                     idx_vec <- unlist(idx_list)[1:N_total]
                                                     
                                                     y_star <- yt[idx_vec, , drop=FALSE]
                                                     
                                                     # 计算 Bootstrap 统计量 - 使用 statL_max
                                                     for (j in 1:nlag) {
                                                         lag_v <- lags[j]
                                                         val <- statL_max(y_star, m_val, lag_v, alphaE, 
                                                                          weight_list_per_lag[[as.character(lag_v)]], kernal, delta0)
                                                         T_boot[b, j] <- if(is.finite(val)) val else -Inf
                                                     }
                                                 }
                                                 
                                                 # [E] Fisher 组合检验逻辑
                                                 
                                                 # 1. 计算 P_obs
                                                 p_values_obs <- numeric(nlag)
                                                 for (j in 1:nlag) {
                                                     p_values_obs[j] <- (sum(T_boot[, j] >= T_orig[j]) + 1) / (M + 1)
                                                 }
                                                 
                                                 # 2. 计算 W_obs
                                                 W_vec_obs <- calc_combined_stats_vec(p_values_obs, CombiningMethod, lag_weights)
                                                 if(length(W_vec_obs) == 0 || is.na(W_vec_obs[1])) return(NA)
                                                 W_obs <- W_vec_obs[1]
                                                 
                                                 # 3. 构建 W 的 Null 分布 (Double Bootstrap Proxy)
                                                 # 利用 Bootstrap 样本在 Bootstrap 分布中的 Rank 来模拟 P 值
                                                 p_mat_boot <- matrix(0, nrow = M, ncol = nlag)
                                                 for (j in 1:nlag) {
                                                     # 将 T_orig 加入一起排序，确保逻辑一致性
                                                     all_vals <- c(T_orig[j], T_boot[, j])
                                                     ranks <- rank(-all_vals, ties.method = "max") 
                                                     # 取出对应 T_boot 的 rank (去掉第一个 T_orig)
                                                     p_mat_boot[, j] <- ranks[-1] / (M + 1)
                                                 }
                                                 
                                                 # 4. 计算 W_boot 分布
                                                 # apply 这里会返回一个向量 (因为 CombiningMethod 长度为 1)
                                                 W_boot_dist <- apply(p_mat_boot, 1, function(p_row) {
                                                     calc_combined_stats_vec(p_row, CombiningMethod, lag_weights)[1]
                                                 })
                                                 
                                                 # 5. 最终 P 值
                                                 p_final <- (sum(W_boot_dist >= W_obs) + 1) / (M + 1)
                                                 
                                                 as.numeric(p_final < TypeI)
                                             }
            
            # 结果聚合
            decisions <- as.numeric(decisions)
            decisions <- decisions[!is.na(decisions)]
            
            if(length(decisions) > 0) {
                rej_rate <- mean(decisions) * 100
            } else {
                rej_rate <- NA 
            }
            
            current_row[[paste0("N", mod)]] <- round(rej_rate, 1)
            cat(sprintf("N%d:%.1f%% ", mod, rej_rate))
        }
        cat("\n")
        final_results <- rbind(final_results, as.data.frame(current_row))
    }
}

stopCluster(cl)

# --- 6. 输出 ---
cat("\n========================================\n")
cat("Simulation Finished.\n")
cat("========================================\n")

col_order <- c("Method", "m", "L", paste0("N", model_list))
final_results <- final_results[, col_order]
print(final_results)

write.csv(final_results, "TypeI_Fisher_Standard_Max.csv", row.names = FALSE)
cat("Results saved to TypeI_Fisher_Standard_Max.csv\n")