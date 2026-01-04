rm(list = ls())
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(foreach)
library(doParallel)

# 检查文件
if(!file.exists("Gendata4PaperNew.R")) stop("Gendata4PaperNew.R missing")
if(!file.exists("statLMax.cpp")) stop("statLMax.cpp missing")
# 确保你创建了新的辅助文件
if(!file.exists("boot_utils.cpp")) stop("Please create boot_utils.cpp first!")

# --- 参数设置 ---
m_list      <- c(100, 200)      
L_list      <- c(1, 2, 3)        
model_list  <- 1:8              
CombiningMethod <- c(1) # Fisher
lags            <- c(0, 1, 2)  
nlag            <- length(lags)
lag_weights     <- c(1.0, 1.0, 1.0) 

gamma    <- 0
nsim     <- 1000      # 正式模拟次数
M        <- 500       # Bootstrap 次数
TypeI    <- 0.05
kernal   <- "gauss"
delta0   <- 1.0
alphaE   <- 1

# 辅助函数
Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1)) }
lambda_compiled <- cmpfun(lambda)

cat("Starting Accelerated Simulation (R DataGen + Cpp Stat + Cpp Utils)\n")

# --- 启动并行 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# --- Worker 环境初始化 ---
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    library(mvtnorm); library(Rcpp); library(RcppArmadillo); library(MTS); library(compiler)
    
    # 1. 加载核心统计量函数 (不动)
    sourceCpp('statLMax.cpp')
    
    # 2. 加载加速辅助函数 (新增优化)
    sourceCpp('boot_utils.cpp') 
    
    # 3. 加载数据生成 (不动)
    source("Gendata4PaperNew.R")
    
    # 4. 本地辅助函数
    Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
    lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1)) }
    lambda_compiled <- cmpfun(lambda)
})

# --- 主循环 ---
final_results <- data.frame()

for (m_val in m_list) {
    for (L_val in L_list) {
        
        current_row <- list(Method = "Fisher", m = m_val, L = L_val)
        cat(sprintf("Computing m=%d, L=%d ... ", m_val, L_val))
        
        # 预计算参数
        N_total <- (1 + L_val) * m_val
        k_vec   <- 1:(N_total - m_val)
        
        weight_list_per_lag <- list()
        for(lag_v in lags) {
            Tm_lag <- m_val - lag_v
            w <- (Tm_lag + k_vec)^2 / (Tm_lag * Qk(k_vec/Tm_lag, gamma)^2)
            weight_list_per_lag[[as.character(lag_v)]] <- w
        }
        
        KT <- max(5, sqrt(log10(m_val)))
        sqrt_log_T <- 2 * sqrt(log10(m_val)/m_val)
        
        for (mod in model_list) {
            
            # .noexport 列表很重要，防止覆盖 Worker 的编译函数
            decisions <- foreach(s = 1:nsim, .combine = c, .errorhandling = "remove",
                                 .noexport = c("statL_max", "gendataPaper", "gen_sb_indices", "calc_fisher_mat"),
                                 .export = c("lags", "nlag", "lag_weights", 
                                             "weight_list_per_lag", "alphaE", "kernal", "delta0",
                                             "M", "KT", "sqrt_log_T", "TypeI", "N_total")) %dopar% {
                                                 
                                                 # 1. 数据生成 (R)
                                                 yt <- tryCatch({
                                                     gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod)
                                                 }, error = function(e) return(NULL))
                                                 
                                                 if(is.null(yt)) return(NA)
                                                 if(!is.matrix(yt)) yt <- as.matrix(yt)
                                                 
                                                 # 2. 计算 T_orig (C++)
                                                 T_orig <- numeric(nlag)
                                                 for (j in 1:nlag) {
                                                     lag_v <- lags[j]
                                                     T_orig[j] <- statL_max(yt, m_val, lag_v, alphaE, 
                                                                            weight_list_per_lag[[as.character(lag_v)]], kernal, delta0)
                                                 }
                                                 
                                                 # 3. 自动 Block Size (R)
                                                 blocksize <- max(1, round(m_val^(1/3))) 
                                                 try({
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
                                                 
                                                 # 4. Full Bootstrap 循环 (优化版)
                                                 T_boot <- matrix(NA, nrow = M, ncol = nlag)
                                                 
                                                 for (b in 1:M) {
                                                     # [优化A] C++ 极速生成索引 (替代 R 的 rgeom 等)
                                                     idx_vec <- gen_sb_indices(N_total, blocksize)
                                                     
                                                     # R 切片 (这一步现在是主要瓶颈，但在不重写 statL_max 的前提下无法消除)
                                                     y_star <- yt[idx_vec, , drop=FALSE]
                                                     
                                                     # C++ 计算统计量
                                                     for (j in 1:nlag) {
                                                         lag_v <- lags[j]
                                                         val <- statL_max(y_star, m_val, lag_v, alphaE, 
                                                                          weight_list_per_lag[[as.character(lag_v)]], kernal, delta0)
                                                         T_boot[b, j] <- if(is.finite(val)) val else -Inf
                                                     }
                                                 }
                                                 
                                                 # 5. Fisher 组合检验 (优化版)
                                                 
                                                 # (1) P_obs
                                                 p_values_obs <- numeric(nlag)
                                                 for (j in 1:nlag) {
                                                     p_values_obs[j] <- (sum(T_boot[, j] >= T_orig[j]) + 1) / (M + 1)
                                                 }
                                                 
                                                 # [优化B] C++ 极速 Fisher 计算
                                                 # 将向量转为 1行矩阵传入
                                                 W_obs <- calc_fisher_mat(matrix(p_values_obs, nrow=1), lag_weights)[1]
                                                 if(is.na(W_obs)) return(NA)
                                                 
                                                 # (2) W_boot Null 分布
                                                 p_mat_boot <- matrix(0, nrow = M, ncol = nlag)
                                                 for (j in 1:nlag) {
                                                     all_vals <- c(T_orig[j], T_boot[, j])
                                                     ranks <- rank(-all_vals, ties.method = "max") 
                                                     p_mat_boot[, j] <- ranks[-1] / (M + 1)
                                                 }
                                                 
                                                 # [优化B] C++ 极速矩阵计算 (替代 apply)
                                                 W_boot_dist <- calc_fisher_mat(p_mat_boot, lag_weights)
                                                 
                                                 # (3) 判决
                                                 p_final <- (sum(W_boot_dist >= W_obs) + 1) / (M + 1)
                                                 as.numeric(p_final < TypeI)
                                             }
            
            decisions <- as.numeric(decisions[!is.na(decisions)])
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
print(final_results)
write.csv(final_results, "Fast_FullBoot_Results.csv", row.names = FALSE)