#########################################################################
# Multi-Lag Combination Test Simulation (Full Bootstrap)
# Methods: Fisher (1), Tippett (2), Stouffer (3)
# Complexity: O(Nsim * B) - High Computational Load
# Kernel: QuadExp | Bandwidth: 1.0 * d
#########################################################################

rm(list = ls())
# 请根据实际路径修改
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)

# --- 1. 检查依赖 ---
if(!file.exists("Gendata4PaperNew.R")) stop("Gendata4PaperNew.R missing")
if(!file.exists("statLMul.cpp")) stop("statLMul.cpp missing") 

# --- 2. 核心组合统计量函数 (支持向量化输出) ---
# 1=Fisher, 2=Tippett, 3=Stouffer
calc_combined_stats_map <- function(pvals, methods_indices, weights) {
    # pvals: numeric vector of p-values for different lags
    if(any(is.na(pvals))) return(rep(NA, length(methods_indices)))
    
    # 预处理
    pvals <- pmax(pvals, 1e-10)
    pvals <- pmin(pvals, 1 - 1e-10)
    
    res_list <- list()
    
    # Method 1: Fisher
    if(1 %in% methods_indices) {
        fisher_stat <- -2 * sum(weights * log(pvals))
        res_list[["Fisher"]] <- fisher_stat
    }
    # Method 2: Tippett
    if(2 %in% methods_indices) {
        tippett_stat <- 1 - min(pvals)
        res_list[["Tippett"]] <- tippett_stat
    }
    # Method 3: Stouffer
    if(3 %in% methods_indices) {
        z_scores <- qnorm(1 - pvals)
        # 归一化权重
        w_norm <- weights / sqrt(sum(weights^2))
        stouffer_stat <- sum(weights * z_scores) / sqrt(sum(weights^2))
        res_list[["Stouffer"]] <- stouffer_stat
    }
    
    return(unlist(res_list))
}

# --- 3. 参数配置 ---
m_list      <- c(100, 200)       # 训练样本大小
L_list      <- c(1, 2, 3)        # 监控周期倍数
model_list  <- 1:8               # H0 模型 (Type I Error)

# 组合检验参数
lags            <- c(0, 1, 2)    
nlag            <- length(lags)
lag_weights     <- rep(1.0, nlag) 
CombiningMethod <- c(1, 2, 3)     # 同时计算 1=Fisher, 2=Tippett, 3=Stouffer

# 模拟参数
nsim        <- 1000               # 外层模拟次数 (注意：Full Boot 很慢，测试时设小一点)
B           <- 500               # 内层 Bootstrap 次数
TypeI       <- 0.05
kernal      <- "quadexp"         
alphaE      <- 1
gamma       <- 0.0
delta0      <- 1.0               

# --- 4. 辅助：生成结果列名 ---
method_names <- c("Fisher", "Tippett", "Stouffer")
active_methods <- method_names[CombiningMethod]

# --- 5. 启动并行集群 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
if(num_cores < 1) num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# --- 6. 环境推送到子节点 ---
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')
    library(mvtnorm); library(Rcpp); library(RcppArmadillo); library(compiler)
    
    source("Gendata4PaperNew.R")
    sourceCpp("statLMul.cpp") 
    
    Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
    lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1)) }
    
    Qk_compiled     <- cmpfun(Qk)
    lambda_compiled <- cmpfun(lambda)
})

# 导出变量
clusterExport(cl, c("alphaE", "kernal", "gamma", "lags", "nlag", 
                    "lag_weights", "CombiningMethod", "calc_combined_stats_map", "B", "TypeI"))

# --- 7. 主循环 ---
cat("==========================================================\n")
cat("Start Full Bootstrap Combination Test\n")
cat(sprintf("Methods: %s | Nsim: %d | Boot(B): %d\n", paste(active_methods, collapse=","), nsim, B))
cat("==========================================================\n\n")

final_results <- data.frame()

for (m_val in m_list) {
    for (L_val in L_list) {
        
        Tm_base <- m_val 
        N_total <- (1 + L_val) * m_val
        
        # 初始化当前行数据 (包含 Method 标识)
        current_rows <- list()
        for(met in active_methods) {
            current_rows[[met]] <- list(Method = met, m = m_val, L = L_val)
        }
        
        cat(sprintf("[Scenario] m=%d, L=%d ... ", m_val, L_val))
        
        for (mod_id in model_list) {
            
            # === 并行 Full Bootstrap 模拟 ===
            # 每个核心负责 s=1..nsim 中的一次完整模拟（含B次Bootstrap）
            # 返回值：一个矩阵，行=nsim，列=length(active_methods)
            # 每一列代表该次模拟中，对应方法是否拒绝 H0 (0/1)
            
            decisions_mat <- foreach(s = 1:nsim, .combine = rbind, 
                                     .packages=c('mvtnorm','Rcpp','RcppArmadillo')) %dopar% {
                                         
                                         # --- 1. 生成观测数据 ---
                                         yt <- tryCatch({
                                             gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod_id)
                                         }, error = function(e) return(NULL))
                                         if(is.null(yt)) return(rep(NA, length(active_methods)))
                                         if(!is.matrix(yt)) yt <- as.matrix(yt)
                                         
                                         d_dim <- ncol(yt)
                                         delta_dyn <- 1.0 * d_dim
                                         
                                         # --- 2. 计算观测统计量 T_obs (向量, 长度 nlag) ---
                                         T_obs <- numeric(nlag)
                                         for(j in 1:nlag) {
                                             curr_lag <- lags[j]
                                             Tm <- m_val - curr_lag
                                             k_vec <- 1:(N_total - m_val)
                                             vaha_vec <- (Tm + k_vec)^2 / (Tm * Qk_compiled(k_vec/Tm, 0.0)^2)
                                             T_obs[j] <- max(statL(yt, m_val, curr_lag, alphaE, vaha_vec, kernal, delta_dyn))
                                         }
                                         
                                         # --- 3. Full Bootstrap 循环 ---
                                         # 生成 B 个 Bootstrap 样本，计算 T_boot 矩阵 (B x nlag)
                                         yt_train <- yt[1:m_val, , drop=FALSE]
                                         
                                         # Block Size 选择
                                         KT <- max(5, sqrt(log10(m_val)))
                                         log_T_limit <- 2 * sqrt(log10(m_val) / m_val)
                                         Block_sizes <- numeric(d_dim)
                                         for(dim_i in 1:d_dim) {
                                             ut <- yt_train[, dim_i]
                                             acf_res <- acf(ut, lag.max = 2*KT, plot = FALSE)$acf
                                             R <- as.vector(acf_res)
                                             tmp <- which(abs(R[1:KT]) < log_T_limit)
                                             M <- if (length(tmp) > 0) 2 * (tmp[1] - 1) else 2 * (KT - 1)
                                             idx <- -M:M
                                             w <- lambda_compiled(idx / M)
                                             ghat <- sum(w * R[abs(idx) + 1])
                                             Ghat <- sum(w * R[abs(idx) + 1] * abs(idx))
                                             D.SB <- 2 * ghat^2
                                             Block_sizes[dim_i] <- if(D.SB < 1e-10) 1 else max(1, round((2 * Ghat^2 / D.SB)^(1/3) * m_val^(1/3)))
                                         }
                                         blocksize <- max(1, round(mean(Block_sizes)))
                                         
                                         T_boot_mat <- matrix(NA, nrow = B, ncol = nlag)
                                         
                                         # 开始 B 次重采样
                                         for(b in 1:B) {
                                             # 生成 y_star
                                             expected_blocks <- ceiling(N_total / blocksize) + 5
                                             len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
                                             while(sum(len_vec) < N_total) len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
                                             len_vec <- len_vec[1:which(cumsum(len_vec) >= N_total)[1]]
                                             num_blocks <- length(len_vec)
                                             start_vec <- sample.int(m_val, num_blocks, replace = TRUE)
                                             idx_vec <- integer(N_total)
                                             curr_pos <- 1
                                             for(k in 1:num_blocks){
                                                 blk_len <- len_vec[k]
                                                 end_pos <- curr_pos + blk_len - 1
                                                 if(end_pos > N_total) blk_len <- N_total - curr_pos + 1
                                                 indices <- (start_vec[k] + 0:(blk_len-1) - 1) %% m_val + 1
                                                 idx_vec[curr_pos:(curr_pos+blk_len-1)] <- indices
                                                 curr_pos <- curr_pos + blk_len
                                                 if(curr_pos > N_total) break
                                             }
                                             y_star <- yt_train[idx_vec, , drop = FALSE]
                                             
                                             # 计算 stats
                                             for(j in 1:nlag) {
                                                 curr_lag <- lags[j]
                                                 Tm <- m_val - curr_lag
                                                 k_vec <- 1:(N_total - m_val)
                                                 vaha_vec <- (Tm + k_vec)^2 / (Tm * Qk_compiled(k_vec/Tm, 0.0)^2)
                                                 T_boot_mat[b, j] <- max(statL(y_star, m_val, curr_lag, alphaE, vaha_vec, kernal, delta_dyn))
                                             }
                                         }
                                         
                                         # --- 4. 计算 P 值 (Local) ---
                                         # p_obs[j]: T_obs[j] 在 T_boot_mat[, j] 中的位置
                                         p_vals_obs <- numeric(nlag)
                                         for(j in 1:nlag) {
                                             p_vals_obs[j] <- (sum(T_boot_mat[, j] >= T_obs[j]) + 1) / (B + 1)
                                         }
                                         
                                         # 计算观测数据的组合统计量 (Vector: Fisher, Tippett, Stouffer)
                                         W_obs_vec <- calc_combined_stats_map(p_vals_obs, CombiningMethod, lag_weights)
                                         
                                         # --- 5. Double Bootstrap Proxy (Local Threshold) ---
                                         # 计算 Bootstrap 样本自身的 P 值矩阵 (B x nlag)
                                         P_boot_mat <- matrix(NA, nrow = B, ncol = nlag)
                                         for(j in 1:nlag) {
                                             # rank(-x) gives 1 for largest. We want P-value = (count >= val) / (B+1)
                                             # This is equivalent to rank within descending sort
                                             P_boot_mat[, j] <- rank(-T_boot_mat[, j], ties.method = "max") / (B + 1)
                                         }
                                         
                                         # 计算 B 个组合统计量 (Null Distribution)
                                         # Apply across rows: returns Matrix (methods x B)
                                         W_boot_dist_mat <- apply(P_boot_mat, 1, function(p_row) {
                                             calc_combined_stats_map(p_row, CombiningMethod, lag_weights)
                                         })
                                         # 转置为 (B x methods)
                                         if(length(active_methods) > 1) {
                                             W_boot_dist_mat <- t(W_boot_dist_mat)
                                         } else {
                                             W_boot_dist_mat <- matrix(W_boot_dist_mat, ncol=1)
                                         }
                                         
                                         # --- 6. 判决 ---
                                         # 对每种方法分别判断
                                         decisions <- numeric(length(active_methods))
                                         for(k in 1:length(active_methods)) {
                                             # 当前方法的 Null 分布
                                             null_dist <- W_boot_dist_mat[, k]
                                             # 当前方法的观测值
                                             obs_val <- W_obs_vec[k]
                                             # 临界值 (95%)
                                             # Fisher/Stouffer: Larger is significant
                                             # Tippett: Calculation above defined it as (1-min P), so Larger is significant
                                             crit_val <- quantile(null_dist, 1 - TypeI, na.rm=TRUE)
                                             
                                             if(obs_val > crit_val) decisions[k] <- 1 else decisions[k] <- 0
                                         }
                                         
                                         return(decisions)
                                     }
            
            # === 汇总拒绝率 ===
            # decisions_mat: nsim x n_methods
            if(nsim == 1) decisions_mat <- t(decisions_mat)
            
            for(k in 1:length(active_methods)) {
                met_name <- active_methods[k]
                rej_rate <- mean(decisions_mat[, k], na.rm=TRUE) * 100
                current_rows[[met_name]][[paste0("N", mod_id)]] <- round(rej_rate, 2)
            }
            
            cat(sprintf("N%d Done ", mod_id))
            
        } # end models
        cat("\n")
        
        # 将结果添加到总表
        for(met in active_methods) {
            final_results <- rbind(final_results, as.data.frame(current_rows[[met]]))
        }
        
    } # end L
} # end m

stopCluster(cl)

# --- 8. 输出结果 ---
cat("\n========================================\n")
cat("Simulation Finished.\n")
print(final_results)

# 1. 设置列的顺序 (Method, m, L, N1...N8)
col_order <- c("Method", "m", "L", paste0("N", model_list))
final_results <- final_results[, col_order]
# 2. 排序 (Sort)
# 逻辑：先按 Method 分组 (Fisher -> Tippett -> Stouffer)，组内按 m 升序，再按 L 升序
# 技巧：将 Method 转为 factor 并指定 level 顺序，以此控制排序先后
final_results$Method <- factor(final_results$Method, levels = c("Fisher", "Tippett", "Stouffer"))
final_results <- final_results[order(final_results$Method, final_results$m, final_results$L), ]

# 3. 打印预览
print(final_results)


fileSaveName <- paste0("TypeI_FullBoot_Comb_statL_", kernal)
write.csv(final_results, file = paste0(fileSaveName, ".csv"), row.names = FALSE)
cat(sprintf("\nResult saved to %s.csv\n", fileSaveName))
