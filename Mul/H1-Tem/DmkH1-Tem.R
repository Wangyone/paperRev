#########################################################################
# Structural Change (Model 301-306) Power & MDD Simulation
# Method: Wrap-Speed Bootstrap 
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


# --- 3. 参数配置 ---
lag_list    <- c(0, 1, 2)        # 滞后阶数
m_list      <- c(100, 200)       # 训练样本大小
L_list      <- c(1, 2, 3)        # 监控周期倍数
model_list  <- 301:306           # 结构变化模型
cpt_loc     <- 0.2               # 变点位置

nsim        <- 1000              # 模拟次数
TypeI       <- 0.05
kernal      <- "gauss"         # 核函数 gauss quadexp energy
alphaE      <- 1
gamma       <- 0.0

# --- 4. 初始化结果容器 (Base R Pre-allocation) ---

# 4.1 生成列名
# 第一列是 lag
col_names <- c("lag", "m", "L")
for(mod in model_list) {
    col_names <- c(col_names, paste0("N", mod, "-Power"), paste0("N", mod, "-MDD"))
}

# 4.2 计算总行数
total_rows <- length(lag_list) * length(m_list) * length(L_list)

# 4.3 创建空数据框 (预分配内存)
final_results_wide <- as.data.frame(matrix(NA, nrow = total_rows, ncol = length(col_names)))
colnames(final_results_wide) <- col_names

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

# 导出通用变量
clusterExport(cl, c("alphaE", "kernal", "gamma", "cpt_loc"))

# --- 7. 主循环 ---
cat("==========================================================\n")
cat("Start Structural Change Simulation (Multi-Lag Wrap-Speed)\n")
cat(sprintf("Models: 301-306 | Lag: %s | Nsim: %d\n", paste(lag_list, collapse=","), nsim))
cat("==========================================================\n\n")

row_idx <- 0 

# [Loop 1] 遍历 Lag
for (lag_val in lag_list) {
    
    cat(sprintf("\n>>> Processing Lag = %d\n", lag_val))
    
    # [Loop 2] 遍历 m
    for (m_val in m_list) {
        
        # [Loop 3] 遍历 L
        for (L_val in L_list) {
            
            row_idx <- row_idx + 1 
            
            # 填充基础列
            final_results_wide[row_idx, "lag"] <- lag_val
            final_results_wide[row_idx, "m"]   <- m_val
            final_results_wide[row_idx, "L"]   <- L_val
            
            # 参数准备
            Tm <- m_val - lag_val
            N_total <- (1 + L_val) * m_val
            true_cp_time <- m_val + floor(L_val * m_val * cpt_loc)
            
            cat(sprintf("[Row %d/%d] lag=%d, m=%d, L=%d ... ", row_idx, total_rows, lag_val, m_val, L_val))
            
            # [Loop 4] 遍历模型
            for (mod_id in model_list) {
                
                # === 并行 Wrap-Speed 模拟 ===
                results_list <- foreach(s = 1:nsim, 
                                        .packages=c('mvtnorm','Rcpp','RcppArmadillo'),
                                        .export=c('m_val', 'L_val', 'mod_id', 'Tm', 'N_total', 'lag_val')) %dopar% {
                                            
                                            # 子节点计算权重
                                            k_vec <- 1:(N_total - m_val)
                                            vaha_vec <- (Tm + k_vec)^2 / (Tm * Qk_compiled(k_vec/Tm, 0.0)^2)
                                            
                                            # 1. 生成数据
                                            yt <- gendataPaper(monitoringPeriod = L_val, trainSize = m_val, 
                                                               model = mod_id, cptLocation = cpt_loc)
                                            if(!is.matrix(yt)) yt <- as.matrix(yt)
                                            
                                            d_dim <- ncol(yt)
                                            delta_dynamic <- 1.0 * d_dim 
                                            
                                            # 2. 计算观测统计量 (传入 lag_val)
                                            stats_obs <- statL(yt, m_val, lag_val, alphaE, vaha_vec, kernal, delta_dynamic)
                                            
                                            # 3. Bootstrap
                                            yt_train <- yt[1:m_val, , drop=FALSE]
                                            
                                            # (A) Block Size
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
                                            
                                            # (B) Stationary Bootstrap
                                            expected_blocks <- ceiling((N_total + m_val) / blocksize) + 5
                                            len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
                                            while(sum(len_vec) < N_total) len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
                                            len_vec <- len_vec[1:which(cumsum(len_vec) >= N_total)[1]]
                                            
                                            num_blocks <- length(len_vec)
                                            start_vec <- sample.int(m_val, num_blocks, replace = TRUE)
                                            
                                            idx_vec <- integer(N_total)
                                            curr_pos <- 1
                                            for(j in 1:num_blocks){
                                                blk_len <- len_vec[j]
                                                end_pos <- curr_pos + blk_len - 1
                                                if(end_pos > N_total) blk_len <- N_total - curr_pos + 1
                                                indices <- (start_vec[j] + 0:(blk_len-1) - 1) %% m_val + 1
                                                idx_vec[curr_pos:(curr_pos+blk_len-1)] <- indices
                                                curr_pos <- curr_pos + blk_len
                                                if(curr_pos > N_total) break
                                            }
                                            y_star <- yt_train[idx_vec, , drop = FALSE]
                                            
                                            # 4. Bootstrap Max (传入 lag_val)
                                            stats_boot <- statL(y_star, m_val, lag_val, alphaE, vaha_vec, kernal, delta_dynamic)
                                            max_boot <- max(stats_boot)
                                            
                                            list(obs_vec = stats_obs, boot_max = max_boot)
                                        }
                
                # === 汇总计算 ===
                all_boot_maxs <- sapply(results_list, function(x) x$boot_max)
                threshold <- quantile(all_boot_maxs, 1 - TypeI, na.rm=TRUE)
                
                eval_metrics <- t(sapply(results_list, function(res) {
                    obs_stats <- res$obs_vec
                    alarm_indices <- which(obs_stats > threshold)
                    
                    if (length(alarm_indices) > 0) {
                        first_k <- alarm_indices[1]
                        alarm_time <- m_val + first_k
                        delay <- alarm_time - true_cp_time
                        valid_delay <- if(delay > 0) delay else NA
                        return(c(1, valid_delay))
                    } else {
                        return(c(0, NA))
                    }
                }))
                
                power_val <- mean(eval_metrics[, 1]) * 100
                valid_delays <- eval_metrics[, 2]
                mdd_val <- mean(valid_delays[!is.na(valid_delays)], na.rm=TRUE)
                
                # === 填表 ===
                col_pow <- paste0("N", mod_id, "-Power")
                col_mdd <- paste0("N", mod_id, "-MDD")
                
                final_results_wide[row_idx, col_pow] <- round(power_val, 2)
                final_results_wide[row_idx, col_mdd] <- round(mdd_val, 2)
                
                cat(sprintf("N%d ", mod_id))
            } # end models loop
            
            cat("\n")
            # 实时保存（防止崩溃）
            # write.csv(final_results_wide, "temp_results.csv", row.names=FALSE)
            
        } # end L loop
    } # end m loop
} # end lag loop

stopCluster(cl)

# --- 8. 结果保存与展示 ---
cat("\n========================================\n")
cat("Simulation Finished.\n")
cat("========================================\n")

# 输出预览（前20行）
print(head(final_results_wide, 20))

fileSaveName <- paste0("Temporal_Power_MDD_", kernal)
write.csv(final_results_wide, file = paste0(fileSaveName, ".csv"), row.names = FALSE)
cat(sprintf("\nResult saved to %s.csv\n", fileSaveName))