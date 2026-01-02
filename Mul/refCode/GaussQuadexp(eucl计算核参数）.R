rm(list = ls())
######################### Multi-Model Type I Error (Gauss vs Quadexp) #########################

# --- 1. 环境设置 ---
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(compiler)
library(progress)
library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(foreach)
library(doParallel)

# 检查必要文件
if(!file.exists("Gendata4PaperNew.R")) stop("找不到 Gendata4PaperNew.R")

# --- 2. 模拟参数 ---
m_list <- c(100, 200)       # 训练样本大小
L_list <- c(1, 2, 3)        # 监控周期倍数
model_list <- 1:8           # 模型 N1 - N8
kernel_list <- c("gauss", "quadexp") # 遍历两种核

gamma <- 0
nsim <- 500
TypeI <- 0.05
lag <- 0
alphaE <- 1

# --- 3. 启动并行集群 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("启动并行计算，核心数: %d\n", num_cores))

# --- 4. 初始化节点环境 ---
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    library(mvtnorm); library(Rcpp); library(RcppArmadillo); library(MTS); library(compiler)
    
    # 加载 C++ 和 R 源码
    sourceCpp('statLMul.cpp')
    sourceCpp('eucLMulNew.cpp') # 必须加载，用于 quadexp 的 delta 计算
    source("Gendata4PaperNew.R")
    
    # 辅助函数
    Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
    lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1)) }
    
    lambda_compiled <- cmpfun(lambda)
    Qk_compiled     <- cmpfun(Qk)
    statL_compiled  <- cmpfun(statL) 
})

# --- 5. 主循环 ---
final_table <- data.frame()
total_scenarios <- length(kernel_list) * length(m_list) * length(L_list)
current_scenario <- 0

cat("==========================================================\n")
cat("Starting Simulation: Gauss (fixed) vs Quadexp (median heuristic)\n")
cat("==========================================================\n\n")

for (kern_name in kernel_list) {
    for (m_val in m_list) {
        for (L_val in L_list) {
            
            current_scenario <- current_scenario + 1
            cat(sprintf("[Scenario %d/%d] Kernel: %s | m=%d | L=%d... ", 
                        current_scenario, total_scenarios, kern_name, m_val, L_val))
            
            current_row <- list(Kernel = kern_name, m = m_val, L = L_val)
            
            for (mod in model_list) {
                
                # 预计算通用参数
                Tm <- m_val - lag
                N_total <- (1 + L_val) * m_val
                KT <- max(5, sqrt(log10(m_val)))
                log_T_limit <- 2 * sqrt(log10(m_val) / m_val)
                
                # === 并行模拟 ===
                sim_results <- foreach(s = 1:nsim, .combine = rbind, 
                                       .packages = c("mvtnorm", "Rcpp", "RcppArmadillo", "MTS"),
                                       .export = c("m_val", "L_val", "mod", "lag", "alphaE", "gamma", 
                                                   "kern_name", "KT", "log_T_limit", "N_total")) %dopar% {
                                                       
                                                       # 1. 生成数据
                                                       yt <- gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod)
                                                       if (!is.matrix(yt)) yt <- as.matrix(yt)
                                                       d <- ncol(yt)
                                                       
                                                       # 2. Block Size 选择 (保持不变)
                                                       Block_sizes <- numeric(d)
                                                       for(dim_i in 1:d) {
                                                           ut <- yt[1:m_val, dim_i]
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
                                                       
                                                       # 3. Bootstrap 索引
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
                                                           if(curr_pos + blk_len - 1 > N_total) blk_len <- N_total - curr_pos + 1
                                                           idx_vec[curr_pos:(curr_pos+blk_len-1)] <- (start_vec[j] + 0:(blk_len-1) - 1) %% m_val + 1
                                                           curr_pos <- curr_pos + blk_len
                                                           if(curr_pos > N_total) break
                                                       }
                                                       y.star <- yt[idx_vec, , drop = FALSE]
                                                       
                                                       # 4. 计算权重 vaha (节点内计算)
                                                       Tm_local <- m_val - lag
                                                       k_vec_local <- 1:(N_total - m_val)
                                                       vaha_local <- (Tm_local + k_vec_local)^2 / (Tm_local * Qk_compiled(k_vec_local/Tm_local, gamma)^2)
                                                       
                                                       # ============================================
                                                       # 5. 计算 delta0 (根据核函数类型)
                                                       # ============================================
                                                       curr_delta <- 1.0 # 默认为 1.0 (Gauss)
                                                       
                                                       if (kern_name == "quadexp") {
                                                           # 对 quadexp 使用训练集数据的中位数距离
                                                           # eucL 返回的是距离向量
                                                           E_vec <- eucL(yt[1:m_val, , drop=FALSE], lag)
                                                           valid_E <- E_vec[E_vec > 0]
                                                           
                                                           if (length(valid_E) > 0) {
                                                               curr_delta <- median(valid_E)
                                                           } else {
                                                               curr_delta <- 1.0 # 防止全0异常
                                                           }
                                                       }
                                                       
                                                       # 6. 计算统计量
                                                       Tn <- max(statL_compiled(yt, m_val, m = lag, alphaE, vaha_local, kern_name, curr_delta))
                                                       Tn_star <- max(statL_compiled(y.star, m_val, m = lag, alphaE, vaha_local, kern_name, curr_delta))
                                                       
                                                       c(Tn, Tn_star)
                                                   }
                
                # 聚合与计算
                Tn_vec <- sim_results[, 1]
                Tn_star_vec <- sim_results[, 2]
                crit_val <- quantile(Tn_star_vec, 1 - TypeI, na.rm=TRUE)
                rej_rate <- mean(Tn_vec > crit_val, na.rm=TRUE)
                
                current_row[[paste0("N", mod)]] <- round(rej_rate * 100, 1)
                cat(sprintf("N%d:%.1f%% ", mod, rej_rate*100))
            }
            cat("\n")
            final_table <- rbind(final_table, as.data.frame(current_row))
        }
    }
}

stopCluster(cl)

# --- 6. 结果保存 ---
colnames(final_table) <- c("Kernel", "m", "L", paste0("N", model_list))
print(final_table)

write.csv(final_table, "TypeI_Results_Gauss_vs_Quadexp.csv", row.names = FALSE)
cat("\n结果已保存至 'TypeI_Results_Gauss_vs_Quadexp.csv'。\n")