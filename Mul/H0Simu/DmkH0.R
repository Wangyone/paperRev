rm(list = ls())
######################### Multi-Model Type I Error Comparison (Fixed) #########################

# --- 1. 环境设置 ---
# 请根据实际情况修改路径
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

# --- 2. 准备数据生成函数 ---
if(!file.exists("Gendata4PaperNew.R")){
    stop("错误: 找不到 Gendata4PaperNew.R，请确保文件在工作目录下。")
}

# --- 3. 模拟参数配置 ---
m_list <- c(100, 200)       # 训练样本大小 (m)
L_list <- c(1, 2, 3)        # 监控周期倍数 (L)
model_list <- 1:8           # 模型编号 N1 - N8

gamma <- 0
nsim <- 1000                # 模拟次数
TypeI <- 0.05
kernal <- "quadexp"           # 核函数 gauss quadexp energy
lag <- 0
alphaE <- 1

# --- 4. 并行集群与环境初始化 ---

num_cores <- parallel::detectCores(logical = FALSE) - 1
# 如果核心数不够，防止报错
if(num_cores < 1) num_cores <- 1 

cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("已启动并行计算，使用核心数: %d\n", num_cores))

# 预加载包和函数到各个节点
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    library(mvtnorm)
    library(Rcpp)
    library(RcppArmadillo)
    library(MTS)
    library(compiler)
    
    # 加载 C++ 和 R 源码
    sourceCpp('statLMul.cpp')
    # sourceCpp('eucLMulNew.cpp') # 如果不需要欧氏距离对比，可以注释掉以节省时间
    source("Gendata4PaperNew.R")
    
    # 辅助函数
    Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
    lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1)) }
    
    # 编译优化
    lambda_compiled <- cmpfun(lambda)
    Qk_compiled     <- cmpfun(Qk)
    statL_compiled  <- cmpfun(statL) 
})

# --- 5. 主循环 ---

final_table <- data.frame()

cat("==========================================================\n")
cat("开始多模型 Type I Error 模拟 (Bandwidth = 1.0 * d)\n")
cat(sprintf("Models: N1-N8 | m: %s | L: %s | Sim: %d\n", 
            paste(m_list, collapse=","), paste(L_list, collapse=","), nsim))
cat("==========================================================\n\n")

total_scenarios <- length(m_list) * length(L_list)
current_scenario <- 0

for (m_val in m_list) {
    for (L_val in L_list) {
        
        current_scenario <- current_scenario + 1
        cat(sprintf("[Scenario %d/%d] Computing m=%d, L=%d... ", 
                    current_scenario, total_scenarios, m_val, L_val))
        
        current_row <- list(m = m_val, L = L_val)
        
        for (mod in model_list) {
            
            # 计算参数
            Tm <- m_val - lag
            N_total <- (1 + L_val) * m_val
            
            # Block Bootstrap 参数
            KT <- max(5, sqrt(log10(m_val)))
            log_T_limit <- 2 * sqrt(log10(m_val) / m_val)
            
            # === 并行模拟 ===
            sim_results <- foreach(s = 1:nsim, .combine = rbind, 
                                   .packages = c("mvtnorm", "Rcpp", "RcppArmadillo", "MTS"),
                                   .export = c("m_val", "L_val", "mod", "lag", "alphaE", "gamma", 
                                               "kernal", "KT", "log_T_limit", "N_total")) %dopar% {
                                                   
                                                   # 1. 生成数据
                                                   yt <- gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod)
                                                   if (!is.matrix(yt)) yt <- as.matrix(yt)
                                                   d <- ncol(yt) # 获取维度 d
                                                   
                                                   # 2. Block Size 选择 (自动最优带宽)
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
                                                   
                                                   # 3. Stationary Bootstrap 生成 y.star
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
                                                   y.star <- yt[idx_vec, , drop = FALSE]
                                                   
                                                   # 4. 权重向量
                                                   Tm_local <- m_val - lag
                                                   k_vec_local <- 1:(N_total - m_val)
                                                   vaha_local <- (Tm_local + k_vec_local)^2 / (Tm_local * Qk_compiled(k_vec_local/Tm_local, gamma)^2)
                                                   
                                                   # 5. 计算统计量 [关键修改处]
                                                   # 将 delta0 参数设置为 1.0 * d
                                                   delta_param <- 1.0 * d 
                                                   
                                                   Tn <- max(statL_compiled(yt, m_val, m = lag, alphaE, vaha_local, kernal, delta_param))
                                                   Tn_star <- max(statL_compiled(y.star, m_val, m = lag, alphaE, vaha_local, kernal, delta_param))
                                                   
                                                   c(Tn, Tn_star)
                                               }
            
            # 聚合结果
            Tn_vec <- sim_results[, 1]
            Tn_star_vec <- sim_results[, 2]
            
            # 计算 Type I Error
            crit_val <- quantile(Tn_star_vec, 1 - TypeI, na.rm=TRUE)
            rej_rate <- mean(Tn_vec > crit_val, na.rm=TRUE)
            
            current_row[[paste0("N", mod)]] <- round(rej_rate * 100, 1)
            cat(sprintf("N%d:%.1f%% ", mod, rej_rate*100))
        }
        cat("\n")
        final_table <- rbind(final_table, as.data.frame(current_row))
    }
}

stopCluster(cl)

# --- 6. 结果输出 ---
cat("\n========================================\n")
cat("Simulation Finished.\n")
cat("========================================\n")

col_order <- c("m", "L", paste0("N", model_list))
final_table <- final_table[, col_order]

print(final_table)

# 文件名包含 kernal 以区分
fileSaveName <- paste0("TypeI_Results_Dmk_", kernal)
write.csv(final_table, file = paste0(fileSaveName, ".csv"), row.names = FALSE)

cat(sprintf("\n结果已保存至 '%s.csv'。\n", fileSaveName))
