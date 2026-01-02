rm(list = ls())
######################### Multi-Model Comparison (Auto-Dimension + Kernel Loop) #########################

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

# 确保辅助文件存在
if(!file.exists("Gendata4PaperNew.R")){
    stop("错误: 找不到 Gendata4PaperNew.R，请确保文件在工作目录下。")
}

# --- 2. 模拟参数配置 ---
m_list <- c(100, 200)       # 训练样本大小 (m)
L_list <- c(1, 2, 3)        # 监控周期倍数 (L)
model_list <- 1:8           # 模型编号 N1 - N8
kernel_list <- c("gauss", "quadexp") # 【新增】遍历的核函数

# 固定参数
gamma <- 0
nsim <- 100                # 模拟次数
TypeI <- 0.05
lag <- 0
alphaE <- 1
fixed_delta <- 1.0          # 核函数参数统一固定为 1.0

# --- 3. 启动并行集群与环境初始化 ---

num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("已启动并行计算，使用核心数: %d\n", num_cores))

# 让每个 Worker 节点初始化环境
clusterEvalQ(cl, {
    # 注意：子节点的工作目录也需要设置，否则找不到 .cpp 文件
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    
    library(mvtnorm)
    library(Rcpp)
    library(RcppArmadillo)
    library(MTS)
    library(compiler)
    
    # 加载源码
    sourceCpp('statLMul.cpp')
    sourceCpp('eucLMulNew.cpp')
    source("Gendata4PaperNew.R")
    
    # 定义辅助函数
    Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
    lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1)) }
    
    # 编译函数
    lambda_compiled <- cmpfun(lambda)
    Qk_compiled     <- cmpfun(Qk)
    statL_compiled  <- cmpfun(statL) 
})

# --- 4. 主循环 ---

final_table <- data.frame()

cat("==========================================================\n")
cat("开始多模型 Type I Error 模拟 (Kernel: Gauss & Quadexp)\n")
cat(sprintf("Models: N1-N8 | m: %s | L: %s | Sim: %d\n", 
            paste(m_list, collapse=","), paste(L_list, collapse=","), nsim))
cat("==========================================================\n\n")

# 计算总任务数用于进度显示
total_scenarios <- length(kernel_list) * length(m_list) * length(L_list)
current_scenario <- 0

# 【新增】最外层遍历 Kernel
for (kern_name in kernel_list) {
    
    for (m_val in m_list) {
        for (L_val in L_list) {
            
            current_scenario <- current_scenario + 1
            cat(sprintf("[Scenario %d/%d] Kernel: %s | m=%d | L=%d... ", 
                        current_scenario, total_scenarios, kern_name, m_val, L_val))
            
            # 初始化当前行的结果
            current_row <- list(Kernel = kern_name, m = m_val, L = L_val)
            
            for (mod in model_list) {
                
                # 预计算参数
                Tm <- m_val - lag
                N_total <- (1 + L_val) * m_val
                k_vec <- 1:(N_total - m_val)
                
                KT <- max(5, sqrt(log10(m_val)))
                log_T_limit <- 2 * sqrt(log10(m_val) / m_val)
                
                # === 并行模拟 ===
                # 将当前的 kern_name 和 fixed_delta 传入
                sim_results <- foreach(s = 1:nsim, .combine = rbind, 
                                       .packages = c("mvtnorm", "Rcpp", "RcppArmadillo", "MTS"),
                                       .export = c("m_val", "L_val", "mod", "lag", "alphaE", "gamma", 
                                                   "kern_name", "fixed_delta", "KT", "log_T_limit", "N_total")) %dopar% {
                                                       
                                                       # 1. 生成数据
                                                       yt <- gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod)
                                                       if (!is.matrix(yt)) yt <- as.matrix(yt)
                                                       d <- ncol(yt)
                                                       
                                                       # 2. Block Size 选择
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
                                                       
                                                       # 3. Bootstrap 索引生成
                                                       expected_blocks <- ceiling((N_total + m_val) / blocksize) + 5
                                                       len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
                                                       while(sum(len_vec) < N_total) len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
                                                       len_vec <- len_vec[1:which(cumsum(len_vec) >= N_total)[1]]
                                                       
                                                       num_blocks <- length(len_vec)
                                                       start_vec <- sample.int(m_val, num_blocks, replace = TRUE)
                                                       
                                                       # 构建 idx_vec
                                                       idx_vec <- integer(N_total)
                                                       curr_pos <- 1
                                                       for(j in 1:num_blocks){
                                                           blk_len <- len_vec[j]
                                                           if(curr_pos + blk_len - 1 > N_total) blk_len <- N_total - curr_pos + 1
                                                           indices <- (start_vec[j] + 0:(blk_len-1) - 1) %% m_val + 1
                                                           idx_vec[curr_pos:(curr_pos+blk_len-1)] <- indices
                                                           curr_pos <- curr_pos + blk_len
                                                           if(curr_pos > N_total) break
                                                       }
                                                       y.star <- yt[idx_vec, , drop = FALSE]
                                                       
                                                       # 4. 计算权重 vaha (节点内计算)
                                                       Tm_local <- m_val - lag
                                                       k_vec_local <- 1:(N_total - m_val)
                                                       vaha_local <- (Tm_local + k_vec_local)^2 / (Tm_local * Qk_compiled(k_vec_local/Tm_local, gamma)^2)
                                                       
                                                       # 5. 计算统计量
                                                       # 【关键】：传入 kern_name 和 fixed_delta (1.0)
                                                       Tn <- max(statL_compiled(yt, m_val, m = lag, alphaE, vaha_local, kern_name, fixed_delta))
                                                       Tn_star <- max(statL_compiled(y.star, m_val, m = lag, alphaE, vaha_local, kern_name, fixed_delta))
                                                       
                                                       c(Tn, Tn_star)
                                                   }
                
                # 聚合结果
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

# --- 5. 结果输出 ---

cat("\n========================================\n")
cat("Simulation Finished.\n")
cat("========================================\n")

# 排序列：Kernel, m, L, N1...N8
col_order <- c("Kernel", "m", "L", paste0("N", model_list))
final_table <- final_table[, col_order]

print(final_table)

# 保存文件
file_name <- "TypeI_Results_Gauss_Quadexp.csv"
write.csv(final_table, file_name, row.names = FALSE)
cat(sprintf("\n结果已保存至 '%s'。\n", file_name))