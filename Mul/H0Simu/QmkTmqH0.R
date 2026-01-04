#########################################################################
# Multi-Model Type I Error (Optimized for Speed & Layout)
#########################################################################
rm(list = ls())
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

# --- 加载主进程依赖 ---
library(npcp)
library(mvtnorm)
library(MTS)
library(foreach)
library(doParallel)

# 检查文件
if(!file.exists("Gendata4PaperNew.R")){
    stop("Error: Gendata4PaperNew.R not found!")
}
# 主进程先加载一次，确保后续变量存在
source("Gendata4PaperNew.R") 

# --- 1. 参数设置 ---
trainSize     <- 100      
L_list        <- c(1, 2, 3) 
model_list    <- 1:8      

delta_fixed   <- 1.0      
gamma         <- 0        
alpha         <- 0.05    
nsim          <- 1000    # 恢复为 400 以获得稳定结果
B             <- 500    # 恢复为 200 (B太小会导致并行通信开销占比过大)
b_fixed       <- ceiling(trainSize^(1/3))        

cat("==========================================================\n")
cat("Start Simulation (High Performance Mode)\n")
cat("==========================================================\n\n")

# --- 2. 启动并行 (关键优化点) ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# [环境预热]：这是最高效的方法
# 在这里一次性把所有包和函数加载到工人节点
# 之后 foreach 就不需要再重复传输这些信息了
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    library(npcp); library(mvtnorm); library(MTS); library(Rcpp)
    source("Gendata4PaperNew.R")
})

# --- 3. 主循环 ---
startTime <- Sys.time()

results_Q <- data.frame()
results_T <- data.frame()

for (L_val in L_list) {
    
    current_row_Q <- list(m = trainSize, L = L_val)
    current_row_T <- list(m = trainSize, L = L_val)
    
    cat(sprintf("Computing L=%d (m=%d) ...\n", L_val, trainSize))
    
    for (mod in model_list) {
        
        # [优化]: 移除 .packages 和 .export 中关于 gendataPaper 的部分
        # 仅导出循环内的标量变量
        loop_results <- foreach(s = 1:nsim, .combine = rbind, 
                                .export = c("mod", "L_val", "trainSize", "delta_fixed", 
                                            "gamma", "B", "alpha", "b_fixed")) %dopar% {
                                                
                                                # 生成数据
                                                yt <- gendataPaper(monitoringPeriod = L_val, trainSize = trainSize, 
                                                                   model = mod, cptLocation = NULL, delta = 0)
                                                
                                                train_data   <- as.matrix(yt[1:trainSize, ])
                                                monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
                                                
                                                # npcp 核心计算
                                                traj <- npcp::simClosedEndCpDist(x.learn = train_data, n = nrow(yt), 
                                                                                 gamma = gamma, B = B, method = "mult", 
                                                                                 b = b_fixed)
                                                thresh <- npcp::threshClosedEndCpDist(traj, p = 1, alpha = alpha)
                                                
                                                det_stat <- npcp::detClosedEndCpDist(x.learn = train_data, x = monitor_data, 
                                                                                     gamma = gamma, delta = 1e-4)
                                                
                                                res_qmk <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mc", plot = FALSE)
                                                res_tmq <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mac", plot = FALSE)
                                                
                                                c(as.numeric(res_qmk$alarm), as.numeric(res_tmq$alarm))
                                            }
        
        rate_qmk <- mean(loop_results[, 1], na.rm = TRUE) * 100
        rate_tmq <- mean(loop_results[, 2], na.rm = TRUE) * 100
        
        col_name <- paste0("N", mod)
        current_row_Q[[col_name]] <- round(rate_qmk, 2)
        current_row_T[[col_name]] <- round(rate_tmq, 2)
        
        # 进度条效果
        cat(sprintf("  -> Model %d | Q: %5.2f%% | Tmq: %5.2f%%\n", mod, rate_qmk, rate_tmq))
    }
    
    results_Q <- rbind(results_Q, as.data.frame(current_row_Q))
    results_T <- rbind(results_T, as.data.frame(current_row_T))
}

stopCluster(cl)

# --- 4. 格式化输出 (保持你的改进版格式) ---
col_order <- c("m", "L", paste0("N", model_list))
results_Q <- results_Q[, col_order]
results_T <- results_T[, col_order]

cat("\n================ [Result Block: Q Statistic] ================\n")
print(results_Q)

cat("\n================ [Result Block: Tmq Statistic] ================\n")
print(results_T)

total_time <- difftime(Sys.time(), startTime, units="mins")
cat(sprintf("\nTotal Time: %.2f min\n", as.numeric(total_time)))

# --- 5. 保存 ---
df_out_Q <- cbind(Stat = "Qmk", results_Q)
df_out_T <- cbind(Stat = "Tmq", results_T)
df_final <- rbind(df_out_Q, df_out_T)

write.csv(df_final, "TypeI_Qmk_Tmq_100.csv", row.names = FALSE)
cat("\nResults saved to 'TypeI_Qmk_Tmq_100.csv'.\n")