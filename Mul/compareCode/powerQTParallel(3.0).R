#########################################################################
#  Multi-Model Type I Error (Split Output Format)
#########################################################################
rm(list = ls())
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')
library(npcp)
library(mvtnorm)
library(MTS)
library(foreach)
library(doParallel)
source("D:/Desktop/paperRev/JointCUSUM/Mul/Gendata4PaperNew.R", echo=TRUE)


# --- 1. 参数设置 ---
trainSize        <- 100     
L_list           <- c(1, 2, 3) 
model_list       <- 1:8     

delta_fixed      <- 1.0     
gamma            <- 0       
alpha            <- 0.05    
nsim             <- 400    
B                <- 200     
b_fixed          <- ceiling(trainSize^(1/3))        

cat("==========================================================\n")
cat("开始多模型并行模拟 (Formatted Output)\n")
cat("==========================================================\n\n")

# --- 2. 启动并行 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 环境预热
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    library(npcp); library(mvtnorm); library(MTS); library(Rcpp)
    source("Gendata4PaperNew.R")
})

# --- 3. 主循环 ---
startTime <- Sys.time()

# 初始化两个独立的数据框
results_Q <- data.frame()
results_T <- data.frame()

for (L_val in L_list) {
    
    # 初始化当前行的基础信息
    current_row_Q <- list(m = trainSize, L = L_val)
    current_row_T <- list(m = trainSize, L = L_val)
    
    cat(sprintf("正在计算 L=%d (m=%d) ...\n", L_val, trainSize))
    
    for (mod in model_list) {
        
        # 并行计算 (返回 Q 和 T 的报警状态)
        loop_results <- foreach(s = 1:nsim, .combine = rbind, 
                                .export = c("mod", "L_val", "trainSize", "delta_fixed", 
                                            "gamma", "B", "alpha", "b_fixed")) %dopar% {
                                                
                                                yt <- gendataPaper(monitoringPeriod = L_val, trainSize = trainSize, 
                                                                   model = mod, cptLocation = NULL, delta = 0)
                                                
                                                train_data   <- as.matrix(yt[1:trainSize, ])
                                                monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
                                                
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
        
        # 分别计算拒绝率
        rate_qmk <- mean(loop_results[, 1]) * 100
        rate_tmq <- mean(loop_results[, 2]) * 100
        
        # 分别存入对应的列表
        current_row_Q[[paste0("N", mod)]] <- round(rate_qmk, 2)
        current_row_T[[paste0("N", mod)]] <- round(rate_tmq, 2)
        
        cat(sprintf("  -> Model %d | Q: %.1f%% | T: %.1f%%\n", mod, rate_qmk, rate_tmq))
    }
    
    # 将当前行的结果追加到对应的总表
    results_Q <- rbind(results_Q, as.data.frame(current_row_Q))
    results_T <- rbind(results_T, as.data.frame(current_row_T))
}

stopCluster(cl)

# --- 4. 格式化输出 ---

# 添加统计量标识列
results_Q <- cbind(Stat = "Q", results_Q)
results_T <- cbind(Stat = "T", results_T)

# 上下合并
final_df <- rbind(results_Q, results_T)

# 调整列顺序: Stat, m, L, N1...N8
col_order <- c("Stat", "m", "L", paste0("N", model_list))
final_df <- final_df[, col_order]

total_time <- difftime(Sys.time(), startTime, units="mins")
cat(sprintf("\n总耗时: %.2f min\n", as.numeric(total_time)))

cat("\n================ 最终结果表 ================\n")
print(final_df)

# 保存 CSV
write.csv(final_df, "TypeI_Results_QmkTmq.csv", row.names = FALSE)
cat("\n表格已保存至: TypeI_Results_QmkTmq.csv\n")