#########################################################################
#  Multi-Model Type I Error (Both Qmk & Tmq)
#########################################################################
rm(list = ls())

library(npcp)
library(mvtnorm)
library(MTS)
library(foreach)
library(doParallel)

if(!file.exists("Gendata4PaperNew.R")){
    stop("错误: 找不到 Gendata4PaperNew.R")
}

# --- 1. 参数设置 ---
trainSize        <- 100     
L_list           <- c(1, 2, 3) 
model_list       <- 1:8     

delta_fixed      <- 1.0     
gamma            <- 0       
alpha            <- 0.05    
nsim             <- 300    
B                <- 50     
b_fixed          <- 5       

cat("==========================================================\n")
cat("开始多模型并行模拟 (Qmk + Tmq)\n")
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
final_results_df <- data.frame()

for (L_val in L_list) {
    
    # 初始化当前行
    current_row <- list(m = trainSize, L = L_val)
    cat(sprintf("正在计算 L=%d (m=%d) ...\n", L_val, trainSize))
    
    for (mod in model_list) {
        
        # 【修改 1】：并行返回两个结果
        loop_results <- foreach(s = 1:nsim, .combine = rbind, 
                                .export = c("mod", "L_val", "trainSize", "delta_fixed", 
                                            "gamma", "B", "alpha", "b_fixed")) %dopar% {
                                                
                                                # (A) 生成数据
                                                yt <- gendataPaper(monitoringPeriod = L_val, trainSize = trainSize, 
                                                                   model = mod, cptLocation = NULL, delta = 0)
                                                
                                                # (B) 划分
                                                train_data   <- as.matrix(yt[1:trainSize, ])
                                                monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
                                                
                                                # (C) 训练 (固定 b)
                                                traj <- npcp::simClosedEndCpDist(x.learn = train_data, n = nrow(yt), 
                                                                                 gamma = gamma, B = B, method = "mult", 
                                                                                 b = b_fixed)
                                                thresh <- npcp::threshClosedEndCpDist(traj, p = 1, alpha = alpha)
                                                
                                                # (D) 统计量
                                                det_stat <- npcp::detClosedEndCpDist(x.learn = train_data, x = monitor_data, 
                                                                                     gamma = gamma, delta = 1e-4)
                                                
                                                # (E) 监控 (同时计算 Qmk 和 Tmq)
                                                res_qmk <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mc", plot = FALSE)
                                                res_tmq <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mac", plot = FALSE)
                                                
                                                # 返回向量: [Qmk结果, Tmq结果]
                                                c(as.numeric(res_qmk$alarm), as.numeric(res_tmq$alarm))
                                            }
        
        # 【修改 2】：分别计算两列的拒绝率
        rate_qmk <- mean(loop_results[, 1]) * 100
        rate_tmq <- mean(loop_results[, 2]) * 100
        
        # 存入表格 (为了区分，列名加上后缀)
        # 格式：N1_Q (Qmk结果), N1_T (Tmq结果)
        current_row[[paste0("N", mod, "_Q")]] <- round(rate_qmk, 2)
        current_row[[paste0("N", mod, "_T")]] <- round(rate_tmq, 2)
        
        cat(sprintf("  -> Model %d | Qmk: %.1f%% | Tmq: %.1f%%\n", mod, rate_qmk, rate_tmq))
    }
    
    final_results_df <- rbind(final_results_df, as.data.frame(current_row))
}

stopCluster(cl)

# --- 4. 结果输出 ---
total_time <- difftime(Sys.time(), startTime, units="mins")
cat(sprintf("\n总耗时: %.2f min\n", as.numeric(total_time)))

print(final_results_df)

# 保存 CSV
write.csv(final_results_df, "TypeI_Results_TmqQmk.csv", row.names = FALSE)
cat("表格已保存至: TypeI_Results_Both_Stats.csv\n")