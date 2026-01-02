#########################################################################
#  Multi-Model Type I Error Simulation (Final Fixed Version)
#########################################################################
rm(list = ls())

# --- 加载必要的包 ---
library(npcp)
library(mvtnorm)
library(MTS)
library(foreach)
library(doParallel)

# 检查脚本是否存在
if(!file.exists("Gendata4PaperNew.R")){
    stop("错误: 找不到 Gendata4PaperNew.R，请检查路径！")
}

# --- 1. 全局参数设置 ---
trainSize        <- 100     # 固定训练样本 m
L_list           <- c(1, 2, 3) # 遍历监控倍数
model_list       <- 1:8     # 遍历模型 N1 - N8

# 核心参数
delta_fixed      <- 1.0     # 【关键变量】均值偏移幅度
gamma            <- 0       
alpha            <- 0.05    
nsim             <- 1000    # 模拟次数
B                <- 500     # Bootstrap 次数

# [优化]: 预计算固定带宽 b ~ m^(1/3)
b_fixed <- ceiling(trainSize^(1/3))  

cat("==========================================================\n")
cat("开始多模型并行模拟 (Final Fixed)\n")
cat(sprintf("Fixed m: %d | Block Size b: %d | Sim: %d\n", trainSize, b_fixed, nsim))
cat("==========================================================\n\n")

# --- 2. 启动并行集群 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat(sprintf("并行核心已启动: %d\n", num_cores))

# 环境预热：让每个 Worker 加载函数库
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') # 确保路径正确
    library(npcp)
    library(mvtnorm)
    library(MTS)
    library(Rcpp)
    source("Gendata4PaperNew.R")
})

# --- 3. 主循环 ---
startTime <- Sys.time()
final_results_df <- data.frame()

# 外层循环：监控周期 L
for (L_val in L_list) {
    
    current_row <- list(m = trainSize, L = L_val)
    cat(sprintf("正在计算 L=%d (m=%d) ...\n", L_val, trainSize))
    
    # 内层循环：模型 Model
    for (mod in model_list) {
        
        # 【修复重点】：显式使用 .export 导出所有必要的变量
        loop_results <- foreach(s = 1:nsim, .combine = rbind, 
                                .export = c("mod", "L_val", "trainSize", "delta_fixed", 
                                            "gamma", "B", "alpha", "b_fixed")) %dopar% {
                                                
                                                # (A) 生成数据
                                                # 这里用到 delta_fixed，之前报错就是因为子进程找不到它
                                                yt <- gendataPaper(monitoringPeriod = L_val, 
                                                                   trainSize = trainSize, 
                                                                   model = mod, 
                                                                   cptLocation = NULL, 
                                                                   delta = 0) # Null模型下delta其实不生效，但函数签名需要
                                                
                                                # (B) 划分
                                                train_data   <- as.matrix(yt[1:trainSize, ])
                                                monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
                                                
                                                # (C) 训练阶段 (使用固定 b 加速)
                                                traj <- npcp::simClosedEndCpDist(x.learn = train_data, n = nrow(yt), 
                                                                                 gamma = gamma, B = B, method = "mult", 
                                                                                 b = b_fixed)
                                                
                                                thresh <- npcp::threshClosedEndCpDist(traj, p = 1, alpha = alpha)
                                                
                                                # (D) 计算统计量
                                                det_stat <- npcp::detClosedEndCpDist(x.learn = train_data, x = monitor_data, 
                                                                                     gamma = gamma, delta = 1e-4)
                                                
                                                # (E) 监控
                                                res_tmq <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mac", plot = FALSE)
                                                
                                                # 只返回 Tmq 结果 (1=报警, 0=未报警)
                                                as.numeric(res_tmq$alarm)
                                            }
        
        # loop_results 现在是一个向量 (因为 .combine=rbind 对单列结果会转为向量) 
        # 或者矩阵。为了保险，强制视为向量计算均值
        rate_tmq <- mean(as.numeric(loop_results)) * 100
        
        # 记录结果
        current_row[[paste0("N", mod)]] <- round(rate_tmq, 2)
        cat(sprintf("  -> Model %d Done. Rate: %.2f%%\n", mod, rate_tmq))
    }
    
    # 合并结果
    final_results_df <- rbind(final_results_df, as.data.frame(current_row))
}

# --- 4. 结束清理 ---
stopCluster(cl)

# --- 5. 结果输出 ---
total_time <- difftime(Sys.time(), startTime, units="mins")

cat("\n================ 最终结果汇总 ================\n")
cat(sprintf("总耗时: %.2f min\n", as.numeric(total_time)))

# 打印表格
print(final_results_df)

# 保存 CSV
filename <- "TypeI_Results_m100_Fixed.csv"
write.csv(final_results_df, filename, row.names = FALSE)
cat(sprintf("\n表格已保存至: %s\n", filename))