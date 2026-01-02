#########################################################################
#  Model Type I (Optimized Version)
#########################################################################
rm(list = ls())

# --- 加载包 ---
library(npcp)
library(mvtnorm)
library(MTS)
library(foreach)
library(doParallel)
# library(doRNG) # 建议安装并使用此包以保证并行随机数的可重复性

# 确保文件存在
if(!file.exists("Gendata4PaperNew.R")){
    stop("错误: 找不到 Gendata4PaperNew.R")
}

# --- 1. 参数设置 ---
DGP              <- 8       
trainSize        <- 100     
monitoringPeriod <- 1       
delta_fixed      <- 1.0     
gamma            <- 0       
alpha            <- 0.05    
nsim             <- 1000    
B                <- 1000     

# [优化点 1]: 预计算固定带宽 b
# 避免在每次循环内部进行耗时的带宽自动寻优
# 经验法则: b ~ m^(1/3)
b_fixed <- ceiling(trainSize^(1/3)) 
# 或者手动指定：对于 m=100, b=5 是一个合理的选择
# b_fixed <- 5 

cat(sprintf("开始并行测试 | 模型: %d | 样本: %d | Bootstrap块大小: %d\n", 
            DGP, trainSize, b_fixed))

# --- 2. 启动并行 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# registerDoRNG(123) # 如果使用 doRNG，在这里设置种子

cat(sprintf("并行核心已启动: %d\n", num_cores))

# [优化点 2]: 使用 clusterEvalQ 初始化环境
# 这比 .export 更稳定，避免传递 NULL 指针错误
clusterEvalQ(cl, {
    # 设置工作目录 (确保子节点能找到 R 脚本)
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    
    library(npcp)
    library(mvtnorm)
    library(MTS)
    library(Rcpp)
    
    # 加载自定义函数
    source("Gendata4PaperNew.R")
})

# --- 3. 模拟主循环 ---
startTime <- Sys.time()

# 如果安装了 doRNG，将 %dopar% 换成 %dorng%
loop_results <- foreach(s = 1:nsim, .combine = rbind, 
                        # 移除 .export，因为已经通过 clusterEvalQ 加载了
                        .export = c("DGP", "trainSize", "monitoringPeriod", 
                                    "delta_fixed", "gamma", "B", "alpha", "b_fixed")) %dopar% {
                                        
                                        # (A) 生成数据
                                        yt <- gendataPaper(monitoringPeriod = monitoringPeriod, 
                                                           trainSize = trainSize, 
                                                           model = DGP, 
                                                           cptLocation = NULL, 
                                                           delta = delta_fixed)
                                        
                                        # (B) 划分
                                        train_data   <- as.matrix(yt[1:trainSize, ])
                                        monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
                                        
                                        # (C) 训练阶段
                                        # [优化点 1 应用]: 传入固定的 b 参数
                                        traj <- npcp::simClosedEndCpDist(x.learn = train_data, n = nrow(yt), 
                                                                         gamma = gamma, B = B, method = "mult", 
                                                                         b = b_fixed) # <--- 关键修改
                                        
                                        thresh <- npcp::threshClosedEndCpDist(traj, p = 1, alpha = alpha)
                                        
                                        # (D) 计算统计量
                                        det_stat <- npcp::detClosedEndCpDist(x.learn = train_data, x = monitor_data, 
                                                                             gamma = gamma, delta = 1e-4)
                                        
                                        # (E) 监控
                                        res_tmq <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mac", plot = FALSE)
                                        res_qmk <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mc", plot = FALSE)
                                        
                                        c(as.numeric(res_qmk$alarm), as.numeric(res_tmq$alarm))
                                    }

# --- 4. 结束 ---
stopCluster(cl)

# --- 5. 结果 ---
alarm_count_Qmk <- sum(loop_results[, 1])
alarm_count_Tmq <- sum(loop_results[, 2])

final_time <- difftime(Sys.time(), startTime, units="mins")

cat("\n================ 最终结果 (Optimized) ================\n")
cat(sprintf("Qmk Type I Error: %.2f%%\n", (alarm_count_Qmk / nsim) * 100))
cat(sprintf("Tmq Type I Error: %.2f%%\n", (alarm_count_Tmq / nsim) * 100))
cat(sprintf("总耗时: %.2f min\n", as.numeric(final_time)))
cat("=====================================================\n")