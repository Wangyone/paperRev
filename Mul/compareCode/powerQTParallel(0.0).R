#########################################################################
#  Model Type I
#########################################################################
rm(list = ls())

library(npcp)
library(mvtnorm)
library(MTS)
library(foreach)
library(doParallel)

# 加载数据生成函数
if(file.exists("Gendata4PaperNew.R")){
    source("Gendata4PaperNew.R") 
} else {
    stop("找不到 Gendata4PaperNew.R，请检查路径！")
}

# --- 1. 参数设置 ---
DGP              <- 8      
trainSize        <- 100     # 训练集长度 m
monitoringPeriod <- 1       # 监控期倍数 L
delta_fixed      <- 1.0     # 均值偏移幅度
gamma            <- 0       # 权重参数
alpha            <- 0.05    # 显著性水平
nsim             <- 1000     # 重复模拟次数
B                <- 500     # Bootstrap 抽样次数

cat(sprintf("开始并行 Power 测试 | 模型: %d | 核心数预备...\n", DGP))

# --- 2. 启动并行核心 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("已启动并行计算，使用核心数: %d\n", num_cores))

# --- 3. 并行模拟循环 ---
startTime <- Sys.time()

# 使用 foreach 并行执行 nsim 次模拟
loop_results <- foreach(s = 1:nsim, .combine = rbind, 
                        .packages = c("npcp", "mvtnorm", "MTS"),
                        .export = c("gendataPaper", "DGP", "trainSize", 
                                    "monitoringPeriod", "delta_fixed", "gamma", "B", "alpha")) %dopar% {
    
    # (A) 生成数据 (cptLocation = NULL 触发随机变点)
    yt <- gendataPaper(monitoringPeriod = monitoringPeriod, 
                       trainSize = trainSize, 
                       model = DGP, 
                       cptLocation = NULL, 
                       delta = delta_fixed)
    
    # (B) 划分数据集
    train_data   <- as.matrix(yt[1:trainSize, ])
    monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
    
    # (C) 训练阶段：计算非参数阈值
    # 在并行环境中，显式调用 npcp:: 确保稳定
    traj <- npcp::simClosedEndCpDist(x.learn = train_data, n = nrow(yt), 
                                     gamma = gamma, B = B, method = "mult")
    thresh <- npcp::threshClosedEndCpDist(traj, p = 1, alpha = alpha)
    
    # (D) 计算检测统计量
    det_stat <- npcp::detClosedEndCpDist(x.learn = train_data, x = monitor_data, 
                                         gamma = gamma, delta = 1e-4)
    
    # (E) 实施监控
    res_tmq <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mac", plot = FALSE)
    res_qmk <- npcp::monClosedEndCpDist(det_stat, thresh, statistic = "mc", plot = FALSE)
    
    # 返回报警状态 (1表示报警，0表示未报警)
    c(as.numeric(res_qmk$alarm), as.numeric(res_tmq$alarm))
                                    }

# --- 4. 关闭并行核心 ---
stopCluster(cl)

# --- 5. 结果汇总 ---
alarm_count_Qmk <- sum(loop_results[, 1])
alarm_count_Tmq <- sum(loop_results[, 2])

final_time <- difftime(Sys.time(), startTime, units="mins")

cat("\n\n================ 最终 Power 测试结果 ================\n")
cat(sprintf("测试模型: %d (BEKK-GARCH)\n", DGP))
cat(sprintf("模拟次数: %d | 总耗时: %.2f min\n", nsim, as.numeric(final_time)))
cat(sprintf("均值偏移: %.1f | 显著性水平: %.2f\n", delta_fixed, alpha))
cat("----------------------------------------------------\n")
cat(sprintf("Qmk (Statistic 'mc')  Power: %.2f%%\n", (alarm_count_Qmk / nsim) * 100))
cat(sprintf("Tmq (Statistic 'mac') Power: %.2f%%\n", (alarm_count_Tmq / nsim) * 100))
cat("====================================================\n")