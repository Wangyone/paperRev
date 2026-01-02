#########################################################################
#  Model 104 (BEKK-GARCH) Power 测试 - 随机变点版
#########################################################################
rm(list = ls())

library(npcp)
library(mvtnorm)
library(MTS)

# 加载数据生成函数 (确保路径正确)
source("Gendata4PaperNew.R") 

# --- 1. 参数设置 ---
DGP          <- 8    
trainSize    <- 100     # 训练集长度 m
monitoringPeriod <- 1   # 监控期倍数 L (总监控长度 = m * L)
delta_fixed  <- 1.0     # 默认均值偏移幅度 (M4模型中 mu2 增加 delta)
gamma        <- 0       # 权重参数
alpha        <- 0.05    # 显著性水平
nsim         <- 500    # 重复模拟次数
B            <- 100     # Bootstrap 抽样次数 (用于确定阈值)

# 结果容器
alarm_count_Qmk <- 0
alarm_count_Tmq <- 0

cat(sprintf("开始 Power 测试 | 模型: %d | 样本量: %d | 模拟次数: %d\n", 
            DGP, trainSize, nsim))
cat("变点位置：监控期内随机生成 (0 < cptLocation < 0.8)\n")

# --- 2. 模拟循环 ---
set.seed(42) # 保证实验可重复性

for (s in 1:nsim) {
    
    # (A) 生成数据 
    # 注意：cptLocation = NULL 会触发 gendataPaper 内部的随机变点逻辑
    yt <- gendataPaper(monitoringPeriod = monitoringPeriod, 
                       trainSize = trainSize, 
                       model = DGP, 
                       cptLocation = NULL, 
                       delta = delta_fixed)
    
    # (B) 划分数据集
    train_data   <- as.matrix(yt[1:trainSize, ])
    monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
    
    # (C) 训练阶段：计算非参数阈值函数
    # 使用乘性自举 (Multiplier Bootstrap) 处理 BEKK 的依赖结构
    traj <- simClosedEndCpDist(x.learn = train_data, n = nrow(yt), 
                               gamma = gamma, B = B, method = "mult")
    thresh <- threshClosedEndCpDist(traj, p = 1, alpha = alpha)
    
    # (D) 计算检测统计量
    det_stat <- detClosedEndCpDist(x.learn = train_data, x = monitor_data, 
                                   gamma = gamma, delta = 1e-4)
    
    # (E) 实施监控并判断是否拒绝原假设 (Power)
    res_Tmq <- monClosedEndCpDist(det_stat, thresh, statistic = "mac")
    res_Qmk <- monClosedEndCpDist(det_stat, thresh, statistic = "mc")
    
    if (res_Tmq$alarm) alarm_count_Tmq <- alarm_count_Tmq + 1
    if (res_Qmk$alarm) alarm_count_Qmk <- alarm_count_Qmk + 1
    
    # 打印进度
    if (s %% 50 == 0) {
        cat(sprintf("\r运行中: %d/%d | 累计 Power Qmk: %.1f%% | Tmq: %.1f%%", 
                    s, nsim, (alarm_count_Qmk/s)*100, (alarm_count_Tmq/s)*100))
    }
}

# --- 3. 结果汇总 ---
cat("\n\n================ 最终 Power 测试结果 ================\n")
cat(sprintf("测试模型: %d (BEKK-GARCH)\n", DGP))
cat(sprintf("均值偏移: %.1f | 显著性水平: %.2f\n", delta_fixed, alpha))
cat("----------------------------------------------------\n")
cat(sprintf("Qmk (Statistic 'mc')  Power: %.2f%%\n", (alarm_count_Qmk / nsim) * 100))
cat(sprintf("Tmq (Statistic 'mac') Power: %.2f%%\n", (alarm_count_Tmq / nsim) * 100))
cat("====================================================\n")