rm(list = ls())
######################### Type1 error / Power Simulation #########################
# 设置工作路径 (请根据实际情况修改)
setwd("D:/Desktop/paperRev/JointCUSUM/Mul/")

# 加载依赖包
library(latex2exp)
library(ggplot2)
library(corpcor)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(npcp)

# 加载数据生成函数
# 请确保此文件存在
if(file.exists("Gendata4PaperNew.R")){
    source("Gendata4PaperNew.R")
} else {
    stop("找不到 Gendata4PaperNew.R，请检查路径！")
}

## --- 参数设置 ---

DGP <- 104

# 基础参数
trainSize <- 100       # 训练样本大小 (m)
T <- trainSize         # 定义 T 防止报错

monitoringPeriod <- 1  # 监控期倍数因子 (L)
tau <- 0.2             # 变点位置 (相对于监控期开始的比例)

p <- 1                 # 阈值函数的步数
gamma <- 0             # 权重参数

# --- 计算相对变点位置 (关键修复步骤) ---
# 监控期总长度 = 100 * 1 = 100
# 变点发生时刻 k* = 100 * 0.2 = 20 (即变点发生在监测期的第20个点)
monitor_len <- trainSize * monitoringPeriod
k_star <- floor(monitor_len * tau) 

cat(sprintf("变点发生在监测期的第 %d 个时刻 (k*)\n", k_star))

# --- 模拟参数 ---
changeAmp <- seq(0, 4, by = 0.5) # 遍历的变点幅度 delta

# 结果存储容器
RQMK <- numeric(NROW(changeAmp))
RTMQ <- numeric(NROW(changeAmp))
EDQMK <- numeric(NROW(changeAmp))
EDTMQ <- numeric(NROW(changeAmp))

nsim <- 500  # 模拟次数
B <- 100     # Bootstrap 次数

cat("Starting Simulation...\n")

########################### 主循环：遍历 Delta ###########################
for (KK in 1:NROW(changeAmp)) {
    deltachange <- changeAmp[KK]
    
    # 初始化统计量存储矩阵
    Tn_Mac <- matrix(0, nrow=nsim, ncol=1) # 记录是否拒绝 (0/1)
    Tn_Mc  <- matrix(0, nrow=nsim, ncol=1)
    
    # 初始化为 NA，这样未报警或误报的实验在计算 EDD 时会被 na.omit 自动剔除
    ARL_Mac <- matrix(NA, nrow=nsim, ncol=1) # 记录延迟 (Delay)
    ARL_Mc  <- matrix(NA, nrow=nsim, ncol=1)
    
    startTime <- Sys.time()
    
    # --- 内部循环：蒙特卡洛模拟 ---
    for (s in 1:nsim) {
        # 1. 生成数据
        yt <- gendataPaper(monitoringPeriod = monitoringPeriod, trainSize = trainSize, 
                           model = DGP, cptLocation = tau, delta = deltachange)
        
        # 2. 区分训练集与监控集
        m <- trainSize                
        n <- NROW(yt)                 # 总长度 (m + 监控长度)
        
        train <- as.matrix(yt[1:m, ])      # 历史数据
        newXt <- as.matrix(yt[-c(1:m), ])  # 新数据
        
        # 3. npcp 步骤
        # Step 1: 模拟 Null 分布
        traj.mult <- simClosedEndCpDist(x.learn = train, n = n, gamma = gamma, 
                                        B = B, method = "mult", b = NULL)
        
        # Step 2: 计算阈值函数
        tf.mult <- threshClosedEndCpDist(traj.mult, p = p, alpha = 0.05)
        
        # Step 3: 计算检测器统计量
        det <- detClosedEndCpDist(x.learn = train, gamma = gamma, x = newXt, delta = 1e-4)
        
        # Step 4: 实施监控
        result1 <- monClosedEndCpDist(det, tf.mult, statistic = "mac", plot = FALSE) # Tmq
        result2 <- monClosedEndCpDist(det, tf.mult, statistic = "mc", plot = FALSE)  # Qmk
        
        # 4. 记录结果
        Tn_Mac[s] <- as.numeric(result1$alarm)
        Tn_Mc[s]  <- as.numeric(result2$alarm)
        
        # === [核心修复] EDD 计算逻辑 ===
        # 逻辑：只有当 (报警时间 > k_star) 时，才算有效延迟
        # 计算公式：Delay = 报警时刻 - k_star
        
        # Tmq 统计量
        if (result1$alarm) {
            if (result1$time.alarm > k_star) {
                ARL_Mac[s] <- result1$time.alarm - k_star
            } 
            # 如果 <= k_star，视为误报(False Alarm)，保持为 NA，不计入 EDD
        }
        
        # Qmk 统计量
        if (result2$alarm) {
            if (result2$time.alarm > k_star) {
                ARL_Mc[s]  <- result2$time.alarm - k_star
            }
        }
        
        # --- 进度打印 ---
        if (s %% 50 == 0) {
            curr_rej0 <- sum(Tn_Mc[1:s])
            curr_rej1 <- sum(Tn_Mac[1:s])
            cat(sprintf("\r[Delta=%.1f] Sim %d/%d | Rej(Qmk): %.1f%% | Rej(Tmq): %.1f%%", 
                        deltachange, s, nsim, (curr_rej0/s)*100, (curr_rej1/s)*100))
        }
    } # end nsim loop
    
    cat("\n") 
    
    # --- 计算最终指标 ---
    
    # 1. 拒绝率 (Power)
    RQMK[KK] <- round(100 * mean(Tn_Mc), 1)
    RTMQ[KK] <- round(100 * mean(Tn_Mac), 1)
    
    # 2. 平均检测延迟 (EDD)
    # 使用 na.omit 自动排除：1.未报警的情况 2.误报的情况(<=k_star)
    EDQMK[KK] <- mean(na.omit(ARL_Mc))
    EDTMQ[KK] <- mean(na.omit(ARL_Mac))
    
    # 防止 NaN (当 Power=0 或全为误报时)
    if(is.nan(EDQMK[KK])) EDQMK[KK] <- 0
    if(is.nan(EDTMQ[KK])) EDTMQ[KK] <- 0
    
} # end Delta loop

# --- 结果打印 ---
print(data.frame(Delta=changeAmp, Power_Qmk=RQMK, EDD_Qmk=EDQMK, EDD_Tmq=EDTMQ))

# --- 绘图部分 ---

# 1. Power Curve
plot(changeAmp, RQMK, type="l", lty=2, pch=1, ylim=c(0,100), xlab=TeX("$\\delta$"), 
     ylab="Rejection percentage", main="Power Curve")
lines(changeAmp, RTMQ, lty=1, pch=2, col="blue")
legend("bottomright", legend=c(TeX("$Q_{m,k}$"), TeX("$T_{m,q}$")), 
       col=c("black", "blue"), lty=c(2, 1), pch=c(1, 2))

# 2. EDD Curve
# 动态调整 Y 轴范围
ylim_range <- range(c(EDQMK, EDTMQ), na.rm=TRUE)
if(max(ylim_range) == 0) ylim_range <- c(0, 10)

plot(changeAmp, EDQMK, type="l", lty=2, pch=1, xlab=TeX("$\\delta$"), 
     ylab="Mean detection delay", ylim=ylim_range, main="EDD Curve")
lines(changeAmp, EDTMQ, lty=1, pch=2, col="blue")
legend("topright", legend=c(TeX("$Q_{m,k}$"), TeX("$T_{m,q}$")), 
       col=c("black", "blue"), lty=c(2, 1), pch=c(1, 2))

# --- 保存结果 ---
P103 <- rbind(changeAmp, RQMK, RTMQ, EDQMK, EDTMQ)
# write.csv(P103, file = "P103QMK_Result_Fixed.csv", row.names = TRUE)

cat("Simulation Finished Successfully.\n")