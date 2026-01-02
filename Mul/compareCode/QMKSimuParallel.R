rm(list = ls())
setwd("D:/Desktop/paperRev/JointCUSUM/Mul/")

# --- 1. 加载必要的包 ---
library(latex2exp)
library(ggplot2)
library(corpcor)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(npcp)
library(foreach)
library(doParallel)

# 加载数据生成函数
# 注意：在并行计算中，确保此文件路径正确，或将函数定义直接写在脚本里
if(file.exists("Gendata4PaperNew.R")){
    source("Gendata4PaperNew.R")
} else {
    stop("找不到 Gendata4PaperNew.R，请检查路径！")
}

# --- 2. 参数设置 ---
DGP <- 306
trainSize <- 100        
T <- trainSize          
monitoringPeriod <- 1   
tau <- 0.2              
p <- 1                  
gamma <- 0              

# 模拟参数
changeAmp <- seq(1, 1, by = 1) # 建议范围稍大一点看全貌
nsim <- 1000   
B <- 500      

# 计算相对变点位置 (相对于监测期开始的索引)
# 监测期长度 = T * L = 100 * 1 = 100
# 变点发生时刻 k* = 100 * 0.2 = 20
monitor_len <- trainSize * monitoringPeriod
k_star <- floor(monitor_len * tau) 

# --- 3. 启动并行核心 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("已启动并行计算，使用核心数: %d\n", num_cores))

# --- 4. 主循环 ---
results_list <- list()

cat("Starting Simulation...\n")
cat(sprintf("变点位置 k* = %d\n", k_star))

for (KK in 1:NROW(changeAmp)) {
    deltachange <- changeAmp[KK]
    startTime <- Sys.time()
    
    # === 并行化内部循环 ===
    loop_result <- foreach(s = 1:nsim, .combine = rbind, 
                           .packages = c("npcp", "mvtnorm", "MTS", "Rcpp"),
                           .export = c("gendataPaper", "k_star", "trainSize", "monitoringPeriod", "DGP", "tau")) %dopar% {
                               
   # 重新加载依赖 (部分系统环境需要)
   # source("Gendata4PaperNew.R") 
   
   # 1. 生成数据
   yt <- gendataPaper(monitoringPeriod = monitoringPeriod, trainSize = trainSize, 
                      model = DGP, cptLocation = tau, delta = deltachange)
   
   m <- trainSize
   n_total <- NROW(yt) 
   
   train <- as.matrix(yt[1:m, ])      
   newXt <- as.matrix(yt[-c(1:m), ])  
   
   # 2. npcp 核心步骤
   # b = NULL 会自动寻优带宽，速度较慢。如果想极速，可固定 b=5 (视数据依赖性而定)
   traj.mult <- npcp::simClosedEndCpDist(x.learn = train, n = n_total, gamma = gamma, 
                                         B = B, method = "mult", b = NULL)
   
   tf.mult <- npcp::threshClosedEndCpDist(traj.mult, p = p, alpha = 0.05)
   det <- npcp::detClosedEndCpDist(x.learn = train, gamma = gamma, x = newXt, delta = 1e-4)
   
   res_tmq <- npcp::monClosedEndCpDist(det, tf.mult, statistic = "mac", plot = FALSE)
   res_qmk <- npcp::monClosedEndCpDist(det, tf.mult, statistic = "mc", plot = FALSE)
   
   # === [FIX] 延迟计算逻辑 ===
   # 逻辑：只有当 (报警时间 > k_star) 时，才算有效检测延迟
   # 否则视为误报或未报警，EDD 记为 NA
   
   # Qmk 延迟
   delay_qmk <- NA
   if (res_qmk$alarm && res_qmk$time.alarm > k_star) {
       delay_qmk <- res_qmk$time.alarm - k_star
   }
   
   # Tmq 延迟
   delay_tmq <- NA
   if (res_tmq$alarm && res_tmq$time.alarm > k_star) {
       delay_tmq <- res_tmq$time.alarm - k_star
   }
   
   # 返回向量: c(Qmk_Alarm, Tmq_Alarm, Qmk_Delay, Tmq_Delay)
   c(as.numeric(res_qmk$alarm), 
     as.numeric(res_tmq$alarm), 
     delay_qmk, 
     delay_tmq)
                           }
    
    # === 处理并行结果 ===
    Tn_Mc  <- loop_result[, 1] # Qmk Alarm Status
    Tn_Mac <- loop_result[, 2] # Tmq Alarm Status
    ARL_Mc <- loop_result[, 3] # Qmk Delay (adjusted)
    ARL_Mac <- loop_result[, 4] # Tmq Delay (adjusted)
    
    # 计算指标
    power_qmk <- round(100 * mean(Tn_Mc), 1)
    power_tmq <- round(100 * mean(Tn_Mac), 1)
    
    # EDD 仅计算有效延迟的平均值 (na.omit 自动剔除误报和未报警)
    edd_qmk   <- mean(na.omit(ARL_Mc))
    edd_tmq   <- mean(na.omit(ARL_Mac))
    
    # 填补 NaN (当 Power=0 时可能出现)
    if(is.nan(edd_qmk)) edd_qmk <- 0
    if(is.nan(edd_tmq)) edd_tmq <- 0
    
    # 存储结果
    results_list[[KK]] <- c(deltachange, power_qmk, power_tmq, edd_qmk, edd_tmq)
    
    # 打印进度
    elapsed <- difftime(Sys.time(), startTime, units="secs")
    cat(sprintf("[Delta=%.1f] 完成 | Qmk Power: %.1f%% | EDD: %.2f | 耗时: %.1fs\n", 
                deltachange, power_qmk, edd_qmk, as.numeric(elapsed)))
}

# --- 5. 关闭并行核心 ---
stopCluster(cl)

# --- 6. 整理最终数据 ---
final_matrix <- do.call(rbind, results_list)
colnames(final_matrix) <- c("Delta", "Power_Qmk", "Power_Tmq", "EDD_Qmk", "EDD_Tmq")
res_df <- as.data.frame(final_matrix)

# 打印表格
print(res_df)

# --- 7. 绘图 (Plotting) ---

# 设置保存 PDF (可选)
# pdf("Simulation_Results.pdf", width=10, height=5)

par(mfrow = c(1, 2)) # 一行两图

# === Plot 1: Power Curve ===
plot(res_df$Delta, res_df$Power_Qmk, type="o", lty=2, pch=1, 
     ylim=c(0,100), col="black",
     xlab=TeX("Change Magnitude ($\\delta$)"), 
     ylab="Rejection Rate (%)", main="Power Curve")
lines(res_df$Delta, res_df$Power_Tmq, type="o", lty=1, pch=2, col="blue")
abline(h=5, col="red", lty=3) # 5% Type I error line
legend("bottomright", legend=c(TeX("$Q_{m,k}$"), TeX("$T_{m,q}$")), 
       col=c("black", "blue"), lty=c(2, 1), pch=c(1, 2))

# === Plot 2: EDD Curve ===
# 动态计算 Y 轴上限，留出空间
max_y <- max(res_df$EDD_Qmk, res_df$EDD_Tmq, na.rm=TRUE)
if(max_y == 0) max_y <- 10

plot(res_df$Delta, res_df$EDD_Qmk, type="o", lty=2, pch=1, 
     ylim=c(0, max_y * 1.1), col="black",
     xlab=TeX("Change Magnitude ($\\delta$)"), 
     ylab="Mean Detection Delay", main="EDD Curve")
lines(res_df$Delta, res_df$EDD_Tmq, type="o", lty=1, pch=2, col="blue")
legend("topright", legend=c(TeX("$Q_{m,k}$"), TeX("$T_{m,q}$")), 
       col=c("black", "blue"), lty=c(2, 1), pch=c(1, 2))

par(mfrow = c(1, 1)) # 恢复默认布局

# dev.off() # 如果使用了 pdf 保存，取消注释此行

# 保存数据到 CSV
# write.csv(res_df, file = "P103Compare_Result.csv", row.names = FALSE)
cat("结果已保存至 Simulation_Results_Parallel.csv\n")