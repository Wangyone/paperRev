#########################################################################
# Multi-Model Power & MDD (Q & Tmq - NPCP Version)
# Fixed MDD Calculation: Delay = Absolute Alarm Time - Absolute Change Point
#########################################################################
rm(list = ls())
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(npcp)
library(mvtnorm)
library(MTS)
library(foreach)
library(doParallel)
library(ggplot2)

# --- 1. 检查依赖 ---
if(!file.exists("Gendata4PaperNew.R")){
    stop("Error: Gendata4PaperNew.R not found!")
}
source("Gendata4PaperNew.R") 

# --- 2. 参数设置 ---
trainSize     <- 100      
L_val         <- 1        
model_list    <- 101:104  # 目标模型列表 (M1-M4)
deltas        <- seq(0, 2, by = 0.1) 

# NPCP 参数
gamma         <- 0        
alpha         <- 0.05    
nsim          <- 1000      # 模拟次数 (建议 400-1000)
B             <- 500      # Bootstrap 次数
b_fixed       <- ceiling(trainSize^(1/3)) 

# 【关键】：变点时刻设置
# 监测期从 m+1 开始。长度 = m*L = 100。
# 变点在监测期的 20% 处：m + 0.2*100 = 120。
cp_ratio      <- 0.2
monitor_len   <- floor(L_val * trainSize)
true_cp_time  <- trainSize + floor(cp_ratio * monitor_len) # 120

cat("==========================================================\n")
cat("Start Power/MDD Simulation (Q & Tmq)\n")
cat(sprintf("Config: m=%d | L=%d | Change at absolute time: %d\n", trainSize, L_val, true_cp_time))
cat("==========================================================\n\n")

# --- 3. 启动并行 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 环境预热
clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    library(npcp); library(mvtnorm); library(MTS)
    source("Gendata4PaperNew.R")
})

# --- 4. 主循环 ---
final_results <- data.frame()

for (mod in model_list) {
    
    cat(sprintf("\n>>> Processing Model %d ...\n", mod))
    
    for (d_val in deltas) {
        
        # 并行计算 (返回: Q_Rej, Q_Delay, T_Rej, T_Delay)
        loop_results <- foreach(s = 1:nsim, .combine = rbind, 
                                .export = c("mod", "L_val", "trainSize", "d_val", "cp_ratio", "true_cp_time",
                                            "gamma", "B", "alpha", "b_fixed", "monitor_len")) %dopar% {
                                                
                                                # 1. 生成数据 (H1)
                                                yt <- gendataPaper(monitoringPeriod = L_val, trainSize = trainSize, 
                                                                   model = mod, cptLocation = cp_ratio, delta = d_val)
                                                
                                                train_data   <- as.matrix(yt[1:trainSize, ])
                                                monitor_data <- as.matrix(yt[(trainSize + 1):nrow(yt), ])
                                                
                                                # 2. NPCP 核心流程
                                                traj <- npcp::simClosedEndCpDist(x.learn = train_data, n = nrow(yt), 
                                                                                 gamma = gamma, B = B, method = "mult", 
                                                                                 b = b_fixed)
                                                
                                                tf_mult <- npcp::threshClosedEndCpDist(traj, p = 1, alpha = alpha)
                                                
                                                det_stat <- npcp::detClosedEndCpDist(x.learn = train_data, x = monitor_data, 
                                                                                     gamma = gamma, delta = 1e-4)
                                                
                                                res_qmk <- npcp::monClosedEndCpDist(det_stat, tf_mult, statistic = "mc", plot = FALSE)
                                                res_tmq <- npcp::monClosedEndCpDist(det_stat, tf_mult, statistic = "mac", plot = FALSE)
                                                
                                                # 3. 提取结果 & 计算延迟
                                                # 注意: time.alarm 是绝对索引 (1...n)
                                                
                                                # --- Qmk ---
                                                q_rej <- if(res_qmk$alarm) 1 else 0
                                                q_delay <- NA
                                                if(q_rej) {
                                                    # 直接使用绝对时间
                                                    alarm_abs_time <- res_qmk$time.alarm
                                                    # Delay = 报警时间 - 变点时间
                                                    delay <- alarm_abs_time - true_cp_time
                                                    # 排除误报 (False Alarm)
                                                    if(delay > 0) q_delay <- delay
                                                }
                                                
                                                # --- Tmq ---
                                                t_rej <- if(res_tmq$alarm) 1 else 0
                                                t_delay <- NA
                                                if(t_rej) {
                                                    alarm_abs_time <- res_tmq$time.alarm
                                                    delay <- alarm_abs_time - true_cp_time
                                                    if(delay > 0) t_delay <- delay
                                                }
                                                
                                                c(q_rej, q_delay, t_rej, t_delay)
                                            }
        
        # --- 统计指标 ---
        
        # Power (包含误报)
        q_power <- mean(loop_results[, 1], na.rm = TRUE) * 100
        t_power <- mean(loop_results[, 3], na.rm = TRUE) * 100
        
        # MDD (Mean Detection Delay)
        # 排除未拒绝 (NA) 和 误报 (已经在 loop 中处理为 NA)
        q_delays_vec <- loop_results[, 2]
        q_mdd <- mean(q_delays_vec[!is.na(q_delays_vec)], na.rm = TRUE)
        if(is.nan(q_mdd)) q_mdd <- NA
        
        t_delays_vec <- loop_results[, 4]
        t_mdd <- mean(t_delays_vec[!is.na(t_delays_vec)], na.rm = TRUE)
        if(is.nan(t_mdd)) t_mdd <- NA
        
        cat(sprintf("   Delta=%.1f | Q: Pow=%.1f%% MDD=%.2f | Tmq: Pow=%.1f%% MDD=%.2f\n", 
                    d_val, q_power, q_mdd, t_power, t_mdd))
        
        # 保存结果
        final_results <- rbind(final_results, data.frame(
            Model = mod,
            Delta = d_val,
            Method = "Qmk",
            Power = q_power,
            MDD = q_mdd
        ))
        final_results <- rbind(final_results, data.frame(
            Model = mod,
            Delta = d_val,
            Method = "Tmq",
            Power = t_power,
            MDD = t_mdd
        ))
    }
}

stopCluster(cl)

# --- 5. 结果处理 ---
cat("\n========================================\n")
cat("Simulation Finished.\n")
print(head(final_results))


fileSaveName <- paste0("Power_MDD_Results_NPCP(", head(model_list, 1), "-", tail(model_list, 1), ")")
write.csv(final_results, file = paste0(fileSaveName, ".csv"), row.names = FALSE)



# 简单绘图
if(nrow(final_results) > 0) {
    p1 <- ggplot(final_results, aes(x=Delta, y=Power, color=Method, linetype=factor(Model))) +
        geom_line() + geom_point() +
        labs(title="Power Comparison (Q vs Tmq)", y="Power (%)") +
        theme_minimal()
    
    p2 <- ggplot(final_results, aes(x=Delta, y=MDD, color=Method, linetype=factor(Model))) +
        geom_line() + geom_point() +
        labs(title="MDD Comparison (Q vs Tmq)", y="Mean Detection Delay") +
        theme_minimal()
    
    print(p1)
    print(p2)
    # ggsave("NPCP_Power.png", p1, width=7, height=5)
    # ggsave("NPCP_MDD.png", p2, width=7, height=5)
}