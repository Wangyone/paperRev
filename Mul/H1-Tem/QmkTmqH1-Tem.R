#########################################################################
# NPCP Structural Change Simulation (Model 301-306)
# Statistics: Qmk & Tmq
# Output Format: Two Wide Tables (One for Qmk, One for Tmq)
#########################################################################

rm(list = ls())
# 请根据实际路径修改
setwd("D:/Desktop/paperRev/JointCUSUM/Mul/")

# --- 加载依赖包 ---
library(latex2exp)
library(mvtnorm)
library(npcp)
library(foreach)
library(doParallel)

# --- 1. 检查依赖文件 ---
if(!file.exists("Gendata4PaperNew.R")) stop("Gendata4PaperNew.R missing")
source("Gendata4PaperNew.R")

# --- 2. 参数配置 ---
m_list      <- c(100, 200)       # 训练样本大小
L_list      <- c(1, 2, 3)        # 监控周期倍数
model_list  <- 301:306           # 结构变化模型
cpt_loc     <- 0.2               # 变点位置

nsim        <- 100               # 模拟次数
B           <- 10               # Bootstrap 次数
p           <- 1                 # 阈值函数参数
gamma       <- 0                 # 权重参数

# --- 3. 初始化结果容器 ---
# 我们需要两个表: result_Qmk 和 result_Tmq

# 生成列名: m, L, N301-Power, N301-MDD, N302-Power ...
base_cols <- c("m", "L")
model_cols <- c()
for(mod in model_list) {
    model_cols <- c(model_cols, paste0("N", mod, "-Power"), paste0("N", mod, "-MDD"))
}
col_names <- c(base_cols, model_cols)

total_rows <- length(m_list) * length(L_list)

# 创建两个空数据框
result_Qmk <- as.data.frame(matrix(NA, nrow = total_rows, ncol = length(col_names)))
result_Tmq <- as.data.frame(matrix(NA, nrow = total_rows, ncol = length(col_names)))

colnames(result_Qmk) <- col_names
colnames(result_Tmq) <- col_names

# --- 4. 启动并行集群 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
if(num_cores < 1) num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 将环境推送到子节点
clusterEvalQ(cl, {
    setwd("D:/Desktop/paperRev/JointCUSUM/Mul/")
    library(mvtnorm)
    library(npcp)
    source("Gendata4PaperNew.R")
})

# 导出变量
clusterExport(cl, c("cpt_loc", "B", "p", "gamma"))

# --- 5. 主循环 ---
cat("==========================================================\n")
cat("Start NPCP Simulation (Models 301-306, Qmk & Tmq)\n")
cat(sprintf("Nsim: %d | Bootstrap: %d\n", nsim, B))
cat("==========================================================\n\n")

row_idx <- 0 

for (m_val in m_list) {
    for (L_val in L_list) {
        
        row_idx <- row_idx + 1 
        
        # 填充基础列
        result_Qmk[row_idx, "m"] <- m_val
        result_Qmk[row_idx, "L"] <- L_val
        result_Tmq[row_idx, "m"] <- m_val
        result_Tmq[row_idx, "L"] <- L_val
        
        # 计算当前场景下的关键时间点
        monitor_len  <- m_val * L_val
        k_star       <- floor(monitor_len * cpt_loc) # 相对监控起点的变点时刻
        
        cat(sprintf("[Scenario %d/%d] m=%d, L=%d (Change at k*=%d) ... ", 
                    row_idx, total_rows, m_val, L_val, k_star))
        
        for (mod_id in model_list) {
            
            # === 并行模拟核心 ===
            sim_res <- foreach(s = 1:nsim, .packages = c('mvtnorm', 'npcp'), 
                               .export = c('m_val', 'L_val', 'mod_id', 'monitor_len', 'k_star')) %dopar% {
                                   
                                   # 1. 生成数据
                                   yt <- gendataPaper(monitoringPeriod = L_val, trainSize = m_val, 
                                                      model = mod_id, cptLocation = cpt_loc)
                                   
                                   # 2. 区分训练集与监控集
                                   train <- as.matrix(yt[1:m_val, ])
                                   newXt <- as.matrix(yt[-c(1:m_val), ])
                                   n_total <- nrow(yt)
                                   
                                   # 3. NPCP 步骤
                                   traj.mult <- simClosedEndCpDist(x.learn = train, n = n_total, gamma = gamma, 
                                                                   B = B, method = "mult", b = NULL)
                                   
                                   tf.mult <- threshClosedEndCpDist(traj.mult, p = p, alpha = 0.05)
                                   
                                   det <- detClosedEndCpDist(x.learn = train, gamma = gamma, x = newXt, delta = 1e-4)
                                   
                                   # 实施监控
                                   res_Q <- monClosedEndCpDist(det, tf.mult, statistic = "mc", plot = FALSE)
                                   res_T <- monClosedEndCpDist(det, tf.mult, statistic = "mac", plot = FALSE)
                                   
                                   # 4. 提取结果 (Qmk)
                                   alarm_Q <- 0; delay_Q <- NA
                                   if (res_Q$alarm) {
                                       alarm_Q <- 1
                                       if (res_Q$time.alarm > k_star) {
                                           delay_Q <- res_Q$time.alarm - k_star
                                       }
                                   }
                                   
                                   # 5. 提取结果 (Tmq)
                                   alarm_T <- 0; delay_T <- NA
                                   if (res_T$alarm) {
                                       alarm_T <- 1
                                       if (res_T$time.alarm > k_star) {
                                           delay_T <- res_T$time.alarm - k_star
                                       }
                                   }
                                   
                                   c(alarm_Q, delay_Q, alarm_T, delay_T)
                               }
            
            # === 汇总结果 ===
            res_matrix <- do.call(rbind, sim_res)
            
            # Qmk 指标
            pow_Q <- mean(res_matrix[, 1], na.rm=TRUE) * 100
            delays_Q <- res_matrix[, 2]
            mdd_Q <- mean(delays_Q[!is.na(delays_Q)], na.rm=TRUE)
            
            # Tmq 指标
            pow_T <- mean(res_matrix[, 3], na.rm=TRUE) * 100
            delays_T <- res_matrix[, 4]
            mdd_T <- mean(delays_T[!is.na(delays_T)], na.rm=TRUE)
            
            # 填表
            col_pow <- paste0("N", mod_id, "-Power")
            col_mdd <- paste0("N", mod_id, "-MDD")
            
            result_Qmk[row_idx, col_pow] <- round(pow_Q, 2)
            result_Qmk[row_idx, col_mdd] <- round(mdd_Q, 2)
            
            result_Tmq[row_idx, col_pow] <- round(pow_T, 2)
            result_Tmq[row_idx, col_mdd] <- round(mdd_T, 2)
            
            cat(sprintf("N%d ", mod_id))
            
        } # end models
        cat("\n")
        
    } # end L
} # end m

stopCluster(cl)

# --- 6. 结果保存与展示 ---
cat("\n========================================\n")
cat("Simulation Finished.\n")
cat("========================================\n")

cat("\n>>> Qmk Results (Preview):\n")
print(head(result_Qmk))

cat("\n>>> Tmq Results (Preview):\n")
print(head(result_Tmq))

# 保存为两个分开的文件
write.csv(result_Qmk, file = "Temporal_Qmk_Power_MDD.csv", row.names = FALSE)
write.csv(result_Tmq, file = "Temporal_Tmq_Power_MDD.csv", row.names = FALSE)

cat("\nResults saved to 'Temporal_Qmk_Power_MDD.csv' and 'Temporal_Tmq_Power_MDD.csv'.\n")
