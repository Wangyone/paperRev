#########################################################################
# Multi-Model Power & MDD (statL - Wrap-Speed Version)
# MDD = Mean Detection Delay (Alarm Time - Change Point)
# Update: Bandwidth delta0 = 1.0 * d (Dimension Adaptive) 
#########################################################################
rm(list = ls())
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(foreach)
library(doParallel)
library(ggplot2)

# --- 1. 检查依赖文件 ---
if(!file.exists("Gendata4PaperNew.R")) stop("Gendata4PaperNew.R missing")
if(!file.exists("statLMul.cpp")) stop("statLMul.cpp missing") 

# --- 2. 参数设置 ---
m_val       <- 100       # 训练样本 (m)
L_val       <- 1         # 监控倍数 (L)
true_cp_idx <- 20        # 变点相对索引 (20th point in monitoring)
true_cp_time<- m_val + true_cp_idx # 绝对变点时刻: 100 + 20 = 120

model_list  <- 101:104   # 均值变点模型
deltas      <- seq(0, 2, by = 0.1)
# model_list  <- 201:203   # 相关系数变点模型
# deltas      <- seq(0, 0.9, by = 0.05) 
# model_list  <- 211:214   # 方差变点模型
# deltas      <- seq(0, 4, by = 0.1) 

nsim        <- 2000      # 模拟次数
TypeI       <- 0.05
kernal      <- "quadexp" # gauss quadexp energy
alphaE      <- 1
m0          <- 0          

# 辅助函数
Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
Tm <- m_val - m0
N_total <- (1 + L_val) * m_val
k_vec <- 1:(N_total - m_val)
# 预计算权重
gamma_boundary <- 0.0
vaha_vec <- (Tm + k_vec)^2 / (Tm * Qk(k_vec/Tm, gamma_boundary)^2)

cat("==========================================================\n")
cat("Start Power/MDD Simulation (statL)\n")
cat(sprintf("Config: m=%d | Kernel: %s | Bandwidth: 1.0*d\n", m_val, kernal))
cat("==========================================================\n\n")

# --- 3. 启动并行 ---
num_cores <- parallel::detectCores(logical = FALSE) - 1
if(num_cores < 1) num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')
    library(mvtnorm); library(Rcpp); library(RcppArmadillo); library(MTS)
    source("Gendata4PaperNew.R")
    sourceCpp("statLMul.cpp") 
})

# --- 4. 主循环 ---
final_results_all <- data.frame()

for (mod_id in model_list) {
    cat(sprintf("\n>>> Processing Model %d ...\n", mod_id))
    
    # ---------------------------------------------------------
    # Phase 1: 阈值校准 (H0)
    # ---------------------------------------------------------
    cat("   [Phase 1] Calibrating Threshold (Dynamic Delta)...\n")
    null_max_stats <- foreach(s = 1:nsim, .combine = c) %dopar% {
        # 生成无变点数据
        yt <- gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod_id, delta = 0)
        if(!is.matrix(yt)) yt <- as.matrix(yt)
        
        # [NEW] 动态计算带宽参数
        d_dim <- ncol(yt)
        delta_dynamic <- 1.0 * d_dim
        
        # 计算统计量
        stats <- statL(yt, m_val, m0, alphaE, vaha_vec, kernal, delta_dynamic)
        max(stats)
    }
    Threshold <- quantile(null_max_stats, 1 - TypeI, na.rm=TRUE)
    cat(sprintf("   -> Threshold (95%%): %.4f\n", Threshold))
    
    # ---------------------------------------------------------
    # Phase 2: 计算 Power & MDD (H1)
    # ---------------------------------------------------------
    cat("   [Phase 2] Evaluating Deltas...\n")
    
    for (d_val in deltas) {
        
        # 返回矩阵列: [Is_Rejected, Valid_Delay]
        res_mat <- foreach(s = 1:nsim, .combine = rbind) %dopar% {
            
            yt <- gendataPaper(monitoringPeriod = L_val, trainSize = m_val, model = mod_id, 
                               cptLocation = 0.2, delta = d_val)
            if(!is.matrix(yt)) yt <- as.matrix(yt)
            
            # 核函数参数 1*dim
            d_dim <- ncol(yt)
            delta_dynamic <- 1.0 * d_dim
            
            stats <- statL(yt, m_val, m0, alphaE, vaha_vec, kernal, delta_dynamic)
            alarm_indices <- which(stats > Threshold)
            
            if (length(alarm_indices) > 0) {
                # statL 返回的索引 k 对应绝对时刻 m + k
                first_alarm_k <- alarm_indices[1]
                alarm_time <- m_val + first_alarm_k
                
                # 计算延迟 = 报警时刻 - 真实变点时刻
                delay <- alarm_time - true_cp_time
                
                # 仅记录有效延迟 (报警发生在变点之后)
                # 误报 (delay <= 0) 记为 NA，不计入 MDD
                valid_delay <- if(delay > 0) delay else NA
                
                return(c(1, valid_delay))
            } else {
                # 未拒绝
                return(c(0, NA))
            }
        }
        
        # 1. Power (拒绝率)
        rejection_rate <- mean(res_mat[, 1]) * 100
        
        # 2. MDD (排除 NA)
        valid_delays <- res_mat[, 2]
        mdd <- mean(valid_delays[!is.na(valid_delays)], na.rm=TRUE)
        if(is.nan(mdd)) mdd <- NA
        
        cat(sprintf("      Delta=%.2f | Power: %.1f%% | MDD: %.2f\n", d_val, rejection_rate, mdd))
        
        final_results_all <- rbind(final_results_all, data.frame(
            Model = mod_id,
            Delta = d_val,
            Method = paste0("statL_", kernal),
            Power = rejection_rate,
            MDD = mdd
        ))
    }
}

stopCluster(cl)

# --- 5. 保存与绘图 ---
cat("\n========================================\n")
cat("Simulation Finished.\n")
print(head(final_results_all))

# 保存文件名增加 adaptive 标识
fileSaveName <- paste0("Power_MDD_", kernal, "(", head(model_list, 1), "-", tail(model_list, 1), ")")
write.csv(final_results_all, file = paste0(fileSaveName, ".csv"), row.names = FALSE)

# 简单绘图
if(nrow(final_results_all) > 0) {
    p1 <- ggplot(final_results_all, aes(x=Delta, y=Power, color=factor(Model))) +
        geom_line() + geom_point() +
        labs(title=paste0("Power (", kernal, ", delta=1.0*d)"), y="Power (%)", color="Model") +
        theme_minimal()
    
    p2 <- ggplot(final_results_all, aes(x=Delta, y=MDD, color=factor(Model))) +
        geom_line() + geom_point() +
        labs(title=paste0("MDD (", kernal, ", delta=1.0*d)"), y="Mean Detection Delay", color="Model") +
        theme_minimal()
    
    print(p1)
    print(p2)
    ggsave(paste0("statL_", kernal, "(", head(model_list, 1), "-", tail(model_list, 1), ")", "_Power_adaptive.png"), p1, width=6, height=4)
    ggsave(paste0("statL_", kernal, "(", head(model_list, 1), "-", tail(model_list, 1), ")", "_MDD_adaptive.png"), p2, width=6, height=4)
}
