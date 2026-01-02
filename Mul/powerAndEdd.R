rm(list = ls())
######################### 模拟：拒绝率 (Power) 与 平均检测延迟 (EDD) #########################
setwd("D:/Desktop/paperRev/JointCUSUM/Mul/")

# --- 0. 环境与依赖 ---
library(compiler)
library(progress)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(latex2exp) 

# 加载源码 (请确保路径正确)
sourceCpp('statLMul.cpp')
sourceCpp('eucLMulNew.cpp')
source("Gendata4PaperNew.R")

# --- 1. 参数设置 ---

# 边界函数
Qk <- function(s0, gamma0){
    (1 + s0) * (s0 / (1 + s0))^gamma0
}

# 基础参数
DGP <- 3          # 模型 ID
m0 <- 2           # 滞后阶数 (Lag)
tunta <- 1        # 核函数类型: 1=高斯 (gauss)
delta0 <- 1       # 核参数
gamma <- 0        # 边界函数参数
trainSize <- 200  # 训练样本大小 (T)
monitoringPeriod <- 1 # 监控期倍数因子 (L)，总监控长度 = L * T
TypeI <- 0.05     # 显著性水平 (Alpha)

# 变点设置 (Change Point Settings)
tau <- 0.2          # 变点发生位置 (相对于监控期开始的比例)

# 计算变点在监控期内的相对索引 (k*)
# 监控长度 = 1 * 200 = 200
# 变点相对位置 k* = 200 * 0.2 = 40 (即变点发生在进入监控期后的第 40 个时刻)
monitor_len <- trainSize * monitoringPeriod
k_star <- floor(monitor_len * tau) 

# 模拟序列 delta (变化幅度)
delta_seq <- seq(0, 0.4, by = 0.2) 
nsim <- 1000      # 模拟次数

# --- 2. 预计算与函数编译 ---

N <- (1 + monitoringPeriod) * trainSize # 总样本量
Tm <- trainSize
k_vec <- 1:(N - trainSize)
# 计算权重向量
vaha <- (Tm + k_vec)^2 / (Tm * Qk(k_vec/Tm, gamma)^2) 
alphaE <- 1
# 确定核函数名称
kernal <- switch(tunta, "1"="gauss", "2"="energy", "3"="quadexp")

# Bootstrap 参数
KT <- max(5, sqrt(log10(trainSize)))  
sqrt_log_T <- 2 * sqrt(log10(trainSize)/trainSize) 

# 平顶滞后窗口函数 (Flat-top lag-window function)
lambda <- function(t) {
    abst <- abs(t)
    (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}

# 编译函数以提升运行速度
lambda_compiled <- cmpfun(lambda)
statL_compiled <- cmpfun(statL) 

# --- 数据生成包装函数 (根据 delta 和 k* 注入变点) ---
genData_Sim <- function(T_size, L_factor, model_id, delta_mag, k_star_idx) {
    # 1. 生成基础 H0 数据 (无变点)
    yt <- gendataPaper(monitoringPeriod = L_factor, trainSize = T_size, model = model_id)
    
    # 2. 注入变点 (均值漂移)
    if (delta_mag != 0) {
        N_total <- nrow(yt)
        # 计算变点的绝对索引：训练集大小 + 监控期相对索引 + 1
        cp_abs_index <- T_size + k_star_idx + 1
        
        if (cp_abs_index <= N_total) {
            # 从变点开始直到结束，加上偏移量 delta
            yt[cp_abs_index:N_total, ] <- yt[cp_abs_index:N_total, ] + delta_mag
        }
    }
    return(yt)
}

# --- 3. 主循环：遍历 Delta ---

results_df <- data.frame(
    Delta = delta_seq,
    Power = NA,
    EDD = NA
)

cat("================================================\n")
cat("开始模拟 (Starting Simulation)\n")
cat(sprintf("变点相对位置 k* (Relative Change Point): %d\n", k_star))
cat("================================================\n")

for (idx_d in 1:length(delta_seq)) {
    
    curr_delta <- delta_seq[idx_d]
    
    # 初始化存储容器
    Tn_Stat_List <- vector("list", nsim) # 存储每次模拟的完整统计量序列
    Tn_max <- numeric(nsim)              # 观测数据的最大统计量
    Tn_star_max <- numeric(nsim)         # Bootstrap 数据的最大统计量
    
    # 进度条设置
    pb <- progress_bar$new(
        format = sprintf(" Delta %.2f [:bar] :percent | ETA: :eta", curr_delta),
        total = nsim, clear = FALSE, width = 60
    )
    
    for (s in 1:nsim) {
        
        # (A) 生成在 k_star 处发生变点的数据
        yt <- genData_Sim(trainSize, monitoringPeriod, DGP, curr_delta, k_star)
        d <- ncol(yt)
        
        # (B) 自动选择 Block Size (保持原有逻辑)
        Block_sizes <- sapply(1:d, function(i) {
            ut <- yt[1:trainSize, i]
            acf_out <- acf(ut, lag.max=2*KT, plot=FALSE, demean=TRUE) 
            R <- drop(acf_out$acf)
            tmp <- which(abs(R[1:KT]) < sqrt_log_T)
            M <- if (length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
            M_range <- -M:M
            R_vals <- R[abs(M_range) + 1]
            lam_vals <- lambda_compiled(M_range/M)
            ghat <- sum(lam_vals * R_vals)
            Ghat <- sum(lam_vals * R_vals * abs(M_range))
            D.SB <- 2 * ghat^2
            if(D.SB <= 1e-8) return(1) 
            max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainSize^(1/3)))
        })
        blocksize <- max(1, round(mean(Block_sizes)))
        
        # (C) 平稳 Bootstrap (Stationary Bootstrap)
        expected_blocks <- ceiling(N / blocksize) * 2
        len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
        while(sum(len_vec) < N) len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
        cumsum_len <- cumsum(len_vec)
        cutoff_idx <- which(cumsum_len >= N)[1]
        len_vec <- len_vec[1:cutoff_idx]
        num_blocks <- length(len_vec)
        start_vec <- sample.int(trainSize, num_blocks, replace = TRUE)
        idx_list <- lapply(1:num_blocks, function(j) {
            (start_vec[j] + 0:(len_vec[j]-1) - 1) %% trainSize + 1
        })
        bootstrap_idx <- unlist(idx_list)[1:N]
        y.star <- yt[bootstrap_idx, , drop=FALSE]
        
        # (D) 计算统计量
        # statL 返回监控期长度 (N-T) 的向量
        # 索引 1 对应 T+1, 索引 k 对应 T+k
        stats_seq <- statL_compiled(yt, trainSize, m=m0, alphaE, vaha, kernal, delta0)
        Tn_Stat_List[[s]] <- stats_seq
        Tn_max[s] <- max(stats_seq)
        
        # Bootstrap 统计量 (仅需最大值用于定阈值)
        star_seq <- statL_compiled(y.star, trainSize, m=m0, alphaE, vaha, kernal, delta0)
        Tn_star_max[s] <- max(star_seq)
        
        pb$tick()
    }
    
    # (E) 计算指标 (Metrics)
    
    # 1. 确定临界值 (Critical Value)
    crit_val <- quantile(Tn_star_max, 1 - TypeI, na.rm=TRUE)
    
    # 2. 计算拒绝率 (Power)
    is_rejected <- Tn_max > crit_val
    power_est <- mean(is_rejected)
    
    # --- [修正] EDD (平均检测延迟) 计算逻辑 ---
    rejected_indices <- which(is_rejected) # 获取所有拒绝 H0 的实验索引
    
    if (length(rejected_indices) > 0) {
        delays <- sapply(rejected_indices, function(idx) {
            seq_vals <- Tn_Stat_List[[idx]]
            
            # 找到统计量首次超过临界值的索引 (first_crossing)
            # 注意：此索引是相对于监控期开始的 (例如 1 代表 T+1)
            first_crossing <- which(seq_vals > crit_val)[1]
            
            # 逻辑说明：
            # 1. 如果报警发生在变点之后 (first_crossing > k_star)：
            #    延迟 = 报警时间 - 变点时间
            # 2. 如果报警发生在变点之前或刚好在变点 (first_crossing <= k_star)：
            #    这属于在变点发生前的“误报”(False Alarm)，不应计入检测延迟。
            #    返回 NA 将其剔除。
            
            if (first_crossing > k_star) {
                return(first_crossing - k_star)
            } else {
                return(NA) # 误报，剔除
            }
        })
        # 计算平均值，自动忽略 NA (即忽略误报)
        edd_est <- mean(na.omit(delays))
    } else {
        edd_est <- NA # 无任何拒绝
    }
    
    # 存储结果
    results_df$Power[idx_d] <- power_est * 100 
    results_df$EDD[idx_d] <- edd_est
    
    # 打印当前 Delta 的结果
    cat(sprintf("\n结果: Delta=%.1f | Power=%.1f%% | EDD=%.2f\n", 
                curr_delta, power_est * 100, ifelse(is.na(edd_est), 0, edd_est)))
    
    # 清理内存
    rm(Tn_Stat_List, Tn_max, Tn_star_max)
    gc(verbose = FALSE)
}

# --- 4. 绘图 (Plotting) ---

print(results_df)

par(mfrow = c(1, 2)) # 设置一行两图

# 图 1: Power 曲线
plot(results_df$Delta, results_df$Power, 
     type = "o", pch = 19, col = "blue", lwd = 2,
     ylim = c(0, 100),
     xlab = TeX("变化幅度 ($\\delta$)"),
     ylab = "拒绝率 (%)",
     main = "Power Curve (功效曲线)")
abline(h = 5, lty = 2, col = "red") # 5% 水平线
grid()

# 图 2: EDD 曲线
plot_df <- na.omit(results_df) # 绘图时移除 NA 数据
# 动态设置 Y 轴上限
y_max <- if(nrow(plot_df) > 0) max(plot_df$EDD) * 1.2 else 10

plot(plot_df$Delta, plot_df$EDD, 
     type = "o", pch = 19, col = "darkgreen", lwd = 2,
     ylim = c(0, y_max),
     xlab = TeX("变化幅度 ($\\delta$)"),
     ylab = "平均检测延迟 (EDD)",
     main = "EDD Curve (延迟曲线)")
grid()

par(mfrow = c(1, 1)) # 恢复默认画板