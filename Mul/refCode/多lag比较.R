rm(list = ls())
### 多Lag参数比较: Type I Error or Power ###
#########################################################

# --- 0. 环境与依赖 ---
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(compiler)     # 字节编译
library(progress)     # 进度条
library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)

# 加载源码 (确保 Gendata4PaperNew.R 是之前修复过无Bug的版本)
sourceCpp('statLMul.cpp')
sourceCpp('eucLMulNew.cpp')
source("Gendata4PaperNew.R") 

# --- 1. 全局参数设置 ---

# Boundary function
Qk <- function(s0, gamma0){
  (1 + s0) * (s0 / (1 + s0))^gamma0
}

# Lag window function
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}

# 编译关键函数
lambda_compiled <- cmpfun(lambda)
statL_compiled <- cmpfun(statL)
Qk_compiled <- cmpfun(Qk)

# --- 模型参数 ---
DGP <- 101           # 101=Mean Change (Power), 1=Null (Type I)
lags <- c(0, 1, 2, 5, 10, 20, 30) # 待测试的 lags (注意 m < T)
trainSize <- 100     # T
monitoringPeriod <- 1 # L
gamma <- 0
tunta <- 1           # 1=gauss
TypeI <- 0.05
nsim <- 1000         # 模拟次数

# --- 关键修复：动态探测维度 d ---
# 在循环开始前，生成一次数据来确定 d，防止硬编码错误
cat("正在初始化参数与探测数据维度...\n")
temp_yt <- gendataPaper(monitoringPeriod, trainSize, model = DGP)
d <- ncol(temp_yt)
cat(sprintf("检测到模型 DGP=%d 的维度 d=%d\n", DGP, d))

# Kernel 参数
delta0 <- 2 * d      # Heuristic parameter for Kernel
if (tunta == 1) { kernal <- "gauss" } else 
  if (tunta == 2) { kernal <- "energy" } else 
    if (tunta == 3) { kernal <- "quadexp" }

# Bootstrap 参数预计算
KT <- max(5, sqrt(log10(trainSize))) 
sqrt_log_T <- 2 * sqrt(log10(trainSize)/trainSize)

# --- 2. 核心模拟函数 (针对单个 lag) ---
run_simulation_for_lag <- function(curr_lag, nsim, trainSize, monitoringPeriod, DGP, d, delta0, kernal) {
  
  # 【关键】：Lag 改变导致 Tm 改变，权重 vaha 必须在函数内重算
  Tm <- trainSize - curr_lag
  
  # 安全检查：防止 lag 过大导致 Tm <= 0
  if (Tm <= 10) {
    warning(paste("Lag", curr_lag, "is too large for trainSize", trainSize, ". Skipping."))
    return(NA)
  }
  
  N <- (1 + monitoringPeriod) * trainSize
  k_vec <- 1:(N - trainSize)
  # 使用编译后的 Qk 函数
  vaha <- (Tm + k_vec)^2 / (Tm * Qk_compiled(k_vec/Tm, gamma)^2)
  alphaE <- 1
  
  # 容器
  Tn <- numeric(nsim)
  Tn.star <- numeric(nsim)
  
  # 进度条
  pb <- progress_bar$new(
    format = paste0(" Lag ", curr_lag, " [:bar] :percent | ETA: :eta"),
    total = nsim, clear = FALSE, width = 60
  )
  
  for (s in 1:nsim) {
    # (A) 生成数据
    yt <- gendataPaper(monitoringPeriod = monitoringPeriod, trainSize = trainSize, model = DGP)
    
    # (B) 自动选择 Block Size (向量化)
    Block_sizes <- sapply(1:d, function(i) {
      ut <- yt[1:trainSize, i]
      # 【修复】：lag.max 必须为整数
      acf_out <- acf(ut, lag.max=as.integer(2*KT), plot=FALSE, demean=TRUE)
      R <- drop(acf_out$acf)
      
      # 【修复】：索引必须为整数
      tmp <- which(abs(R[1:as.integer(KT)]) < sqrt_log_T)
      M <- if (length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
      
      M_range <- -M:M
      R_vals <- R[abs(M_range) + 1]
      lam_vals <- lambda_compiled(M_range/M)
      
      ghat <- sum(lam_vals * R_vals)
      Ghat <- sum(lam_vals * R_vals * abs(M_range))
      D.SB <- 2 * ghat^2
      
      if(is.na(D.SB) || D.SB <= 1e-8) return(1)
      max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainSize^(1/3)))
    })
    blocksize <- max(1, round(mean(Block_sizes)))
    
    # (C) Stationary Bootstrap (向量化索引)
    expected_blocks <- ceiling(N / blocksize) * 2 + 5
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
    
    # 【修复】：drop=FALSE 保护矩阵维度
    y.star <- yt[bootstrap_idx, , drop=FALSE]
    
    # (D) 计算统计量 (使用传入的 curr_lag)
    Tn[s] <- max(statL_compiled(yt, trainSize, m=curr_lag, alphaE, vaha, kernal, delta0))
    Tn.star[s] <- max(statL_compiled(y.star, trainSize, m=curr_lag, alphaE, vaha, kernal, delta0))
    
    pb$tick()
  }
  
  # 计算拒绝率
  final_crit <- quantile(Tn.star, 1 - TypeI, na.rm=TRUE)
  final_rej <- mean(Tn > final_crit, na.rm=TRUE)
  
  return(final_rej)
}

# --- 3. 主执行循环 ---

cat(sprintf("\nStarting Simulation for %d lags: %s\n", length(lags), paste(lags, collapse=",")))
cat(sprintf("Model: %d | T: %d | Kernel: %s | d: %d\n", DGP, trainSize, kernal, d))

results_vec <- numeric(length(lags))
names(results_vec) <- paste0("lag", lags)

total_start_time <- Sys.time()

for (i in seq_along(lags)) {
  curr_lag <- lags[i]
  
  # 运行单个 lag 的模拟
  res <- run_simulation_for_lag(curr_lag, nsim, trainSize, monitoringPeriod, DGP, d, delta0, kernal)
  
  # 记录结果 (处理 NA 情况)
  if (!is.na(res)) {
    results_vec[i] <- round(res * 100, 2)
    cat(sprintf("\n[Lag %2d] Finished. Rejection Rate: %.2f%%\n", curr_lag, results_vec[i]))
  } else {
    results_vec[i] <- NA
    cat(sprintf("\n[Lag %2d] Skipped (Lag too large for T).\n", curr_lag))
  }
  
  # 内存清理
  gc(verbose = FALSE)
}

total_end_time <- Sys.time()
cat("\nTotal Time:", round(difftime(total_end_time, total_start_time, units="mins"), 2), "minutes.\n")

# --- 4. 结果处理与可视化 ---

print(results_vec)

# 保存数据框
results_df <- data.frame(
  Lag = factor(lags, levels = lags), # 保持顺序
  RejectionRate = results_vec
)

# 保存 CSV
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
filename <- paste0("SimResult_LagComp_DGP", DGP, "_", timestamp, ".csv")
write.csv(results_df, filename, row.names = FALSE)
cat("Results saved to:", filename, "\n")

# 绘图
# 根据 DGP 判断是 Power 还是 Type I Error，设置不同的标线
y_label <- if (DGP <= 8 || DGP == 306) "Type I Error (%)" else "Power (%)"
ref_line <- if (DGP <= 8 || DGP == 306) 5 else NA # 仅在 Type I 时画 5% 线

p <- ggplot(results_df, aes(x = Lag, y = RejectionRate)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = RejectionRate), vjust = -0.5, size = 3.5) +
  labs(title = paste("Simulation Result by Lag (Model:", DGP, ")"),
       subtitle = paste("T =", trainSize, ", L =", monitoringPeriod, ", Kernel =", kernal),
       x = "Lag Parameter (m)", 
       y = y_label) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# 只有在测试 Type I Error 时才画 5% 参考线
if (!is.na(ref_line)) {
  p <- p + geom_hline(yintercept = ref_line, linetype = "dashed", color = "red")
}

print(p)