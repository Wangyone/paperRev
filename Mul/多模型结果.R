rm(list = ls())
######################### Multi-Model Comparison (Auto-Dimension) #########################

# --- 1. 环境与依赖 ---
setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

library(compiler)
library(progress)
library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)

sourceCpp('statLMul.cpp')
sourceCpp('eucLMulNew.cpp')
source("Gendata4PaperNew.R") # 这里会加载上面更新过的文件

# --- 2. 辅助函数 ---
Qk <- function(s0, gamma0){ (1 + s0) * (s0 / (1 + s0))^gamma0 }
lambda <- function(t) { abst <- abs(t); (abst <= .5) + 2 * (1 - abst) * ((abst > .5) & (abst <= 1)) }
lambda_compiled <- cmpfun(lambda)
Qk_compiled <- cmpfun(Qk)
statL_compiled <- cmpfun(statL)

# --- 3. 参数设置 ---

# 模型列表 (Type I Error 测试: DGP 1-5 c(1,2, 3, 4, 5) )
# 如果要测 Power，可以改为 c(101,102, 103, 104)  等
DGP_list <- c(1,2, 3, 4, 5)

kernal <- "gauss"    
trainingSize <- 200  
monitoringPeriod <- 1 
gamma <- 0           
nsim <- 1000         
TypeI <- 0.05        
lag <- 0             

Tm <- trainingSize - lag
N <- (1 + monitoringPeriod) * trainingSize
k_vec <- 1:(N - trainingSize)
vaha <- (Tm + k_vec)^2 / (Tm * Qk_compiled(k_vec/Tm, gamma)^2)
alphaE <- 1
KT <- max(5, sqrt(log10(trainingSize)))
log_T_limit <- 2 * sqrt(log10(trainingSize) / trainingSize)

cat("========================================\n")
cat("开始模拟 (自动检测维数)...\n")
cat("T:", trainingSize, "| Sim:", nsim, "\n")
cat("Models:", paste(DGP_list, collapse=", "), "\n")
cat("========================================\n\n")

final_results <- data.frame(
  Model = DGP_list,
  Dim = integer(length(DGP_list)), 
  RejectionRate = numeric(length(DGP_list))
)

for (i in seq_along(DGP_list)) {
  curr_model <- DGP_list[i]
  
  Tn <- numeric(nsim)
  Tn.star <- numeric(nsim)
  detected_d <- NA
  
  pb <- progress_bar$new(
    format = paste0("  Model ", curr_model, " [:bar] :percent | Dim: :dim | Rej: :rej"),
    total = nsim, clear = FALSE, width = 70
  )
  
  for (s in 1:nsim) {
    # 1. 生成数据
    yt <- gendataPaper(monitoringPeriod, trainingSize, model = curr_model)
    
    # 2. 自动检测维数
    if (!is.matrix(yt)) yt <- as.matrix(yt)
    d <- ncol(yt)
    
    if (s == 1) {
      detected_d <- d
      if (kernal == "gauss" || kernal == "energy") {
        curr_delta0 <- 1 
        # 暂定为1，后续我手动调整2 * d
      } else {
        E <- eucL(yt[1:trainingSize, , drop=FALSE], lag)
        curr_delta0 <- median(E[E > 0])
      }
    }
    
    # 3. Block Size (使用检测到的 d)
    Block_sizes <- sapply(1:d, function(dim_i) {
      ut <- yt[1:trainingSize, dim_i]
      R <- as.vector(acf(ut, lag.max = 2*KT, plot = FALSE)$acf)
      tmp <- which(abs(R[1:KT]) < log_T_limit)
      M <- if (length(tmp) > 0) 2 * (tmp[1] - 1) else 2 * (KT - 1)
      idx <- -M:M
      vals <- R[abs(idx) + 1]
      w <- lambda_compiled(idx / M)
      ghat <- sum(w * vals)
      Ghat <- sum(w * vals * abs(idx))
      D.SB <- 2 * ghat^2
      if (D.SB < 1e-10) return(1)
      max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainingSize^(1/3)))
    })
    blocksize <- max(1, round(mean(Block_sizes)))
    
    # 4. Stationary Bootstrap
    expected_blocks <- ceiling((N + trainingSize) / blocksize) + 5
    len_vec <- rgeom(expected_blocks, 1/blocksize) + 1
    while(sum(len_vec) < N) len_vec <- c(len_vec, rgeom(10, 1/blocksize) + 1)
    
    cumsum_len <- cumsum(len_vec)
    cutoff <- which(cumsum_len >= N)[1]
    len_vec <- len_vec[1:cutoff]
    num_blocks <- length(len_vec)
    start_vec <- sample.int(trainingSize, num_blocks, replace = TRUE)
    
    idx_list <- lapply(1:num_blocks, function(j) {
      (start_vec[j] + 0:(len_vec[j]-1) - 1) %% trainingSize + 1
    })
    bootstrap_idx <- unlist(idx_list)[1:N]
    y.star <- yt[bootstrap_idx, , drop = FALSE]
    
    # 5. 计算统计量
    Tn[s] <- max(statL_compiled(yt, trainingSize, m = lag, alphaE, vaha, kernal, curr_delta0))
    Tn.star[s] <- max(statL_compiled(y.star, trainingSize, m = lag, alphaE, vaha, kernal, curr_delta0))
    
    if (s %% 50 == 0 && s > 1) {
      crit <- quantile(Tn.star[1:s], 1 - TypeI, na.rm=TRUE)
      rej_rate <- mean(Tn[1:s] > crit) * 100
      curr_rej_str <- sprintf("%.1f%%", rej_rate)
    } else {
      curr_rej_str <- "..."
    }
    pb$tick(tokens = list(dim = detected_d, rej = curr_rej_str))
  }
  
  final_crit <- quantile(Tn.star, 1 - TypeI)
  final_rej <- mean(Tn > final_crit) * 100
  
  final_results$Dim[i] <- detected_d
  final_results$RejectionRate[i] <- final_rej
  
  cat(sprintf("\n[Model %d] Dim=%d | Rejection Rate: %.2f%%\n", curr_model, detected_d, final_rej))
}

cat("\n==============================\n")
print(final_results)

# 实现 多模型+trainingSize={100,200}，monitoringPeriod ={1,2,3}的结果模拟，
# 结果按照以下样式保持为csv：
# \toprule
# & $m$   & $L$ & $N_1$   & $N_2$   & $N_3$   & $N_4$   & $N_5$  
#   & $N_6$    & $N_7$   & $N_8$    & $N_9$    \\ \midrule
