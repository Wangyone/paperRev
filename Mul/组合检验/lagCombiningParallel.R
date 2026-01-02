rm(list = ls())
# setwd("D:/Desktop/paperRev/JointCUSUM/Mul/") # 设置路径

library(compiler)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(foreach)
library(doParallel)

# --- 1. 参数设置 ---

# 边界函数
Qk <- function(s0, gamma0){
    (1 + s0) * (s0 / (1 + s0))^gamma0
}

DGP <- 306           
lags <- c(0, 1, 2) 
nlag <- length(lags)
tunta <- 1         
delta0 <- 1
gamma <- 0
trainSize <- 100   
monitoringPeriod <- 1 
TypeI <- 0.05      
alphaE <- 1

nsim <- 800        
M <- 299           

# --- 2. 预计算 ---

N <- (1 + monitoringPeriod) * trainSize
kernal <- switch(tunta, "1"="gauss", "2"="energy", "3"="quadexp")

weight_list <- list()
for(lag_val in lags) {
    Tm_lag <- trainSize - lag_val
    k_vec <- 1:(N - trainSize)
    w <- (Tm_lag + k_vec)^2 / (Tm_lag * Qk(k_vec/Tm_lag, gamma)^2)
    weight_list[[as.character(lag_val)]] <- w
}

KT <- max(5, sqrt(log10(trainSize)))  
sqrt_log_T <- 2 * sqrt(log10(trainSize)/trainSize) 

calc_fisher <- function(pvals) {
    -2 * sum(log(pmax(pvals, 1e-10)))
}

# --- 3. 启动并行 ---

num_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("已启动并行计算，使用核心数: %d\n", num_cores))

# 使用 clusterEvalQ 配置子进程环境
clusterEvalQ(cl, {
    library(mvtnorm)
    library(Rcpp)
    library(RcppArmadillo)
    library(MTS)
    library(compiler)
    
    # 设置子进程路径 (必须修改为你的实际路径)
    setwd("D:/Desktop/paperRev/JointCUSUM/Mul/") 
    
    # 重新编译函数
    sourceCpp('statLMul.cpp')
    sourceCpp('eucLMulNew.cpp')
    if(file.exists("Gendata4PaperNew.R")) source("Gendata4PaperNew.R")
    
    lambda <- function(t) {
        abst <- abs(t)
        (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
    }
    lambda_compiled <- cmpfun(lambda)
    statL_compiled <- cmpfun(statL) 
})

# --- 4. 主循环 ---

cat("开始模拟 (Robust Parallel Mode)...\n")
startTime <- Sys.time()

# 使用 .errorhandling = "pass" 或 "remove" 可以在单个任务失败时不中断整个程序
# 这里我们用默认，但在内部做好异常捕获
final_results <- foreach(s = 1:nsim, .combine = rbind, 
                         # 移除 .export 中重复的变量以消除警告
                         .export = c("weight_list", "lags", "nlag", "Qk",
                                     "trainSize", "monitoringPeriod", "DGP", "N", "alphaE", 
                                     "kernal", "delta0", "M", "KT", "sqrt_log_T", "calc_fisher")) %dopar% {
                                         
                                         # === [Safety Check 1] 数据生成 ===
                                         yt <- tryCatch({
                                             gendataPaper(monitoringPeriod, trainSize, DGP)
                                         }, error = function(e) return(NULL))
                                         
                                         # 如果生成失败或包含 NA/Inf，返回全 NA 结果，跳过此轮
                                         if(is.null(yt) || any(!is.finite(yt))) {
                                             return(rep(NA, nlag + 1))
                                         }
                                         
                                         d <- ncol(yt)
                                         
                                         # 计算原始统计量
                                         T_orig <- numeric(nlag)
                                         for (j in 1:nlag) {
                                             lag_val <- lags[j]
                                             seq_val <- statL_compiled(yt, trainSize, m=lag_val, alphaE, 
                                                                       weight_list[[as.character(lag_val)]], kernal, delta0)
                                             # === [Safety Check 2] 统计量 ===
                                             # 如果统计量算出 NA (极少见)，设为 -Inf 防止比较出错
                                             val <- max(seq_val)
                                             T_orig[j] <- if(is.finite(val)) val else -Inf 
                                         }
                                         
                                         # === [Safety Check 3] Block Size 计算 (最容易报错的地方) ===
                                         blocksize <- 1 # 默认值
                                         try({
                                             ut <- yt[1:trainSize, 1]
                                             # 检查 variance 是否为 0 (导致 acf 失败)
                                             if(var(ut) > 1e-10) {
                                                 acf_out <- acf(ut, lag.max=2*KT, plot=FALSE, demean=TRUE) 
                                                 R <- drop(acf_out$acf)
                                                 
                                                 # 确保 R 中没有 NA
                                                 if(all(is.finite(R))) {
                                                     tmp <- which(abs(R[1:KT]) < sqrt_log_T)
                                                     M_lag <- if (length(tmp) > 0) 2*(tmp[1] - 1) else 2*(KT - 1)
                                                     M_range <- -M_lag:M_lag
                                                     
                                                     lam_vals <- lambda_compiled(M_range/M_lag)
                                                     ghat <- sum(lam_vals * R[abs(M_range) + 1])
                                                     Ghat <- sum(lam_vals * R[abs(M_range) + 1] * abs(M_range))
                                                     D.SB <- 2 * ghat^2
                                                     
                                                     # 这里的 if 判断最容易报 "NA/NaN argument"
                                                     # 我们先判断 D.SB 是否为有效数值
                                                     if(!is.na(D.SB) && D.SB > 1e-8) {
                                                         blocksize <- max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainSize^(1/3)))
                                                     }
                                                 }
                                             }
                                         }, silent = TRUE)
                                         
                                         # 防止 blocksize 为 NA
                                         if(is.na(blocksize) || blocksize < 1) blocksize <- 5 
                                         
                                         # Bootstrap 循环
                                         T_boot <- matrix(NA, nrow = M, ncol = nlag)
                                         
                                         # 预估需要的块数量
                                         num_blocks_est <- ceiling(N / blocksize) + 10
                                         
                                         for (b in 1:M) {
                                             # 生成几何分布长度
                                             # 1/blocksize 必须有效，blocksize >= 1 保证了这点
                                             len_vec <- rgeom(num_blocks_est, 1/blocksize) + 1
                                             cumsum_len <- cumsum(len_vec)
                                             
                                             # 找到截断点
                                             cut_idx_vec <- which(cumsum_len >= N)
                                             if(length(cut_idx_vec) == 0) {
                                                 # 极罕见情况：生成的长度不够，补救一下
                                                 cut_point <- length(len_vec)
                                             } else {
                                                 cut_point <- cut_idx_vec[1]
                                             }
                                             valid_lens <- len_vec[1:cut_point]
                                             
                                             start_points <- sample.int(trainSize, length(valid_lens), replace = TRUE)
                                             
                                             idx_list <- lapply(seq_along(valid_lens), function(k) {
                                                 (start_points[k] + 0:(valid_lens[k]-1) - 1) %% trainSize + 1
                                             })
                                             idx_vec <- unlist(idx_list)[1:N]
                                             
                                             y_star <- yt[idx_vec, , drop=FALSE]
                                             
                                             for (j in 1:nlag) {
                                                 lag_val <- lags[j]
                                                 seq_star <- statL_compiled(y_star, trainSize, m=lag_val, alphaE, 
                                                                            weight_list[[as.character(lag_val)]], kernal, delta0)
                                                 val <- max(seq_star)
                                                 T_boot[b, j] <- if(is.finite(val)) val else -Inf
                                             }
                                         }
                                         
                                         # P 值计算
                                         T_total <- rbind(T_orig, T_boot)
                                         p_mat <- matrix(0, nrow = M + 1, ncol = nlag)
                                         
                                         for (j in 1:nlag) {
                                             # 这里的 T_total 可能包含 -Inf，rank 处理没问题
                                             ranks <- rank(-T_total[, j], ties.method = "max", na.last = TRUE)
                                             p_mat[, j] <- (ranks + 1) / (M + 2)
                                         }
                                         
                                         # Fisher 组合
                                         W_vec <- apply(p_mat, 1, calc_fisher)
                                         W_orig <- W_vec[1]
                                         W_boot <- W_vec[2:(M+1)]
                                         
                                         # 决策
                                         p_combined <- (sum(W_boot >= W_orig) + 1) / (M + 1)
                                         res_combined <- as.numeric(p_combined < TypeI)
                                         res_individual <- as.numeric(p_mat[1, ] < TypeI)
                                         
                                         c(res_individual, res_combined)
                                     }

# --- 5. 结果处理 ---

stopCluster(cl)

# 移除可能产生的 NA 行 (由数据生成失败导致)
final_results <- na.omit(final_results)

elapsed <- difftime(Sys.time(), startTime, units="secs")
cat(sprintf("\n模拟完成! 有效样本数: %d/%d, 总耗时: %.2f 秒\n", nrow(final_results), nsim, as.numeric(elapsed)))

colnames(final_results) <- c(paste0("Lag_", lags), "Combined")
rejection_rates <- colMeans(final_results) * 100

res_df <- data.frame(
    Test = names(rejection_rates),
    Rejection_Rate_Percent = rejection_rates
)
print(res_df)