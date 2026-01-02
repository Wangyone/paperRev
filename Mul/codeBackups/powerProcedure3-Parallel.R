rm(list = ls())
library(foreach)
library(doParallel)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
library(compiler)

setwd('D:/Desktop/paperRev/JointCUSUM/Mul/')

# --- 主进程预计算参数 ---
# 仅用于计算 delta0 和 vaha 等常量
Rcpp::sourceCpp('statLMul.cpp')
Rcpp::sourceCpp("eucLMulNew.cpp")
source('dataGenerationNew.R')

d <- 2 
lag <- 1
kernelTunta <- 3  
trainingSize <- 100
monitoringPeriod <- 1
gamma <- 0
deltagauss <- 1*d
nsim <- 1000
DGP <- 301
TypeI <- 0.05
alphaE <- 1

# 初始化参数
yt <- dataGeneration(monitoringPeriod, trainingSize, cptModel = DGP, dataDimension = d)
E <- eucL(yt[1:trainingSize,], lag)
medianE <- (median(E[E > 0]))
rm(E)

if (kernelTunta == 1) { kernal <- "gauss"; delta0 <- deltagauss } 
if (kernelTunta == 2) { kernal <- "energy"; delta0 <- deltagauss } 
if (kernelTunta == 3) { kernal <- "quadexp"; delta0 <- medianE }

Qk <- function(s0,gamma0){ (1 + s0)*(s0/(1 + s0))^gamma0 }
Tm <- trainingSize - lag
N <- (1 + monitoringPeriod) * trainingSize
k <- 1:(N - trainingSize)
vaha <- (Tm + k)^2 / (Tm * Qk(k/Tm, gamma)^2)

KT <- max(5, sqrt(log10(trainingSize))) 
sqrt_log_ts <- 2 * sqrt(log10(trainingSize)/trainingSize)

# --- 并行集群设置 ---
num_cores <- parallel::detectCores() - 2
if(num_cores < 1) num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 1. 在每个 Worker 上加载环境 (C++编译)
clusterEvalQ(cl, {
    library(mvtnorm)
    library(Rcpp)
    library(RcppArmadillo)
    library(MTS)
    
    setwd('D:/Desktop/paperRev/JointCUSUM/Mul/') 
    
    # 在 Worker 本地编译 C++
    Rcpp::sourceCpp('statLMul.cpp')
    Rcpp::sourceCpp("eucLMulNew.cpp")
    source('dataGenerationNew.R')
    
    # 辅助函数
    lambda <- function(t) {
        abst <- abs(t)
        (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
    }
})

# 2. 导出变量 (不包含 C++ 函数名)
clusterExport(cl, c("d", "lag", "trainingSize", "monitoringPeriod", "DGP", 
                    "sqrt_log_ts", "KT", "N", "alphaE", "vaha", "kernal", "delta0"))

# 3. 运行并行循环
cat("Running simulation...\n")

# 注意 .noexport 的使用
results <- foreach(s = 1:nsim, .combine = rbind, 
                   .noexport = c("statL", "eucL")) %dopar% {
                       
                       # --- Loop Body ---
                       yt <- dataGeneration(monitoringPeriod, trainingSize, cptModel = DGP, dataDimension = d)
                       
                       Block <- sapply(1:d, function(i) {
                           ut <- yt[1:trainingSize, i]
                           R <- drop(acf(ut, lag.max=2*KT, plot=FALSE, demean=FALSE)$acf)
                           tmp <- which.max(abs(R[1:KT]) >= sqrt_log_ts)
                           M <- if(tmp > 1) 2*(tmp - 1) else 2*(KT - 1)
                           M_vals <- -M:M
                           R_vals <- R[abs(M_vals) + 1]
                           lambda_vals <- lambda(M_vals/M) # 使用 worker 内的 lambda
                           ghat <- sum(lambda_vals * R_vals)
                           Ghat <- sum(lambda_vals * R_vals * abs(M_vals))
                           max(1, round((Ghat^2/ghat^2)^(1/3) * trainingSize^(1/3)))
                       })
                       
                       blocksize <- max(1, round(mean(Block)))
                       
                       num_blocks_est <- ceiling(N / blocksize) * 2
                       lengths <- rgeom(num_blocks_est, 1/blocksize) + 1
                       while (sum(lengths) < N) {
                           needed <- N - sum(lengths)
                           add_blocks <- ceiling(needed / blocksize) + 5
                           lengths <- c(lengths, rgeom(add_blocks, 1/blocksize) + 1)
                       }
                       
                       total_blocks_final <- length(lengths)
                       starts <- sample.int(trainingSize, total_blocks_final, replace=TRUE)
                       idx_list <- lapply(1:total_blocks_final, function(j) {
                           (starts[j] + 0:(lengths[j]-1) - 1) %% trainingSize + 1
                       })
                       
                       y.star <- do.call(rbind, lapply(idx_list, function(idx) {
                           yt[idx, , drop=FALSE]
                       }))[1:N, ]
                       
                       # 调用 Worker 本地的 C++ 函数
                       v1 <- max(statL(yt, trainingSize, m=lag, alphaE, vaha, kernal, delta0))
                       v2 <- max(statL(y.star, trainingSize, m=lag, alphaE, vaha, kernal, delta0))
                       
                       c(v1, v2)
                   }

stopCluster(cl)

# --- 结果分析 ---
Tn <- results[, 1]
Tn.star <- results[, 2]

critvals_final <- quantile(Tn.star, 1-TypeI, na.rm=TRUE)
rej_final <- mean(Tn > critvals_final, na.rm=TRUE)

cat('\nFinal Results:\n')
print(paste("Kernel:", kernal, "| Delta0:", delta0))
print(paste("Final Rejection Rate:", round(100 * rej_final, 2), "%"))