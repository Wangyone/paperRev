rm(list = ls())
source("X:/Desktop/JointCUSUM/Mul/dataBlock.R", echo = TRUE)

# Parameters
d <- 3
AllData <- 200
phi <- seq(0, 0.9, by = 0.1)  # Autocorrelation coefficients
nsim <- 1000                    # Number of simulations

# Flat-top lag-window function
lambda <- function(t) {
  abst <- abs(t)
  (abst <= 0.5) + 2 * (1 - abst) * ((abst > 0.5) & (abst <= 1))
}

# Optimized block size calculation
blockYt <- function(Xt) {
  d <- NCOL(Xt)
  trainingSize <- NROW(Xt)
  KT <- max(5, sqrt(log10(trainingSize)))  # Optimal lag window
  
  # Calculate block size for each dimension
  blocks <- sapply(1:d, function(i) {
    ut <- Xt[, i]
    R <- as.vector(acf(ut, lag.max = 2 * KT, plot = FALSE)$acf)
    
    # Find M: first lag where |R| < 2*SE
    se <- 2 * sqrt(log10(trainingSize) / trainingSize)
    tmp <- which(abs(R[1:KT]) < se)
    M <- if (length(tmp) > 0) 2 * (tmp[1] - 1) else 2 * (KT - 1)
    # Compute Ghat and D.SB
    lags <- -M:M
    weights <- lambda(lags / M)
    ghat <- sum(weights * R[abs(lags) + 1])
    Ghat <- sum(weights * R[abs(lags) + 1] * abs(lags))
    D.SB <- 2 * ghat^2
    
    # Return block size
    max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainingSize^(1/3)))
  })
  
  round(mean(blocks))  # Final block size
}

## test
# yt <- dataBlock(AllData, 0.9, d)
# blockYt(yt)
# rgeom(1, 1/blockYt(yt)) + 1  

# Initialize result matrix
BlockMatrix <- matrix(NA, nrow = nsim, ncol = length(phi))

# Simulation loop
for (i in 1:nsim) {
  for (j in seq_along(phi)) {
    beta_y <- phi[j]
    yt <- dataBlock(AllData, beta_y, d)
    BlockMatrix[i, j] <- rgeom(1, 1/blockYt(yt)) + 1  
  }
}

# Output results
#print(BlockMatrix)

BlockMatrix <- as.data.frame(BlockMatrix)
colnames(BlockMatrix) <- phi


#### plot 1 箱线图 #### 
boxplot(BlockMatrix[, 1:NROW(phi)], 
        main = " ",
        xlab = TeX("$\\beta$"),
        ylab = TeX("$b^{*}$"),
        col = "skyblue",
        border = "black",
        las = 1)  # las=1使y轴标签水平显示


#### plot 1 误差线图（Error Bar Plot）/ 均值-标准差图（Mean ± SD Plot） #### 
# ErrorBarPlot
# 计算每列的均值和标准差
mean_values <- sapply(BlockMatrix, mean)
sd_values <- sapply(BlockMatrix, sd)

# 设置图形参数
par(mar = c(5, 5, 4, 4) + 0.1)  # 调整边距

# 绘制空图形框架
plot(
  x = phi, y = mean_values, type = "n",
  xlab = expression(paste("  ", beta)),
  ylab = TeX("$b^{*}$"), # b[n]: b_n
  main = " ",
  xlim = c(0, 0.9), ylim = c(0, max(mean_values + sd_values) * 1.1),
  las = 1, cex.lab = 1.2, cex.main = 1.4,
  xaxt = "n"  # 先不绘制默认的X轴
)

# 自定义X轴刻度，按0.1的间隔显示
axis(1, at = seq(0, 0.9, by = 0.1), labels = seq(0, 0.9, by = 0.1))

# 添加误差线（用arrows模拟）
arrows(
  x0 = phi, y0 = mean_values - sd_values,
  x1 = phi, y1 = mean_values + sd_values,
  angle = 90, code = 3, length = 0.05, col = "red"
)

# 添加折线和点
lines(phi, mean_values, col = "blue", lwd = 2)
points(phi, mean_values, pch = 19, col = "blue", cex = 1.5)
# 添加网格线
# grid(nx = NA, ny = NULL, lty = 2, col = "gray")

######### plot 3 ggplot ##################
# 加载必要的包
library(ggplot2)
library(latex2exp)

# 创建绘图数据框
plot_data <- data.frame(
  beta = as.numeric(names(BlockMatrix)),
  mean = mean_values,
  sd = sd_values
)

# 绘制图形
ggplot(plot_data, aes(x = beta, y = mean)) +
  geom_point(size = 3, color = "blue") +          # 绘制均值点
  geom_line(color = "blue", linewidth = 1) +      # 绘制折线
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),     # 误差线表示标准差
    width = 0.02, color = "red"
  ) +
  scale_x_continuous(                             # 设置 X 轴刻度为 0.1 间隔
    breaks = seq(0, 1, by = 0.1),              # 0, 0.1, 0.2, ..., 0.9
    limits = c(0, 0.95)                           # 限制 X 轴范围
  ) +
  scale_y_continuous(                             # 调整 Y 轴范围，避免太满
    limits = c(0, max(plot_data$mean + plot_data$sd) * 1.05)  # 顶部留 5% 空白
  ) +
  labs(
    title = " ",
    x = expression(beta),                         # 更简洁的写法
    y = TeX("$b^{*}$")                           # 使用 LaTeX 数学符号
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)         # 调整图形边距（上、右、下、左）
  )
