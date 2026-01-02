# 移出库加载到全局，避免在函数中重复加载
library(mvtnorm)
library(MTS)
library(parallel)  # 用于并行计算
setwd('D:/Desktop/JointCUSUM/Mul')
rm(list = ls())

dataBlock <- function(sampleSize, beta, dataDimension) {
  ### Data length
  burnin <- ceiling(sampleSize/2)
  Numberdata <- burnin + sampleSize
  
  A <-  diag(dataDimension)
  # 设置上对角线和下对角线为s
  for( i in 1 : (dataDimension - 1) ) {
    A[i, i + 1] <- .5    # 上对角线
    A[i + 1, i] <- 0   # 下对角线
  }
  
  ## 优化协方差矩阵创建
  sigma <- diag(dataDimension)
  # 设置上对角线和下对角线为s
  for(i in 1:(dataDimension - 1)) {
    sigma[i, i+1] <- 0.5    # 上对角线
    sigma[i+1, i] <- 0.5   # 下对角线
  }
  # 生成数据
  x <- ts(VARMAsim(Numberdata, arlags = c(1), 
                   phi = beta * A, #diag(dataDimension),
                   sigma = sigma)$series)
  
  # 移除燃烧期
  yt <- x[-c(1:burnin), ]
  
  return(yt)
}

# Parameters
d <- 3
AllData <- 200
phi <- seq(0, 0.9, by = 0.1)  # Autocorrelation coefficients
nsim <- 2000                  # Number of simulations

# Flat-top lag-window function
lambda <- function(t) {
  abst <- abs(t)
  (abst <= 0.5) + 2 * (1 - abst) * ((abst > 0.5) & (abst <= 1))
}

# 优化blockYt函数，减少不必要的计算和转换
blockYt <- function(Xt) {
  d <- NCOL(Xt)
  trainingSize <- NROW(Xt)
  KT <- max(5, sqrt(log10(trainingSize)))  # Optimal lag window
  
  # 预计算SE值，避免在循环中重复计算
  se <- 2 * sqrt(log10(trainingSize) / trainingSize)
  
  # 计算每个维度的块大小
  blocks <- sapply(1:d, function(i) {
    ut <- Xt[, i]
    # 直接获取ACF值，避免创建acf对象
    R <- acf(ut, lag.max = 2 * KT, plot = FALSE)$acf[, , 1]
    
    # 寻找M值
    tmp <- which(abs(R[1:KT]) < se)
    M <- if (length(tmp) > 0) 2 * (tmp[1] - 1) else 2 * (KT - 1)
    
    # 计算滞后和权重
    lags <- -M:M
    weights <- lambda(lags / M)
    ghat <- sum(weights * R[abs(lags) + 1])
    Ghat <- sum(weights * R[abs(lags) + 1] * abs(lags))
    D.SB <- 2 * ghat^2
    
    # 返回块大小
    max(1, round((2 * Ghat^2 / D.SB)^(1/3) * trainingSize^(1/3)))
  })
  
  round(mean(blocks))  # 最终块大小
}

# 使用并行计算优化模拟循环
# 获取可用核心数
num_cores <- detectCores() - 1  # 保留一个核心给系统
cl <- makeCluster(num_cores)

# 导出必要的变量和函数到集群
clusterExport(cl, c("dataBlock", "blockYt", "lambda", 
                    "AllData", "d", "phi"))
clusterEvalQ(cl, {
  library(MTS)  # 在每个集群节点加载必要的库
})

# 并行执行模拟
BlockMatrix <- parLapply(cl, 1:nsim, function(i) {
  # 对每个phi值计算结果
  sapply(phi, function(beta_y) {
    yt <- dataBlock(AllData, beta_y, d)
    rgeom(1, 1/blockYt(yt)) + 1  
  })
})

# 关闭集群
stopCluster(cl)

# 转换结果为数据框
BlockMatrix <- do.call(rbind, BlockMatrix)
BlockMatrix <- as.data.frame(BlockMatrix)
colnames(BlockMatrix) <- phi

#### 箱线图 #### 
library(latex2exp)
boxplot(BlockMatrix[, 1:length(phi)], 
        main = " ",
        xlab = expression(beta),
        ylab = expression(b^"*"),
        col = "skyblue",
        border = "black",
        las = 1)  # las=1使y轴标签水平显示


######### 箱线图1### 


library(ggplot2)
library(dplyr)    # 加载dplyr包，提供everything()函数
library(tidyr)    # 加载tidyr包，提供pivot_longer()函数

# 将数据转换为长格式以便使用ggplot2
BlockMatrix_long <- tidyr::pivot_longer(BlockMatrix, 
                                        cols = everything(),
                                        names_to = "beta", 
                                        values_to = "b_star")
BlockMatrix_long$beta <- as.numeric(BlockMatrix_long$beta)

# 使用ggplot2创建更美观的箱线图
ggplot(BlockMatrix_long, aes(x = factor(beta), y = b_star, fill = beta)) +
  geom_boxplot(alpha = 0.7,  # 半透明效果
               width = 0.6,  # 箱体宽度
               outlier.size = 1.5,  # 异常值大小
               outlier.color = "#FF6B6B") +  # 异常值颜色
  scale_fill_viridis_c(option = "plasma") +  # 渐变色填充
  labs(title = " ",
       x = TeX("$\\beta$"),
       y = TeX("$b^{*}$")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
    legend.position = "none"  # 移除不必要的图例
  ) +
  ylim(0, max(BlockMatrix_long$b_star) * 1.1)  # 适当调整y轴范围



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

# 计算每列的均值和标准差
mean_values <- sapply(BlockMatrix, mean)
sd_values <- sapply(BlockMatrix, sd)

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









######## 4 优化可视化 ########



# （前提：BlockMatrix已通过文档相关模拟生成，此处延续原计算逻辑）
# 计算每列均值和标准差（保留3位小数，提升数据精度呈现）
mean_values <- sapply(BlockMatrix, function(x) round(mean(x), 3))
sd_values <- sapply(BlockMatrix, function(x) round(sd(x), 3))

# 创建绘图数据框（新增变异系数列，可选用于后续标注）
plot_data <- data.frame(
  beta = as.numeric(names(BlockMatrix)),
  mean = mean_values,
  sd = sd_values,
  cv = round(sd_values / mean_values, 3)  # 变异系数，反映波动相对程度
)

# 绘制优化后的图形
ggplot(plot_data, aes(x = beta, y = mean)) +
  # 1. 误差线优化：加粗线条+半透明，避免与均值点重叠
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.015,          # 误差线横向宽度适配beta间隔（0.1），更纤细美观
    color = "#E74C3C",      # 暗红色误差线，比原红色更柔和且突出
    linewidth = 0.8,        # 加粗误差线，提升可读性
    alpha = 0.8             # 半透明效果，避免遮挡后续元素
  ) +
  # 2. 均值点优化：填充色+边框，增强立体感
  geom_point(
    size = 4,               # 点尺寸放大，适配学术图表阅读距离
    color = "#2980B9",      # 深蓝色边框
    fill = "#3498DB",       # 浅蓝色填充
    shape = 21,             # 带边框的圆形，区分于普通实心点
    stroke = 1.2            # 点边框加粗，增强视觉冲击
  ) +
  # 3. 均值折线优化：平滑线条+渐变效果适配文档中"序列相依性增强"逻辑
  geom_line(
    color = "#2980B9",
    linewidth = 1.2,
    linetype = "solid",     # 实线更符合学术图表严谨性
    alpha = 0.9
  ) +
  # 4. 坐标轴优化：严格匹配文档中beta范围（0-0.9），刻度间隔精准
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.1),  # 仅保留文档中使用的beta值（0,0.1,...,0.9）
    limits = c(-0.02, 0.92),         # 左右留微小空白，避免点/线贴边
    expand = c(0, 0)                 # 取消默认扩展，确保刻度与数据范围一致
  ) +
  scale_y_continuous(
    # 动态调整Y轴范围，顶部留10%空白（比原5%更适配误差线显示）
    limits = c(0, max(plot_data$mean + plot_data$sd) * 1.1),
    expand = c(0, 0),
    # 新增Y轴次要刻度，辅助读取具体数值（尤其适配文档中块大小的量化分析）
    minor_breaks = seq(0, max(plot_data$mean + plot_data$sd) * 1.1, by = 0.5)
  ) +
  # 5. 标签优化：适配文档学术表述，补充单位/说明
  labs(
    x = TeX("$\\beta$ "),  # 补充beta的物理意义（自相关系数）
    y = TeX("$b^{*}$ "),          # 明确b*为"估计块大小"，呼应文档
    caption = TeX("")  # 补充注释，说明误差线含义（文档中nsim=1000）
  ) +
  # 6. 主题优化：适配学术期刊排版，减少冗余元素
  theme_minimal() +
  theme(
    # 坐标轴标签：加粗+合适字号，适配学术图表要求
    axis.title.x = element_text(
      size = 13, 
      face = "bold", 
      margin = margin(t = 8, unit = "pt")  # 顶部留间距，避免与刻度重叠
    ),
    axis.title.y = element_text(
      size = 13, 
      face = "bold", 
      margin = margin(r = 8, unit = "pt")  # 右侧留间距
    ),
    # 坐标轴刻度：清晰可读，避免拥挤
    axis.text.x = element_text(size = 11, color = "#34495E"),
    axis.text.y = element_text(size = 11, color = "#34495E"),
    # 网格线：弱化次要网格，突出数据趋势（符合文档中"相依性越强，块大小越大"的结论）
    panel.grid.major = element_line(color = "#ECF0F1", linewidth = 0.8),
    panel.grid.minor = element_line(color = "#F8F9FA", linewidth = 0.5),
    # 图注：小字号+左对齐，补充关键信息不喧宾夺主
    plot.caption = element_text(
      size = 10, 
      color = "#7F8C8D", 
      hjust = 0, 
      margin = margin(t = 10, unit = "pt")
    ),
    # 边距：适配后续插入文档/报告，避免裁剪
    plot.margin = margin(15, 15, 15, 15, unit = "pt")
  ) 
# （可选）添加均值数值标注，直接呈现关键结果（文档中需突出块大小均值随beta的变化）
# +geom_text(
#     aes(label = mean),          # 标注均值数值
#     vjust = -1.2,               # 位于均值点上方，避免重叠
#     size = 3.5,                 # 字号适配点大小
#     color = "#2C3E50",
#     fontface = "italic"         # 斜体区分于其他文本
#   )





