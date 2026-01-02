# 清理环境
rm(list = ls())

# 加载必要的包
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(progress)

# 设置工作目录
setwd('D:/Desktop/JointCUSUM/Mul')

# 编译函数
Rcpp::sourceCpp('statLMul_Old.cpp')
Rcpp::sourceCpp('statLMul.cpp')
Rcpp::sourceCpp('eucLMulNew.cpp')

# 加载数据生成函数
source('dataGenerationSupply.R')

# 设置参数
d <- 2  # 数据维度
DGP <- 1  # 数据生成过程
trainingSize <- 200  # 训练数据长度
monitoringPeriod <- 1  # 监控期
gamma <- 0
deltagauss <- 1 * d
alphaE <- 1  # energy核参数

# 边界函数
Qk <- function(s0, gamma0) {
  (1 + s0) * (s0 / (1 + s0))^gamma0
}

# 生成数据
set.seed(123)
cat("生成测试数据...\n")
yt <- dataGeneration(monitoringPeriod, trainingSize, model = DGP, dataDimension = d)
cat("数据维度:", dim(yt), "\n")

# 设计测试矩阵
test_cases <- expand.grid(
  lag = c(0, 1, 2),  # 测试不同滞后
  kernel = c("gauss", "energy", "quadexp")  # 测试不同核函数
)

# 存储所有结果
all_results <- list()

cat("\n=========== 开始全面对比实验 ===========\n")

# 遍历所有测试情况
for (case_idx in 1:nrow(test_cases)) {
  
  lag_val <- test_cases$lag[case_idx]
  kernal <- as.character(test_cases$kernel[case_idx])
  
  cat(sprintf("\n\n=== 测试情况 %d/%d: lag=%d, kernel=%s ===\n", 
              case_idx, nrow(test_cases), lag_val, kernal))
  
  # 计算相关参数
  N <- nrow(yt)
  nY <- N - lag_val
  Tm <- trainingSize - lag_val
  
  if (nY <= 0 || Tm <= 0) {
    cat("  跳过: nY或Tm小于等于0\n")
    next
  }
  
  k <- 1:(N - trainingSize)
  vaha <- (Tm + k)^2 / (Tm * Qk(k/Tm, gamma)^2)
  vahaY <- vaha[1:(nY - Tm)]
  
  # 计算delta0
  if (kernal == "gauss" || kernal == "energy") {
    delta0 <- deltagauss
  } else if (kernal == "quadexp") {
    # 计算medianE
    E <- eucL(yt[1:trainingSize,], lag_val)
    medianE <- median(E[E > 0])
    delta0 <- medianE
  }
  
  cat(sprintf("  参数: nY=%d, Tm=%d, delta0=%.4f\n", nY, Tm, delta0))
  
  # 1. 正确性对比
  cat("  1. 正确性对比:\n")
  
  result_old <- tryCatch({
    statL_Old(yt, trainingSize, lag_val, alphaE, vahaY, kernal, delta0)
  }, error = function(e) {
    cat(sprintf("    原函数错误: %s\n", e$message))
    return(NULL)
  })
  
  result_new <- tryCatch({
    statL_New(yt, trainingSize, lag_val, alphaE, vahaY, kernal, delta0)
  }, error = function(e) {
    cat(sprintf("    新函数错误: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(result_old) || is.null(result_new)) {
    cat("    跳过性能测试\n")
    next
  }
  
  # 检查维度
  if (length(result_old) != length(result_new)) {
    cat(sprintf("    警告: 结果长度不一致! 原函数: %d, 新函数: %d\n", 
                length(result_old), length(result_new)))
  } else {
    cat(sprintf("    结果长度: %d\n", length(result_old)))
  }
  
  # 计算差异
  diff <- result_old - result_new
  max_abs_diff <- max(abs(diff))
  mean_abs_diff <- mean(abs(diff))
  max_rel_diff <- max(abs(diff) / (abs(result_old) + 1e-10))
  
  cat(sprintf("    最大绝对误差: %.2e\n", max_abs_diff))
  cat(sprintf("    平均绝对误差: %.2e\n", mean_abs_diff))
  cat(sprintf("    最大相对误差: %.2e\n", max_rel_diff))
  
  # 判断结果是否一致
  if (max_abs_diff < 1e-8) {
    cat("    ✓ 结果完全一致\n")
    consistency <- "完全一致"
  } else if (max_rel_diff < 1e-6) {
    cat("    ✓ 结果在可接受范围内\n")
    consistency <- "可接受"
  } else {
    cat("    ⚠ 结果有显著差异\n")
    consistency <- "有差异"
    
    # 显示差异最大的点
    max_idx <- which.max(abs(diff))
    cat(sprintf("    最大差异点[%d]: 原函数=%.6e, 新函数=%.6e\n", 
                max_idx, result_old[max_idx], result_new[max_idx]))
  }
  
  # 2. 性能对比
  cat("  2. 性能对比:\n")
  
  perf <- tryCatch({
    microbenchmark(
      原函数 = statL_Old(yt, trainingSize, lag_val, alphaE, vahaY, kernal, delta0),
      新函数 = statL_New(yt, trainingSize, lag_val, alphaE, vahaY, kernal, delta0),
      times = 10
    )
  }, error = function(e) {
    cat(sprintf("    性能测试错误: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(perf)) {
    perf_summary <- summary(perf)
    old_median <- perf_summary$median[perf_summary$expr == "原函数"]
    new_median <- perf_summary$median[perf_summary$expr == "新函数"]
    speedup <- old_median / new_median
    
    cat(sprintf("    原函数: %.2f ms\n", old_median/1e6))
    cat(sprintf("    新函数: %.2f ms\n", new_median/1e6))
    cat(sprintf("    性能提升: %.2f 倍\n", speedup))
  } else {
    old_median <- NA
    new_median <- NA
    speedup <- NA
  }
  
  # 3. 结果可视化
  if (length(result_old) == length(result_new) && length(result_old) > 0) {
    # 创建对比图
    n_plot <- min(50, length(result_old))
    idx <- 1:n_plot
    
    png(filename = sprintf("对比图_lag%d_%s.png", lag_val, kernal), 
        width = 800, height = 600)
    
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    
    # 子图1: 统计量序列对比
    y_range <- range(c(result_old[idx], result_new[idx]))
    plot(idx, result_old[idx], type = "l", col = "blue", lwd = 2,
         xlab = "索引", ylab = "统计量值", 
         main = "统计量序列对比", ylim = y_range)
    lines(idx, result_new[idx], col = "red", lwd = 2, lty = 2)
    legend("topright", legend = c("原函数", "新函数"), 
           col = c("blue", "red"), lty = c(1, 2), lwd = 2)
    
    # 子图2: 差值序列
    plot(idx, diff[idx], type = "h", col = "darkgreen", lwd = 2,
         xlab = "索引", ylab = "差值 (原 - 新)",
         main = "差值序列")
    abline(h = 0, col = "gray", lty = 2)
    
    # 子图3: 散点图对比
    plot(result_old, result_new, 
         xlab = "原函数结果", ylab = "新函数结果",
         main = "散点图对比", pch = 19, cex = 0.5)
    abline(0, 1, col = "red", lwd = 2)
    
    # 子图4: 差值分布
    hist(diff, breaks = 30, 
         xlab = "差值", ylab = "频率",
         main = "差值分布", col = "lightblue")
    
    # 添加总标题
    title(sprintf("滞后=%d, 核函数=%s, 性能提升=%.2f倍", 
                  lag_val, kernal, ifelse(is.na(speedup), 0, speedup)), 
          outer = TRUE, cex.main = 1.5)
    
    dev.off()
    cat(sprintf("    已保存对比图: 对比图_lag%d_%s.png\n", lag_val, kernal))
  }
  
  # 存储结果
  all_results[[paste0("lag", lag_val, "_", kernal)]] <- list(
    case = test_cases[case_idx, ],
    params = list(nY = nY, Tm = Tm, delta0 = delta0),
    results = list(old = result_old, new = result_new),
    differences = list(
      max_abs_diff = max_abs_diff,
      mean_abs_diff = mean_abs_diff,
      max_rel_diff = max_rel_diff
    ),
    performance = list(
      old_time_ms = ifelse(is.na(old_median), NA, old_median/1e6),
      new_time_ms = ifelse(is.na(new_median), NA, new_median/1e6),
      speedup = speedup
    ),
    consistency = consistency
  )
}

# 生成总结报告
cat("\n\n=========== 实验总结报告 ===========\n")

# 创建总结数据框
summary_df <- data.frame()
for (case_name in names(all_results)) {
  case <- all_results[[case_name]]
  
  summary_df <- rbind(summary_df, data.frame(
    Case = case_name,
    Lag = case$case$lag,
    Kernel = case$case$kernel,
    nY = case$params$nY,
    Tm = case$params$Tm,
    Delta0 = case$params$delta0,
    MaxAbsDiff = case$differences$max_abs_diff,
    MeanAbsDiff = case$differences$mean_abs_diff,
    MaxRelDiff = case$differences$max_rel_diff,
    OldTime_ms = case$performance$old_time_ms,
    NewTime_ms = case$performance$new_time_ms,
    Speedup = case$performance$speedup,
    Consistency = case$consistency,
    stringsAsFactors = FALSE
  ))
}

# 打印总结表格
print(summary_df)

# 写入CSV文件
write.csv(summary_df, "实验结果总结.csv", row.names = FALSE)
cat("\n已保存详细结果到: 实验结果总结.csv\n")

# 分析性能提升
cat("\n=========== 性能提升分析 ===========\n")

# 按核函数分组统计
if (nrow(summary_df) > 0) {
  # 平均性能提升
  avg_speedup <- mean(summary_df$Speedup, na.rm = TRUE)
  cat(sprintf("平均性能提升: %.2f 倍\n", avg_speedup))
  
  # 按核函数分析
  cat("\n按核函数统计:\n")
  kernel_stats <- aggregate(Speedup ~ Kernel, data = summary_df, 
                            FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                min = min(x, na.rm = TRUE), 
                                                max = max(x, na.rm = TRUE)))
  print(kernel_stats)
  
  # 按滞后分析
  cat("\n按滞后统计:\n")
  lag_stats <- aggregate(Speedup ~ Lag, data = summary_df, 
                         FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                             min = min(x, na.rm = TRUE), 
                                             max = max(x, na.rm = TRUE)))
  print(lag_stats)
  
  # 结果一致性分析
  cat("\n结果一致性统计:\n")
  consistency_stats <- table(summary_df$Consistency)
  print(consistency_stats)
}

# 生成可视化报告
cat("\n=========== 生成可视化报告 ===========\n")

if (nrow(summary_df) > 0) {
  png("性能对比总结图.png", width = 1000, height = 800)
  
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
  
  # 图1: 不同核函数的性能提升
  if (length(unique(summary_df$Kernel)) > 1) {
    boxplot(Speedup ~ Kernel, data = summary_df,
            col = c("lightblue", "lightgreen", "lightcoral"),
            ylab = "性能提升倍数", xlab = "核函数",
            main = "不同核函数的性能提升",
            ylim = c(0, max(summary_df$Speedup, na.rm = TRUE) * 1.1))
    abline(h = 1, col = "red", lty = 2, lwd = 2)
  }
  
  # 图2: 不同滞后的性能提升
  if (length(unique(summary_df$Lag)) > 1) {
    boxplot(Speedup ~ Lag, data = summary_df,
            col = "lightyellow",
            ylab = "性能提升倍数", xlab = "滞后参数",
            main = "不同滞后的性能提升",
            ylim = c(0, max(summary_df$Speedup, na.rm = TRUE) * 1.1))
    abline(h = 1, col = "red", lty = 2, lwd = 2)
  }
  
  # 图3: 运行时间对比
  time_data <- rbind(
    data.frame(Time = summary_df$OldTime_ms, Type = "原函数", 
               Case = paste(summary_df$Lag, summary_df$Kernel, sep = "_")),
    data.frame(Time = summary_df$NewTime_ms, Type = "新函数",
               Case = paste(summary_df$Lag, summary_df$Kernel, sep = "_"))
  )
  
  bar_colors <- ifelse(time_data$Type == "原函数", "lightblue", "lightgreen")
  barplot(time_data$Time, names.arg = time_data$Case, 
          col = bar_colors, las = 2, cex.names = 0.7,
          ylab = "运行时间 (ms)", main = "运行时间对比")
  legend("topright", legend = c("原函数", "新函数"), 
         fill = c("lightblue", "lightgreen"))
  
  # 图4: 最大绝对误差分布
  hist(summary_df$MaxAbsDiff, breaks = 20, 
       col = ifelse(summary_df$MaxAbsDiff < 1e-8, "green", 
                    ifelse(summary_df$MaxAbsDiff < 1e-6, "yellow", "red")),
       xlab = "最大绝对误差", ylab = "频率",
       main = "最大绝对误差分布")
  abline(v = 1e-8, col = "blue", lty = 2, lwd = 2)
  abline(v = 1e-6, col = "orange", lty = 2, lwd = 2)
  
  dev.off()
  cat("已保存性能对比总结图: 性能对比总结图.png\n")
}

cat("\n=========== 实验完成 ===========\n")
cat("总结:\n")
cat("1. 共测试了", nrow(test_cases), "种情况\n")
cat("2. 成功测试了", nrow(summary_df), "种情况\n")
cat("3. 平均性能提升:", round(avg_speedup, 2), "倍\n")
cat("4. 结果保存在: 实验结果总结.csv 和 PNG图片中\n")