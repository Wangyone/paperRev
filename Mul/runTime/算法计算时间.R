rm(list = ls())
######################### Performance Testing #########################
setwd('D:/Desktop/JointCUSUM/Mul')
library(ggplot2)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)

sourceCpp('statLMul.cpp')
source('dataGenerationSupply.R')

## boundary function 
Qk <- function(s0, gamma0) {
  u = (1 + s0) * (s0 / (1 + s0))^gamma0
  return(u)
}

############################ Robust Performance Test Function
test_statL_performance <- function(dimensions = c(2, 5, 10, 20, 30), 
                                   n_runs = 500, 
                                   training_size = 100,
                                   monitoring_period = 1,
                                   lag = 0,
                                   kernel_type = 1) {
  
  # Set kernel parameters
  if (kernel_type == 1) {
    kernal <- "gauss"
    delta0 <- 1  # fixed delta for gauss
  } else if (kernel_type == 2) {
    kernal <- "energy"
    delta0 <- 1  # fixed delta for energy
  } else if (kernel_type == 3) {
    kernal <- "quadexp"
    delta0 <- 1  # fixed delta for quadexp
  }
  
  # Prepare weight function
  Tm <- training_size - lag
  N <- (1 + monitoring_period) * training_size
  k <- 1:(N - training_size)
  gamma <- 0
  vaha <- (Tm + k)^2 / (Tm * Qk(k / Tm, gamma)^2)
  
  results <- data.frame()
  
  for (d in dimensions) {
    cat("Testing dimension:", d, "\n")
    
    # Manual timing without microbenchmark
    times <- numeric(n_runs)
    
    for (i in 1:n_runs) {
      # Generate test data
      yt <- dataGeneration(monitoring_period, training_size, model = 1, dataDimension = d)
      start_time <- Sys.time()
      test_result <- statL(yt, training_size, m = lag, alphaV = 1, 
                           vahaY = vaha, kern = kernal, delta = delta0)
      end_time <- Sys.time()
      times[i] <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000  # convert to ms
    }
    
    # Calculate statistics manually
    time_stats <- list(
      mean = mean(times),
      median = median(times),
      min = min(times),
      max = max(times),
      sd = sd(times)
    )
    
    # Store results
    dim_result <- data.frame(
      dimension = d,
      mean_time_ms = time_stats$mean,
      median_time_ms = time_stats$median,
      min_time_ms = time_stats$min,
      max_time_ms = time_stats$max,
      sd_time_ms = time_stats$sd,
      n_runs = n_runs
    )
    
    results <- rbind(results, dim_result)
    
    # Print progress
    cat(sprintf("  Mean time: %.2f ms, Median: %.2f ms\n", 
                dim_result$mean_time_ms, dim_result$median_time_ms))
  }
  
  return(results)
}

############################ Quick Test Function (for larger dimensions)
quick_performance_test <- function(dimensions = c(2, 10, 30, 50, 100), 
                                   n_runs = 30,
                                   training_size = 50) {  
  # smaller training size for high dimensions
  
  cat("=== Quick Performance Test ===\n")
  results <- data.frame()
  
  for (d in dimensions) {
    cat("Testing dimension:", d, "\n")
    
    # Generate smaller dataset for high dimensions to avoid memory issues
    current_training_size <- ifelse(d > 30, 30, training_size)
    yt <- dataGeneration(1, current_training_size, model = 1, dataDimension = d)
    
    Tm <- current_training_size - 0
    N <- (1 + 1) * current_training_size
    k <- 1:(N - current_training_size)
    gamma <- 0
    vaha <- (Tm + k)^2 / (Tm * Qk(k / Tm, gamma)^2)
    
    # Time a single run first to check feasibility
    start_time <- Sys.time()
    test_result <- statL(yt, current_training_size, m = 0, alphaV = 1, 
                         vahaY = vaha, kern = "gauss", delta = 1)
    single_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs")) * 1000
    
    cat(sprintf("  Single run time: %.2f ms\n", single_time))
    
    # If single run takes too long, reduce number of runs
    current_runs <- ifelse(single_time > 1000, 10, 
                           ifelse(single_time > 500, 20, n_runs))
    
    if (current_runs < n_runs) {
      cat(sprintf("  Reduced to %d runs due to long computation time\n", current_runs))
    }
    
    # Manual timing
    if (current_runs > 1) {
      bench_times <- numeric(current_runs)
      for (i in 1:current_runs) {
        start_time <- Sys.time()
        statL(yt, current_training_size, m = 0, alphaV = 1, 
              vahaY = vaha, kern = "gauss", delta = 1)
        bench_times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs")) * 1000
      }
      
      dim_result <- data.frame(
        dimension = d,
        training_size = current_training_size,
        mean_time_ms = mean(bench_times),
        median_time_ms = median(bench_times),
        min_time_ms = min(bench_times),
        max_time_ms = max(bench_times),
        sd_time_ms = sd(bench_times),
        n_runs = current_runs
      )
    } else {
      dim_result <- data.frame(
        dimension = d,
        training_size = current_training_size,
        mean_time_ms = single_time,
        median_time_ms = single_time,
        min_time_ms = single_time,
        max_time_ms = single_time,
        sd_time_ms = 0,
        n_runs = 1
      )
    }
    
    results <- rbind(results, dim_result)
  }
  
  return(results)
}

############################ Simple Test Function (最小化依赖)
simple_performance_test <- function(dimensions = c(2, 5, 10), 
                                    n_runs = 20,
                                    training_size = 100) {
  
  cat("=== Simple Performance Test ===\n")
  results <- data.frame()
  
  for (d in dimensions) {
    cat("Testing dimension:", d, "\n")
    
    # Generate test data
    yt <- dataGeneration(1, training_size, model = 1, dataDimension = d)
    
    # Prepare parameters
    Tm <- training_size - 0
    N <- 2 * training_size
    k <- 1:(N - training_size)
    gamma <- 0
    vaha <- (Tm + k)^2 / (Tm * Qk(k / Tm, gamma)^2)
    
    # Manual timing
    times <- numeric(n_runs)
    
    for (i in 1:n_runs) {
      start_time <- Sys.time()
      test_result <- statL(yt, training_size, m = 0, alphaV = 1, 
                           vahaY = vaha, kern = "gauss", delta = 1)
      end_time <- Sys.time()
      times[i] <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
    }
    
    # Create result row
    dim_result <- data.frame(
      dimension = d,
      mean_time_ms = mean(times),
      median_time_ms = median(times),
      min_time_ms = min(times),
      max_time_ms = max(times),
      sd_time_ms = sd(times),
      n_runs = n_runs
    )
    
    results <- rbind(results, dim_result)
    cat(sprintf("  Completed: %.2f ms average\n", mean(times)))
  }
  
  return(results)
}

############################ Visualization Function
plot_performance_results <- function(results) {
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("ggplot2 not available, creating basic plot\n")
    plot(results$dimension, results$mean_time_ms, type = "b", 
         xlab = "维数", ylab = "平均计算时间(ms)",
         main = "statL Performance by Dimension")
    return()
  }
  
  p <- ggplot2::ggplot(results, ggplot2::aes(x = dimension, y = mean_time_ms)) +
    ggplot2::geom_line(size = 1, color = "steelblue") +
    ggplot2::geom_point(size = 2, color = "steelblue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_time_ms - sd_time_ms, 
                                      ymax = mean_time_ms + sd_time_ms), 
                         alpha = 0.2, fill = "steelblue") +
    ggplot2::labs(
      title = " ",
      x = "维数(d)",
      y = "平均计算时间(ms)",
      subtitle = " "
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = unique(results$dimension)) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  
  print(p)
  return(p)
}

# 可选：创建更专业的学术图表风格
plot_performance_academic <- function(results) {
  # 使用基础绘图系统创建学术风格的图表
  par(mar = c(5, 6, 4, 2) + 0.1, family = "serif")  # 调整边距，使用衬线字体
  
  # 计算合适的y轴范围
  y_max <- max(results$mean_time_ms + results$sd_time_ms, na.rm = TRUE) * 1.15
  
  # 绘制空图形框架
  plot(
    x = results$dimension, y = results$mean_time_ms, type = "n",
    xlab = "维数 (d)",
    ylab = "",
    main = "",
    xlim = c(min(results$dimension) - 0.5, max(results$dimension) + 0.5),
    ylim = c(0, y_max),
    las = 1, cex.lab = 1.4, cex.axis = 1.2,
    xaxt = "n", bty = "l"  # 只有左和下边框
  )
  
  # 添加Y轴标签
  mtext("平均计算时间 (ms)", side = 2, line = 4, cex = 1.4)
  
  # 自定义X轴刻度
  axis(1, at = results$dimension, labels = results$dimension, cex.axis = 1.2)
  
  # 添加浅色网格线
  grid(nx = NA, ny = NULL, lty = 3, col = "gray80")
  abline(h = pretty(c(0, y_max)), lty = 3, col = "gray80")
  
  # 添加误差线
  # arrows(
  #   x0 = results$dimension, y0 = results$mean_time_ms - results$sd_time_ms,
  #   x1 = results$dimension, y1 = results$mean_time_ms + results$sd_time_ms,
  #   angle = 90, code = 3, length = 0.08, col = "#E74C3C", lwd = 2
  # )
  
  # 添加折线
  lines(results$dimension, results$mean_time_ms, col = "#2980B9", lwd = 3)
  
  # 添加数据点
  points(results$dimension, results$mean_time_ms, pch = 21, 
         bg = "white", col = "#2980B9", cex = 1.5, lwd = 1)
  
  # # 可选：添加图例
  # legend("topleft", 
  #        legend = c("均值", "±标准差"),
  #        col = c("#2980B9", "#E74C3C"),
  #        lty = c(1, 1), lwd = c(3, 2),
  #        pch = c(21, NA),
  #        pt.bg = "white",
  #        bty = "n", cex = 1.1)
}


############################ Main Execution
# 首先运行简单测试确保基础功能正常
cat("=== Running Simple Test First ===\n")
simple_results <- simple_performance_test(
  dimensions = c(2, 3, 5),  # 从小的维度开始
  n_runs = 10,              # 较少的运行次数
  training_size = 50        # 较小的训练集
)

print(simple_results)

# 如果简单测试成功，运行更全面的测试

  cat("\n=== Running Standard Performance Test ===\n")
  standard_results <- test_statL_performance(
    dimensions = c(2, 3, 5, 8, 10, 15, 20, 30, 50),
    n_runs = 500,
    training_size = 300,
    kernel_type = 1  # gauss kernel
  )
  
  print(standard_results)
  
  # 可视化
  plot_performance_results(standard_results)
  plot_performance_academic(standard_results) 
  
  # 保存结果
  write.csv(standard_results, "300_performance_results.csv", row.names = FALSE)
  cat("Results saved to 100_performance_results.csv\n")
  

# 打印系统信息
cat("\n=== System Information ===\n")
cat("R version:", R.version$version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Number of cores:", parallel::detectCores(), "\n")