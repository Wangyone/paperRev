#' 生成多变量时间序列数据
#' 
#' @param monitoringPeriod 监控期长度参数 (L)
#' @param trainingSize 训练集长度 (T)
#' @param cptModel 模型编号
#' @param dataDimension 数据维度
#' @param cptLocation 变点位置 (0-1), NULL则随机
#' @export
dataGeneration <- function(monitoringPeriod, trainingSize, cptModel, dataDimension, cptLocation = NULL) {
  
  # --- 1. 参数检查 ---
  if (!is.numeric(monitoringPeriod) || monitoringPeriod <= 0) stop("参数monitoringPeriod必须为正数")
  if (!is.numeric(trainingSize) || trainingSize <= 0) stop("参数trainingSize必须为正数")
  if (!is.numeric(cptModel) || cptModel < 1) stop("参数cptModel必须为大于0的整数")
  if (!is.numeric(dataDimension) || dataDimension < 1) stop("参数dataDimension必须为正整数")
  if (!is.null(cptLocation) && (cptLocation <= 0 || cptLocation >= 1)) {
    warning("cptLocation参数应在(0,1)范围内，将使用随机生成")
    cptLocation <- NULL
  }
  
  # --- 2. 加载包 ---
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("请安装mvtnorm包: install.packages('mvtnorm')")
  }
  library(mvtnorm)
  
  # --- 3. 基础设置 ---
  # T = trainingSize, L = monitoringPeriod
  burnin <- ceiling(trainingSize/2)
  N <- (1 + monitoringPeriod) * trainingSize
  obs <- burnin + N
  
  rho0 <- 0
  mu0 <- rep(0, dataDimension)
  mu1 <- c(rep(1, ceiling(dataDimension/5)), rep(0, (dataDimension - ceiling(dataDimension/5))))
  
  # --- 4. 辅助函数 ---
  createCormat <- function(d, s) {
    if (d == 1) return(matrix(1))
    mat <- matrix(NA, nrow = d, ncol = d)
    matrix(ifelse(row(mat) == col(mat), 1, s), nrow = d)
  }
  
  createPHI <- function(diagonal, off_diagonal) {
    mat <- diag(diagonal, nrow = dataDimension, ncol = dataDimension)
    if (dataDimension > 1) {
      for (i in 1:(dataDimension - 1)) {
        mat[i, i + 1] <- off_diagonal
        mat[i + 1, i] <- off_diagonal
      }
    }
    return(mat)
  }
  
  reps0 <- function(n, kr) {
    rmvnorm(n, mean = rep(0, dataDimension), sigma = createCormat(dataDimension, kr))
  }
  
  # 模拟VAR过程 (自带内部预热，返回长度为n的数据)
  VarSim <- function(n, phi, sigma) {
    if (n <= 0) return(matrix(0, nrow=0, ncol=ncol(phi))) 
    burnin_sim <- ceiling(n/2)
    total_n <- n + burnin_sim
    d <- ncol(phi)
    x <- matrix(0, nrow = total_n, ncol = d)
    
    eps <- rmvnorm(total_n, mean = rep(0, d), sigma = sigma)
    x[1, ] <- rmvnorm(1, mean = rep(0, d), sigma = sigma)
    
    if (total_n > 1) {
      for (t in 2:total_n) {
        x[t, ] <- phi %*% x[t-1, ] + eps[t, ]
      }
    }
    return(x[(burnin_sim + 1):total_n, , drop = FALSE])
  }
  
  # --- 5. 初始化参数 ---
  sigma0 <- createCormat(dataDimension, rho0)
  A4 <- createPHI(diagonal = 0.1, off_diagonal = 0.04)
  A5 <- createPHI(diagonal = 0.5, off_diagonal = 0.1)
  
  # 获取变点位置 (在 obs 尺度上的索引)
  getChangePoint <- function() {
    if (!is.null(cptLocation)) {
      return(floor(burnin + trainingSize * (1 + monitoringPeriod * cptLocation)))
    } else {
      # 对应 old code: floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
      return(floor(burnin + trainingSize * (1 + monitoringPeriod * runif(1, 0, 4/5))))
    }
  }
  
  # --- 6. 数据生成逻辑 ---
  # 先生成长度为 obs 的完整数据 x_full，最后再截取
  x_full <- matrix(0, obs, dataDimension)
  epsN <- reps0(obs, rho0)
  
  # === 无结构变化 (1-5) ===
  if (cptModel == 1) {
    x_full <- epsN
  } else if (cptModel == 2) {
    x_full <- VarSim(obs, phi = diag(0.3, dataDimension), sigma = sigma0)
  } else if (cptModel == 3) {
    x_full <- VarSim(obs, phi = diag(0.5, dataDimension), sigma = sigma0)
  } else if (cptModel == 4) {
    x_full <- VarSim(obs, phi = A4, sigma = sigma0)
  } else if (cptModel == 5) {
    x_full <- VarSim(obs, phi = A5, sigma = sigma0)
  }
  
  # === 均值变化 (101-103) ===
  else if (cptModel == 101) {
    delta <- 1
    x_full <- epsN
    change.pt <- getChangePoint()
    if(change.pt < obs) x_full[change.pt:obs, ] <- x_full[change.pt:obs, ] + delta * mu1
  } else if (cptModel == 102) {
    delta <- 1
    x_full <- VarSim(obs, phi = A4, sigma = sigma0)
    change.pt <- getChangePoint()
    if(change.pt < obs) x_full[change.pt:obs, ] <- x_full[change.pt:obs, ] + delta * mu1
  } else if (cptModel == 103) {
    delta <- 1
    x_full <- VarSim(obs, phi = A5, sigma = sigma0)
    change.pt <- getChangePoint()
    if(change.pt < obs) x_full[change.pt:obs, ] <- x_full[change.pt:obs, ] + delta * mu1
  }
  
  # === 方差变化 (201-203) ===
  else if (cptModel == 201) {
    delta <- 1
    x_full <- rmvnorm(obs, mean = mu0, sigma = diag(dataDimension))
    y <- rmvnorm(obs, mean = mu0, sigma = (1 + delta) * diag(dataDimension))
    change.pt <- getChangePoint()
    if(change.pt <= obs) x_full[change.pt:obs, ] <- y[change.pt:obs, ]
  } 
  else if (cptModel == 202 || cptModel == 203) {
    delta <- 1
    phi_curr <- if(cptModel == 202) A4 else A5
    
    x_full <- VarSim(obs, phi = phi_curr, sigma = diag(dataDimension))
    y_full <- VarSim(obs, phi = phi_curr, sigma = (1 + delta) * diag(dataDimension))
    
    change.pt <- getChangePoint()
    # 替换逻辑：从变点位置开始替换
    if (change.pt <= obs) {
      x_full[change.pt:obs, ] <- y_full[change.pt:obs, ]
    }
  }
  
  # === 时间依赖性变化 (301) ===
  else if (cptModel == 301) {
    delta <- 0.8
    change.pt <- getChangePoint()
    
    # 逻辑参考 old code: rbind(xt[1:change.pt,], yt[(change.pt+1):obs,])
    # 第一段长度: change.pt
    # 第二段长度: obs - change.pt
    len1 <- change.pt
    len2 <- obs - change.pt
    
    xt_part <- matrix(0, nrow=0, ncol=dataDimension)
    if(len1 > 0) xt_part <- VarSim(len1, phi = diag(0.1, dataDimension), sigma = sigma0)
    
    yt_part <- matrix(0, nrow=0, ncol=dataDimension)
    if(len2 > 0) yt_part <- VarSim(len2, phi = diag(delta, dataDimension), sigma = sigma0)
    
    x_full <- rbind(xt_part, yt_part)
  }
  
  # === 联合分布变化 (302) ===
  else if (cptModel == 302) {
    change.pt <- getChangePoint()
    
    len1 <- change.pt
    len2 <- obs - change.pt
    
    # --- Part 1: AR(1) rho=0.1 ---
    xt_part <- matrix(0, len1, dataDimension)
    if (len1 > 0) {
      xt_part[1, ] <- rmvnorm(1, mean = mu0, sigma = diag(dataDimension))
      if (len1 > 1) {
        rhoXt <- 0.1
        for (i in 2:len1) {
          xt_part[i, ] <- rhoXt * xt_part[i - 1, ] + 
            rmvnorm(1, mean = mu0, sigma = (1 - rhoXt^2) * diag(dataDimension))
        }
      }
    }
    
    # --- Part 2: AR(1) rho=0.7 ---
    yt_part <- matrix(0, len2, dataDimension)
    if (len2 > 0) {
      yt_part[1, ] <- rmvnorm(1, mean = mu0, sigma = diag(dataDimension))
      if (len2 > 1) {
        rhoYt <- 0.7
        for (i in 2:len2) {
          yt_part[i, ] <- rhoYt * yt_part[i - 1, ] + 
            rmvnorm(1, mean = mu0, sigma = (1 - rhoYt^2) * diag(dataDimension))
        }
      }
    }
    
    x_full <- rbind(xt_part, yt_part)
  }
  
  # === 联合分布变化 (303) ===
  else if (cptModel == 303) {
    change.pt <- getChangePoint()
    
    len1 <- change.pt
    len2 <- obs - change.pt
    
    # --- Part 1: White Noise ---
    xt_part <- matrix(0, len1, dataDimension)
    if (len1 > 0) {
      xt_part <- rmvnorm(len1, mean = mu0, sigma = diag(dataDimension))
    }
    
    # --- Part 2: AR(2) rho=0.7 ---
    yt_part <- matrix(0, len2, dataDimension)
    if (len2 > 0) {
      yt_part[1, ] <- rmvnorm(1, mean = mu0, sigma = diag(dataDimension))
      if (len2 > 1) {
        yt_part[2, ] <- rmvnorm(1, mean = mu0, sigma = diag(dataDimension))
        if (len2 > 2) {
          rhoYt <- 0.7
          for (i in 3:len2) {
            yt_part[i, ] <- rhoYt * yt_part[i - 2, ] + 
              rmvnorm(1, mean = mu0, sigma = (1 - rhoYt^2) * diag(dataDimension))
          }
        }
      }
    }
    
    x_full <- rbind(xt_part, yt_part)
  } 
  else {
    stop("无效的模型编号。请使用1-5, 101-103, 201-203, 301-303范围内的模型。")
  }
  
  # --- 7. 去除预热期并返回 ---
  # 统一在最后一步进行切片，确保逻辑最简单
  return(x_full[-(1:burnin), , drop = FALSE])
}