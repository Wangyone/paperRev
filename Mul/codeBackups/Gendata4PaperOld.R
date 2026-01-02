library(mvtnorm)

#' 生成多元时间序列模拟数据
#'
#' @param monitoringPeriod 监控窗口参数 (Monitoring window factor, L). 
#'   - 这是一个比率因子。监控数据的长度 m = floor(monitoringPeriod * trainSize).
#'   - 例如: 输入 1 表示监控长度与训练长度相等; 输入 0.5 表示监控长度是训练长度的一半.
#' @param trainSize 历史训练样本长度 (Historical/Training sample size, T).
#' @param model 模型编号 (整数)
#' @param cptLocation 变点位置参数 (0, 1) 或 NULL.
#'   - NULL: (默认) 随机变点，范围为监控窗口的前80%。
#'   - (0, 1): 固定变点位置。例如 0.5 表示变点发生在监控窗口的正中间。
#' @param delta 变化幅度参数 (Numeric) 或 NULL.
#'   - NULL: (默认) 使用文档规定的特定模型默认值。
#'     - M1-M3: delta=1
#'     - M4: delta=1
#'     - C1-C3: delta=0.9
#'     - V1-V3: delta=2 (方差变为 3倍)
#'     - V4: delta=1 (方差变为 2倍)
#'   - Numeric: 自定义变化幅度。
#'
#' @return matrix (去除 burn-in 后的数据，行数为 trainSize + m)
gendataPaper <- function(monitoringPeriod, trainSize, model, 
                         cptLocation = NULL, delta = NULL) {
    
    # ==========================================
    # 1. 基础参数与矩阵定义
    # ==========================================
    d <- 2
    burnin <- ceiling(trainSize/2)
    
    # 计算监控数据长度 m = L * T
    len_monitoring <- floor(monitoringPeriod * trainSize)
    N_post <- len_monitoring
    
    # 总样本长度 = Burn-in + Training + Monitoring
    obs <- burnin + trainSize + N_post 
    
    # --- 基础常数 ---
    I2 <- diag(2)
    Zero <- c(0, 0)
    
    # --- 基础参数 (Null Hypothesis) ---
    Sigma_Base <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
    
    # VAR 系数
    Mat_B <- matrix(c(0.1, 0.05, 0.05, 0.1), 2, 2)     # N6 Weak
    Mat_A1 <- matrix(c(0.5, 0.2, 0.2, 0.1), 2, 2)      # N7 Moderate
    
    # 结构变化专用系数 (T组)
    Mat_A2 <- matrix(c(0.8, 0.3, 0.1, 0.7), 2, 2)      # T2
    Mat_A_Flip <- matrix(c(0.5, 0.1, 0.1, 0.5), 2, 2)  # T4/T6
    
    # BEKK-GARCH 参数 (N8)
    C_bekk <- 0.001 * matrix(c(4, 0, 5, 3), 2, 2)
    A_bekk <- matrix(c(0.254, 0.040, -0.004, 0.332), 2, 2)
    B_bekk <- matrix(c(0.941, -0.019, 0.023, 0.864), 2, 2)
    Intercept_bekk <- t(C_bekk) %*% C_bekk
    
    # ==========================================
    # 2. 核心辅助函数
    # ==========================================
    
    # 2.1 变点位置计算
    get_cp <- function() {
        base_idx <- burnin + trainSize # 变点发生在训练期结束之后
        
        if (is.null(cptLocation)) {
            # 默认: 随机分布在监控窗口的前 80%
            return(floor(base_idx + (len_monitoring * runif(1, 0, 0.8))))
        } else {
            # 固定位置
            if(cptLocation <= 0 || cptLocation >= 1) stop("cptLocation must be between 0 and 1.")
            return(floor(base_idx + (len_monitoring * cptLocation)))
        }
    }
    
    # 2.2 静态 VAR 生成器
    sim_var_static <- function(n, phi, sigma, p = 1) {
        k <- nrow(sigma)
        x_out <- matrix(0, n, k)
        noise <- rmvnorm(n, mean = rep(0, k), sigma = sigma)
        x_out[1:p, ] <- noise[1:p, ]
        
        for (t in (p + 1):n) {
            ar_part <- rep(0, k)
            for (j in 1:p) {
                idx_start <- (j - 1) * k + 1
                idx_end   <- j * k
                A_j <- phi[, idx_start:idx_end, drop = FALSE]
                ar_part <- ar_part + A_j %*% x_out[t - j, ]
            }
            x_out[t, ] <- ar_part + noise[t, ]
        }
        return(x_out)
    }
    
    # 2.3 动态 VAR(1) 生成器
    sim_var_dynamic <- function(n, phi_func, sigma_func, intercept_func = function(t) c(0,0)) {
        x_out <- matrix(0, n, d)
        x_out[1, ] <- rmvnorm(1, mean = intercept_func(1), sigma = sigma_func(1))
        
        for (t in 2:n) {
            Phi <- phi_func(t)
            Sigma <- sigma_func(t)
            Intc <- intercept_func(t)
            eps <- as.vector(rmvnorm(1, mean = Zero, sigma = Sigma))
            
            x_out[t, ] <- (Phi %*% x_out[t-1, ]) + eps + Intc
        }
        return(x_out)
    }
    
    # 2.4 动态 VAR(2) 生成器
    sim_var2_dynamic <- function(n, phi1_func, phi2_func, sigma_func) {
        x_out <- matrix(0, n, d)
        x_out[1:2, ] <- rmvnorm(2, mean = Zero, sigma = sigma_func(1))
        
        for (t in 3:n) {
            Phi1 <- phi1_func(t)
            Phi2 <- phi2_func(t)
            Sigma <- sigma_func(t)
            eps <- as.vector(rmvnorm(1, mean = Zero, sigma = Sigma))
            
            x_out[t, ] <- (Phi1 %*% x_out[t-1, ]) + (Phi2 %*% x_out[t-2, ]) + eps
        }
        return(x_out)
    }
    
    # 2.5 BEKK 模拟器
    sim_bekk_process <- function(n, cp = n+1, mu_shift_vec = c(0,0), scale_factor = 1) {
        x_out <- matrix(0, n, 2)
        H <- vector("list", n)
        H[[1]] <- 0.01 * I2
        eps_base <- rmvnorm(n, mean = Zero, sigma = I2)
        x_out[1, ] <- t(chol(H[[1]])) %*% eps_base[1, ]
        
        for (i in 2:n) {
            is_post_change <- (i > cp)
            curr_mu <- if (is_post_change) mu_shift_vec else Zero
            curr_scale <- if (is_post_change) scale_factor else 1
            
            curr_eps <- (eps_base[i, ] * curr_scale) + curr_mu
            
            Ri <- as.matrix(x_out[(i-1), ])
            Inn <- Ri %*% t(Ri)
            H[[i]] <- Intercept_bekk + t(A_bekk) %*% Inn %*% A_bekk + t(B_bekk) %*% H[[(i-1)]] %*% B_bekk
            
            x_out[i, ] <- t(chol(H[[i]])) %*% curr_eps
        }
        return(x_out)
    }
    
    # ==========================================
    # 3. 数据生成逻辑
    # ==========================================
    x <- matrix(0, obs, d)
    
    # ------------------------------------------
    # Group 1: Null Hypothesis
    # ------------------------------------------
    if (model <= 8) {
        if (model == 1) x <- rmvnorm(obs, sigma = Sigma_Base)
        else if (model == 2) x <- rmvt(obs, sigma = (3/5) * Sigma_Base, df = 5)
        else if (model == 3) x <- sim_var_static(obs, phi=0.1*I2, sigma=I2)
        else if (model == 4) x <- sim_var_static(obs, phi=0.3*I2, sigma=I2)
        else if (model == 5) x <- sim_var_static(obs, phi=0.5*I2, sigma=I2)
        else if (model == 6) x <- sim_var_static(obs, phi=Mat_B, sigma=I2)
        else if (model == 7) x <- sim_var_static(obs, phi=Mat_A1, sigma=I2)
        else if (model == 8) x <- sim_bekk_process(obs)
    }
    
    # ------------------------------------------
    # Group 2: Mean Change (M1-M4)
    # ------------------------------------------
    else if (model %in% c(101, 102, 103, 104)) {
        cp <- get_cp()
        d_val <- if(is.null(delta)) 1 else delta
        
        if (model == 101) { # M1
            Mu_Vec <- c(d_val, d_val)
            x_pre <- rmvnorm(cp, mean = Zero, sigma = I2)
            x_post <- rmvnorm(obs-cp, mean = Mu_Vec, sigma = I2)
            x <- rbind(x_pre, x_post)
            
        } else if (model == 102) { # M2
            Mu_Vec <- c(d_val, d_val)
            x <- sim_var_dynamic(obs, 
                                 phi_func = function(t) Mat_B,
                                 sigma_func = function(t) I2,
                                 intercept_func = function(t) if(t > cp) Mu_Vec else Zero)
            
        } else if (model == 103) { # M3
            Mu_Vec <- c(d_val, d_val)
            x <- sim_var_dynamic(obs, 
                                 phi_func = function(t) Mat_A1,
                                 sigma_func = function(t) I2,
                                 intercept_func = function(t) if(t > cp) Mu_Vec else Zero)
            
        } else if (model == 104) { # M4
            Mu_Vec <- c(0, d_val)
            x <- sim_bekk_process(obs, cp = cp, mu_shift_vec = Mu_Vec)
        }
    }
    
    # ------------------------------------------
    # Group 3: Correlation Change (C1-C3)
    # ------------------------------------------
    else if (model %in% c(201, 202, 203)) {
        cp <- get_cp()
        d_val <- if(is.null(delta)) 0.9 else delta
        Sigma_New <- matrix(c(1, d_val, d_val, 1), 2, 2)
        
        if (model == 201) {
            x <- rbind(rmvnorm(cp, sigma = I2), rmvnorm(obs-cp, sigma = Sigma_New))
            
        } else if (model == 202) {
            x <- sim_var_dynamic(obs, 
                                 phi_func = function(t) Mat_B,
                                 sigma_func = function(t) if(t > cp) Sigma_New else I2)
            
        } else if (model == 203) {
            x <- sim_var_dynamic(obs, 
                                 phi_func = function(t) Mat_A1,
                                 sigma_func = function(t) if(t > cp) Sigma_New else I2)
        }
    }
    
    # ------------------------------------------
    # Group 4: Variance Change (V1-V4)
    # ------------------------------------------
    else if (model %in% c(211, 212, 213, 214)) {
        cp <- get_cp()
        
        if (model == 214) {
            # V4: delta默认为1 (Var x 2)
            d_val <- if(is.null(delta)) 1 else delta
            scale_factor <- sqrt(1 + d_val)
            x <- sim_bekk_process(obs, cp = cp, scale_factor = scale_factor)
            
        } else {
            # V1-V3: delta默认为2 (Var x 3)
            d_val <- if(is.null(delta)) 2 else delta
            scale_factor <- sqrt(1 + d_val)
            Sigma_New <- (scale_factor^2) * I2
            
            if (model == 211) {
                x <- rbind(rmvnorm(cp, sigma = I2), rmvnorm(obs-cp, sigma = Sigma_New))
            } else if (model == 212) {
                x <- sim_var_dynamic(obs, 
                                     phi_func = function(t) Mat_B,
                                     sigma_func = function(t) if(t > cp) Sigma_New else I2)
            } else if (model == 213) {
                x <- sim_var_dynamic(obs, 
                                     phi_func = function(t) Mat_A1,
                                     sigma_func = function(t) if(t > cp) Sigma_New else I2)
            }
        }
    }
    
    # ------------------------------------------
    # Group 5: Temporal/Joint Change (T1-T6)
    # ------------------------------------------
    else if (model == 301) { # T1
        cp <- get_cp()
        x <- sim_var_dynamic(obs, 
                             phi_func = function(t) if(t > cp) 0.8*I2 else 0.1*I2,
                             sigma_func = function(t) I2)
        
    } else if (model == 302) { # T2
        cp <- get_cp()
        x <- sim_var_dynamic(obs, 
                             phi_func = function(t) if(t > cp) Mat_A2 else Mat_A1,
                             sigma_func = function(t) I2)
        
    } else if (model == 303) { # T3
        cp <- get_cp()
        get_sigma_t3 <- function(rho) { (1 - rho^2) * I2 }
        x <- sim_var_dynamic(obs,
                             phi_func = function(t) if(t > cp) 0.7*I2 else 0.1*I2,
                             sigma_func = function(t) if(t > cp) get_sigma_t3(0.7) else get_sigma_t3(0.1))
        
    } else if (model == 304) { # T4
        cp <- get_cp()
        x <- sim_var_dynamic(obs,
                             phi_func = function(t) if(t > cp) -Mat_A_Flip else Mat_A_Flip,
                             sigma_func = function(t) I2)
        
    } else if (model == 305) { # T5
        cp <- get_cp()
        get_sigma_t5 <- function(rho) { (1 - rho^2) * I2 }
        x <- sim_var2_dynamic(obs,
                              phi1_func = function(t) 0*I2, 
                              phi2_func = function(t) if(t > cp) 0.7*I2 else 0.1*I2,
                              sigma_func = function(t) if(t > cp) get_sigma_t5(0.7) else get_sigma_t5(0.1))
        
    } else if (model == 306) { # T6
        cp <- get_cp()
        x <- sim_var2_dynamic(obs,
                              phi1_func = function(t) matrix(0, 2, 2),
                              phi2_func = function(t) if(t > cp) -Mat_A_Flip else Mat_A_Flip,
                              sigma_func = function(t) I2)
    } else {
        stop(paste("Error: Model ID", model, "not found."))
    }
    
    # 去除 burn-in
    return(x[-(1:burnin), , drop = FALSE])
}