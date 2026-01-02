library(mvtnorm)

#' 生成多元时间序列模拟数据 (修复一元模型生成逻辑版 - Old架构)
#'
#' @param monitoringPeriod 监控窗口比率 (L)
#' @param trainSize 训练样本长度 (T)
#' @param model 模型编号
#' @param cptLocation 变点位置
#' @param delta 变化幅度
#'
#' @return matrix
gendataPaper <- function(monitoringPeriod, trainSize, model, 
                         cptLocation = NULL, delta = NULL) {
    
    # ==========================================
    # 1. 基础参数与维度判定
    # ==========================================
    
    # 定义一元模型列表
    univariate_models <- c(3, 4, 5, 101, 301, 303, 305)
    
    # 动态设定维度
    if (model %in% univariate_models) {
        d <- 1
    } else {
        d <- 2
    }
    
    burnin <- ceiling(trainSize/2)
    len_monitoring <- floor(monitoringPeriod * trainSize)
    N_post <- len_monitoring
    obs <- burnin + trainSize + N_post 
    
    # --- 基础常数 (动态化) ---
    I_d <- diag(d)        # d=1时为1x1矩阵，d=2时为2x2
    Zero <- rep(0, d)     # d=1时为0，d=2时为c(0,0)
    
    # --- 基础参数 (仅用于二元模型的参数保持不变) ---
    Sigma_Base <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
    
    Mat_B  <- matrix(c(0.1, 0.05, 0.05, 0.1), 2, 2)
    Mat_A1 <- matrix(c(0.5, 0.2, 0.2, 0.1), 2, 2)
    Mat_A2 <- matrix(c(0.8, 0.3, 0.1, 0.7), 2, 2)
    Mat_A_Flip <- matrix(c(0.5, 0.1, 0.1, 0.5), 2, 2)
    
    # BEKK 参数 (仅二元)
    C_bekk <- 0.001 * matrix(c(4, 0, 5, 3), 2, 2)
    A_bekk <- matrix(c(0.254, 0.040, -0.004, 0.332), 2, 2)
    B_bekk <- matrix(c(0.941, -0.019, 0.023, 0.864), 2, 2)
    Intercept_bekk <- t(C_bekk) %*% C_bekk
    
    # ==========================================
    # 2. 核心辅助函数
    # ==========================================
    
    get_cp <- function() {
        base_idx <- burnin + trainSize
        if (is.null(cptLocation)) {
            return(floor(base_idx + (len_monitoring * runif(1, 0, 0.8))))
        } else {
            if(cptLocation <= 0 || cptLocation >= 1) stop("cptLocation error")
            return(floor(base_idx + (len_monitoring * cptLocation)))
        }
    }
    
    # 静态 VAR 生成器 (兼容 d=1)
    sim_var_static <- function(n, phi, sigma, p = 1) {
        k <- nrow(sigma) # 自动检测维度
        x_out <- matrix(0, n, k)
        noise <- rmvnorm(n, mean = rep(0, k), sigma = sigma)
        x_out[1:p, , drop=FALSE] <- noise[1:p, , drop=FALSE]
        
        for (t in (p + 1):n) {
            ar_part <- rep(0, k)
            for (j in 1:p) {
                idx_start <- (j - 1) * k + 1
                idx_end   <- j * k
                A_j <- phi[, idx_start:idx_end, drop = FALSE]
                x_lag <- x_out[t - j, ]
                ar_part <- ar_part + as.vector(A_j %*% x_lag)
            }
            x_out[t, ] <- ar_part + noise[t, ]
        }
        return(x_out)
    }
    
    # 动态 VAR(1) 生成器 (用于参数随时间变化的场景)
    sim_var_dynamic <- function(n, phi_func, sigma_func, intercept_func = function(t) Zero) {
        x_out <- matrix(0, n, d)
        x_out[1, ] <- rmvnorm(1, mean = intercept_func(1), sigma = sigma_func(1))
        
        for (t in 2:n) {
            Phi <- phi_func(t)
            Sigma <- sigma_func(t)
            Intc <- intercept_func(t)
            eps <- as.vector(rmvnorm(1, mean = Zero, sigma = Sigma))
            
            val <- (Phi %*% x_out[t-1, ]) + eps + Intc
            x_out[t, ] <- as.vector(val)
        }
        return(x_out)
    }
    
    # 动态 VAR(2) 生成器
    sim_var2_dynamic <- function(n, phi1_func, phi2_func, sigma_func) {
        x_out <- matrix(0, n, d)
        x_out[1:2, ] <- rmvnorm(2, mean = Zero, sigma = sigma_func(1))
        
        for (t in 3:n) {
            Phi1 <- phi1_func(t)
            Phi2 <- phi2_func(t)
            Sigma <- sigma_func(t)
            eps <- as.vector(rmvnorm(1, mean = Zero, sigma = Sigma))
            
            val <- (Phi1 %*% x_out[t-1, ]) + (Phi2 %*% x_out[t-2, ]) + eps
            x_out[t, ] <- as.vector(val)
        }
        return(x_out)
    }
    
    # BEKK (仅限二元)
    sim_bekk_process <- function(n, cp = n+1, mu_shift_vec = c(0,0), scale_factor = 1) {
        if(d != 2) stop("BEKK only supports d=2")
        x_out <- matrix(0, n, 2)
        H <- vector("list", n)
        H[[1]] <- 0.01 * diag(2)
        eps_base <- rmvnorm(n, mean = c(0,0), sigma = diag(2))
        x_out[1, ] <- t(chol(H[[1]])) %*% eps_base[1, ]
        
        for (i in 2:n) {
            is_post_change <- (i > cp)
            curr_mu <- if (is_post_change) mu_shift_vec else c(0,0)
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
    
    # --- Group 1: Null Hypothesis ---
    if (model <= 8) {
        if (model == 1) x <- rmvnorm(obs, sigma = Sigma_Base)
        else if (model == 2) x <- rmvt(obs, sigma = (3/5) * Sigma_Base, df = 5)
        # N3-N5: 一元模型
        else if (model == 3) x <- sim_var_static(obs, phi=0.1*I_d, sigma=I_d)
        else if (model == 4) x <- sim_var_static(obs, phi=0.3*I_d, sigma=I_d)
        else if (model == 5) x <- sim_var_static(obs, phi=0.5*I_d, sigma=I_d)
        # N6-N7: 二元模型
        else if (model == 6) x <- sim_var_static(obs, phi=Mat_B, sigma=I_d)
        else if (model == 7) x <- sim_var_static(obs, phi=Mat_A1, sigma=I_d)
        else if (model == 8) x <- sim_bekk_process(obs)
    }
    
    # --- Group 2: Mean Change (修改为截断式叠加) ---
    else if (model %in% c(101, 102, 103, 104)) {
        cp <- get_cp()
        d_val <- if(is.null(delta)) 1 else delta
        
        # [修改逻辑]: 先生成零均值基础数据，再叠加漂移
        
        if (model == 101) { # M1 (一元)
            # 1. 生成基础噪声 (N(0,1))
            x_base <- rmvnorm(obs, mean = Zero, sigma = I_d)
            # 2. 叠加均值
            x <- x_base
            x[(cp+1):obs, ] <- x[(cp+1):obs, ] + d_val # d=1, d_val为标量
            
        } else if (model == 102) { # M2 (二元)
            # 1. 生成基础 VAR 过程 (Weak VAR, 零均值)
            x_base <- sim_var_static(obs, phi = Mat_B, sigma = I_d)
            # 2. 叠加均值向量
            Mu_Vec <- rep(d_val, d) 
            x <- x_base
            # 对每一行加上 Mu_Vec
            # t() 确保 broadcasting 正确
            x[(cp+1):obs, ] <- t(t(x[(cp+1):obs, ]) + Mu_Vec)
            
        } else if (model == 103) { # M3 (二元)
            # 1. 生成基础 VAR 过程 (Medium VAR, 零均值)
            x_base <- sim_var_static(obs, phi = Mat_A1, sigma = I_d)
            # 2. 叠加均值向量
            Mu_Vec <- rep(d_val, d)
            x <- x_base
            x[(cp+1):obs, ] <- t(t(x[(cp+1):obs, ]) + Mu_Vec)
            
        } else if (model == 104) { # M4 (二元 BEKK)
            # BEKK 内部逻辑已经包含了 shift，此处维持原样或按需修改
            # 这里的 sim_bekk_process 已经支持在误差项上加 mu_shift_vec
            Mu_Vec <- c(0, d_val)
            x <- sim_bekk_process(obs, cp = cp, mu_shift_vec = Mu_Vec)
        }
    }
    
    # --- Group 3: Correlation Change (全二元) ---
    else if (model %in% c(201, 202, 203)) {
        cp <- get_cp()
        d_val <- if(is.null(delta)) 0.9 else delta
        Sigma_New <- matrix(c(1, d_val, d_val, 1), 2, 2)
        
        if (model == 201) {
            x <- rbind(rmvnorm(cp, sigma = I_d), rmvnorm(obs-cp, sigma = Sigma_New))
        } else if (model == 202) {
            x <- sim_var_dynamic(obs, 
                                 phi_func = function(t) Mat_B,
                                 sigma_func = function(t) if(t > cp) Sigma_New else I_d)
        } else if (model == 203) {
            x <- sim_var_dynamic(obs, 
                                 phi_func = function(t) Mat_A1,
                                 sigma_func = function(t) if(t > cp) Sigma_New else I_d)
        }
    }
    
    # --- Group 4: Variance Change (大部分二元) ---
    else if (model %in% c(211, 212, 213, 214)) {
        cp <- get_cp()
        
        if (model == 214) {
            d_val <- if(is.null(delta)) 1 else delta
            scale_factor <- sqrt(1 + d_val)
            x <- sim_bekk_process(obs, cp = cp, scale_factor = scale_factor)
        } else {
            d_val <- if(is.null(delta)) 2 else delta
            scale_factor <- sqrt(1 + d_val)
            Sigma_New <- (scale_factor^2) * I_d 
            
            if (model == 211) {
                x <- rbind(rmvnorm(cp, sigma = I_d), rmvnorm(obs-cp, sigma = Sigma_New))
            } else if (model == 212) {
                x <- sim_var_dynamic(obs, 
                                     phi_func = function(t) Mat_B,
                                     sigma_func = function(t) if(t > cp) Sigma_New else I_d)
            } else if (model == 213) {
                x <- sim_var_dynamic(obs, 
                                     phi_func = function(t) Mat_A1,
                                     sigma_func = function(t) if(t > cp) Sigma_New else I_d)
            }
        }
    }
    
    # --- Group 5: Structural Change (混合) ---
    else if (model == 301) { # T1 (一元)
        cp <- get_cp()
        x <- sim_var_dynamic(obs, 
                             phi_func = function(t) if(t > cp) 0.8*I_d else 0.1*I_d,
                             sigma_func = function(t) I_d)
        
    } else if (model == 302) { # T2 (二元)
        cp <- get_cp()
        x <- sim_var_dynamic(obs, 
                             phi_func = function(t) if(t > cp) Mat_A2 else Mat_A1,
                             sigma_func = function(t) I_d)
        
    } else if (model == 303) { # T3 (一元)
        cp <- get_cp()
        get_sigma_t3 <- function(rho) { (1 - rho^2) * I_d }
        x <- sim_var_dynamic(obs,
                             phi_func = function(t) if(t > cp) 0.7*I_d else 0.1*I_d,
                             sigma_func = function(t) if(t > cp) get_sigma_t3(0.7) else get_sigma_t3(0.1))
        
    } else if (model == 304) { # T4 (二元)
        cp <- get_cp()
        x <- sim_var_dynamic(obs,
                             phi_func = function(t) if(t > cp) -Mat_A_Flip else Mat_A_Flip,
                             sigma_func = function(t) I_d)
        
    } else if (model == 305) { # T5 (一元)
        cp <- get_cp()
        get_sigma_t5 <- function(rho) { (1 - rho^2) * I_d }
        x <- sim_var2_dynamic(obs,
                              phi1_func = function(t) 0*I_d, 
                              phi2_func = function(t) if(t > cp) 0.7*I_d else 0.1*I_d,
                              sigma_func = function(t) if(t > cp) get_sigma_t5(0.7) else get_sigma_t5(0.1))
        
    } else if (model == 306) { # T6 (二元)
        cp <- get_cp()
        x <- sim_var2_dynamic(obs,
                              phi1_func = function(t) matrix(0, 2, 2),
                              phi2_func = function(t) if(t > cp) -Mat_A_Flip else Mat_A_Flip,
                              sigma_func = function(t) I_d)
    } else {
        stop(paste("Error: Model ID", model, "not found."))
    }
    
    # 去除 burn-in
    return(x[-(1:burnin), , drop = FALSE])
}