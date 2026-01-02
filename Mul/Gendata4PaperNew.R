library(mvtnorm)

#' 生成多元/一元时间序列模拟数据 (均值变化改为截断式叠加)
#'
#' @param monitoringPeriod 监控窗口比率 (L)
#' @param trainSize 训练样本长度 (T)
#' @param model 模型编号
#' @param cptLocation 变点位置 (0-1)
#' @param delta 变化幅度
#'
#' @return matrix (T_total x d)
gendataPaper <- function(monitoringPeriod, trainSize, model, cptLocation = NULL, delta = NULL) {
    
    # ==========================================
    # 1. 基础设置与维度判断
    # ==========================================
    
    univariate_models <- c(3, 4, 5, 101, 301, 303, 305)
    is_univariate <- model %in% univariate_models
    internal_d <- if(is_univariate) 1 else 2
    
    burnin <- ceiling(trainSize/2)
    len_monitoring <- floor(monitoringPeriod * trainSize)
    obs <- burnin + trainSize + len_monitoring
    
    # --- 确定变点位置 ---
    base_idx <- burnin + trainSize
    if (is.null(cptLocation)) {
        cp <- floor(base_idx + (len_monitoring * runif(1, 0, 0.8)))
    } else {
        if(cptLocation <= 0 || cptLocation >= 1) stop("cptLocation must be between 0 and 1.")
        cp <- floor(base_idx + (len_monitoring * cptLocation))
    }
    
    # --- 基础矩阵/向量预定义 ---
    I_d <- diag(internal_d)
    Zero <- rep(0, internal_d)
    
    if (!is_univariate) {
        Sigma_Base <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
        Chol_Base  <- t(chol(Sigma_Base)) 
        
        Mat_B  <- matrix(c(0.1, 0.05, 0.05, 0.1), 2, 2)
        Mat_A1 <- matrix(c(0.5, 0.2, 0.2, 0.1), 2, 2)
        Mat_A2 <- matrix(c(0.8, 0.3, 0.1, 0.7), 2, 2)
        Mat_A_Flip <- matrix(c(0.5, 0.1, 0.1, 0.5), 2, 2)
        
        C_bekk <- 0.001 * matrix(c(4, 0, 5, 3), 2, 2)
        A_bekk <- matrix(c(0.254, 0.040, -0.004, 0.332), 2, 2)
        B_bekk <- matrix(c(0.941, -0.019, 0.023, 0.864), 2, 2)
        Intercept_bekk <- t(C_bekk) %*% C_bekk
    }
    
    # ==========================================
    # 2. 核心生成器
    # ==========================================
    
    sim_var_fast <- function(n, cp, phi_pre, phi_post, 
                             sigma_chol_pre, sigma_chol_post, 
                             mu_pre, mu_post) {
        
        x_out <- matrix(0, n, internal_d)
        Z <- matrix(rnorm(n * internal_d), n, internal_d)
        
        val_init <- (sigma_chol_pre %*% Z[1, ]) + mu_pre
        x_out[1, ] <- as.numeric(val_init)
        
        end_pre <- min(cp, n)
        if (end_pre >= 2) {
            for (t in 2:end_pre) {
                eps <- sigma_chol_pre %*% Z[t, ]
                val <- (phi_pre %*% x_out[t-1, ]) + eps + mu_pre
                x_out[t, ] <- as.numeric(val)
            }
        }
        
        if (cp < n) {
            for (t in (cp + 1):n) {
                eps <- sigma_chol_post %*% Z[t, ]
                val <- (phi_post %*% x_out[t-1, ]) + eps + mu_post
                x_out[t, ] <- as.numeric(val)
            }
        }
        return(x_out)
    }
    
    sim_var2_fast <- function(n, cp, phi1_pre, phi1_post, 
                              phi2_pre, phi2_post, 
                              sigma_chol_pre, sigma_chol_post) {
        x_out <- matrix(0, n, internal_d)
        Z <- matrix(rnorm(n * internal_d), n, internal_d)
        
        x_out[1, ] <- as.numeric(sigma_chol_pre %*% Z[1, ])
        x_out[2, ] <- as.numeric(sigma_chol_pre %*% Z[2, ])
        
        end_pre <- min(cp, n)
        if (end_pre >= 3) {
            for (t in 3:end_pre) {
                eps <- sigma_chol_pre %*% Z[t, ]
                val <- (phi1_pre %*% x_out[t-1, ]) + (phi2_pre %*% x_out[t-2, ]) + eps
                x_out[t, ] <- as.numeric(val)
            }
        }
        
        if (cp < n) {
            for (t in (cp + 1):n) {
                eps <- sigma_chol_post %*% Z[t, ]
                val <- (phi1_post %*% x_out[t-1, ]) + (phi2_post %*% x_out[t-2, ]) + eps
                x_out[t, ] <- as.numeric(val)
            }
        }
        return(x_out)
    }
    
    sim_bekk_fast <- function(n, cp, mu_shift, scale_factor) {
        if(internal_d != 2) stop("BEKK defined for bivariate only")
        x_out <- matrix(0, n, 2)
        Z <- matrix(rnorm(n * 2), n, 2)
        I2 <- diag(2)
        H_curr <- 0.01 * I2 
        x_out[1, ] <- as.numeric(t(chol(H_curr)) %*% Z[1, ])
        for (t in 2:n) {
            if (t > cp) {
                z_t <- (Z[t, ] * scale_factor) + mu_shift
            } else {
                z_t <- Z[t, ]
            }
            prev_x <- matrix(x_out[t-1, ], ncol=1) 
            term_ARCH <- t(A_bekk) %*% (prev_x %*% t(prev_x)) %*% A_bekk
            term_GARCH <- t(B_bekk) %*% H_curr %*% B_bekk
            H_curr <- Intercept_bekk + term_ARCH + term_GARCH
            x_out[t, ] <- as.numeric(t(chol(H_curr)) %*% z_t)
        }
        return(x_out)
    }
    
    return_data <- function(data) {
        data[-(1:burnin), , drop=FALSE]
    }
    
    val_delta <- if(is.null(delta)) 1.0 else delta
    
    # ==========================================
    # 3. 模型生成 (包含截断式均值修改)
    # ==========================================
    
    # --- 1. Null Hypothesis ---
    if (model <= 8) {
        if (model == 1) return(return_data(rmvnorm(obs, sigma = Sigma_Base)))
        else if (model == 2) return(return_data(rmvt(obs, sigma = (3/5) * Sigma_Base, df = 5)))
        else if (model == 8) return(return_data(sim_bekk_fast(obs, cp=obs+1, mu_shift=c(0,0), scale_factor=1)))
        
        phi <- if(model==3) 0.1*I_d else if(model==4) 0.3*I_d else if(model==5) 0.5*I_d else if(model==6) Mat_B else Mat_A1
        
        return(return_data(sim_var_fast(obs, cp=obs, phi_pre=phi, phi_post=phi, 
                                        sigma_chol_pre=I_d, sigma_chol_post=I_d, 
                                        mu_pre=Zero, mu_post=Zero)))
    }
    
    # --- 2. Mean Change (M) - 截断式修改版 ---
    else if (model %in% 101:104) {
        # M1: i.i.d. Univariate Gaussian
        if(model == 101) {
            x <- matrix(rnorm(obs), ncol=1)
            if(cp < obs) x[(cp+1):obs, 1] <- x[(cp+1):obs, 1] + val_delta
            return(return_data(x))
        }
        
        mu_vec <- if(model==104) c(val_delta, val_delta) else c(val_delta, val_delta)
        # c(0, val_delta) else c(val_delta, val_delta)
        # M4: BEKK (BEKK本质上是波动模型，这里保持其内部生成的偏移逻辑)
        if(model == 104) return(return_data(sim_bekk_fast(obs, cp, mu_shift=mu_vec, scale_factor=1)))
        
        # M2 & M3: VAR(1) 截断式生成
        # 逻辑：先生成完整的平稳过程，然后在变点后直接加 delta
        phi <- if(model==102) Mat_B else Mat_A1
        x <- sim_var_fast(obs, cp=obs, phi_pre=phi, phi_post=phi,
                          sigma_chol_pre=I_d, sigma_chol_post=I_d,
                          mu_pre=Zero, mu_post=Zero)
        
        if(cp < obs) {
            # 对变点之后的所有行加上均值向量
            x[(cp+1):obs, ] <- sweep(x[(cp+1):obs, , drop=FALSE], 2, mu_vec, "+")
        }
        return(return_data(x))
    }
    
    # --- 3. Correlation Change (C) ---
    else if (model %in% 201:203) {
        d_corr <- if(is.null(delta)) 0.9 else delta
        Sigma_New <- matrix(c(1, d_corr, d_corr, 1), 2, 2)
        Chol_New <- t(chol(Sigma_New))
        
        if(model == 201) {
            Z <- matrix(rnorm(obs*2), obs, 2)
            x <- matrix(0, obs, 2)
            x[1:cp, ] <- Z[1:cp, ] 
            x_part2 <- t(Chol_New %*% t(Z[(cp+1):obs, ]))
            x[(cp+1):obs, ] <- x_part2
            return(return_data(x))
        }
        
        phi <- if(model==202) Mat_B else Mat_A1
        return(return_data(sim_var_fast(obs, cp, phi_pre=phi, phi_post=phi,
                                        sigma_chol_pre=I_d, sigma_chol_post=Chol_New,
                                        mu_pre=Zero, mu_post=Zero)))
    }
    
    # --- 4. Variance Change (V) ---
    else if (model %in% 211:214) {
        d_var <- if(is.null(delta)) { if(model==214) 1 else 2 } else delta
        scale <- sqrt(1 + d_var)
        
        if(model == 214) return(return_data(sim_bekk_fast(obs, cp, mu_shift=c(0,0), scale_factor=scale)))
        
        Chol_Post <- scale * I_d
        
        if(model == 211) {
            x_pre <- rmvnorm(cp, sigma = I_d)
            x_post <- rmvnorm(obs-cp, sigma = (scale^2)*I_d)
            return(return_data(rbind(x_pre, x_post)))
        }
        
        phi <- if(model==212) Mat_B else Mat_A1
        return(return_data(sim_var_fast(obs, cp, phi_pre=phi, phi_post=phi,
                                        sigma_chol_pre=I_d, sigma_chol_post=Chol_Post,
                                        mu_pre=Zero, mu_post=Zero)))
    }
    
    # --- 5. Structural Change (T) ---
    else if (model %in% 301:306) {
        if(model == 301) return(return_data(sim_var_fast(obs, cp, phi_pre=0.1*I_d, phi_post=0.8*I_d, I_d, I_d, Zero, Zero)))
        if(model == 302) return(return_data(sim_var_fast(obs, cp, phi_pre=Mat_A1, phi_post=Mat_A2, I_d, I_d, Zero, Zero)))
        if(model == 304) return(return_data(sim_var_fast(obs, cp, phi_pre=Mat_A_Flip, phi_post=-Mat_A_Flip, I_d, I_d, Zero, Zero)))
        
        get_chol <- function(rho) t(chol((1-rho^2)*I_d))
        
        if(model == 303) return(return_data(sim_var_fast(obs, cp, 
                                                         phi_pre=0.1*I_d, phi_post=0.7*I_d, 
                                                         sigma_chol_pre=get_chol(0.1), sigma_chol_post=get_chol(0.7), 
                                                         Zero, Zero)))
        
        if(model == 305) return(return_data(sim_var2_fast(obs, cp, 
                                                          0*I_d, 0*I_d, 
                                                          0.1*I_d, 0.7*I_d, 
                                                          get_chol(0.1), get_chol(0.7))))
        
        if(model == 306) return(return_data(sim_var2_fast(obs, cp, 
                                                          matrix(0,2,2), matrix(0,2,2),
                                                          Mat_A_Flip, -Mat_A_Flip,
                                                          I_d, I_d)))
    }
    
    stop(paste0("错误: 找不到模型 ID ", model, "。\n",
                "请输入有效的模型编号范围:\n",
                " - [1-8]: 无结构变化 (Null Hypothesis)\n",
                " - [101-104]: 均值变化 (Mean Change)\n",
                " - [201-203]: 相关系数变化 (Correlation Change)\n",
                " - [211-214]: 方差变化 (Variance Change)\n",
                " - [301-306]: 依赖结构/联合分布变化 (Structural Change)"))
}