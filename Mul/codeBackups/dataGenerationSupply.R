dataGeneration <- function(L,T,model,dataDimension) {
  library(mvtnorm)
  library(MTS) # Using VARMAsim fuction
  ### Data length
  burnin <- ceiling(T/2)
  N <- (1 + L)*T
  obs <- burnin + N
  ## parameter
  rho0 <- 0
  mu0=rep(0,dataDimension)
  # 有20%的分量均值变化
  mu1=c( rep(1,ceiling(dataDimension/5)),rep(0,(dataDimension-ceiling(dataDimension/5))  ) )
  ##  covariance matrix
  createCormat <- function(d, s) {
    mat <- matrix(NA, nrow = d, ncol = d)
    # ifelse 判断是否在对角线上
    matrix(ifelse(row(mat) == col(mat), 1, s), nrow = d)
    }
  ## 系数矩阵PHI
  createPHI <- function(dimensionPHI, s) {
    mat <- matrix(0, nrow = dimensionPHI, ncol = dimensionPHI)
    diag(mat) <- 1
    # 只有当维度大于1时才设置非对角线元素
    if(dimensionPHI > 1) {
      for(i in 1:(dimensionPHI - 1)) {
        mat[i, i + 1] <- s    # 上对角线
        mat[i + 1, i] <- s    # 下对角线
      }
    }
    return(mat)
  }
  
 
  sigma0 <- createCormat(dataDimension,rho0)
  # 生成误差项
  reps0 <-  function(n,kr){
    rmvnorm(n, mean = rep(0,dataDimension),
            sigma = createCormat(dataDimension,kr) )}
  ## VAR(1)模型的系数矩阵
  A4 <- createPHI(dataDimension, 0.4)*0.1
  A5 <- createPHI(dataDimension, 0.2)*0.5
  ## 数据
  x <- matrix(0,obs,dataDimension)
  epsN <- reps0(obs,rho0)
  
  if (model == 1) {           ##### DGP N1 #######
    x <- epsN
  } 
  else if (model == 2) {      ##### DGP N2 #######
    x <- ts(VARMAsim(obs,arlags=c(1),phi=0.3*diag(dataDimension),sigma=sigma0)$series)
  }
  else if (model == 3) {      ##### DGP N3 #######
    x <- ts(VARMAsim(obs,arlags=c(1),phi=0.5*diag(dataDimension),sigma=sigma0)$series)
  }
  else if (model == 4) {      ##### DGP N4 #######  VAR(1)
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A4,sigma=sigma0)$series)
  }
  else if (model == 5) {      ##### DGP N5 #######  VAR(1)
    x <- ts(MTS::VARMAsim(obs,arlags=c(1),phi = A5,sigma = sigma0)$series)
  }
  ######################## mean change ################################
  else if (model == 101) {     ##### DGP M1 δ=1 ####### 
    delta <- 1
    x <- epsN
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  else if (model == 102) {     ##### DGP M2 δ=1 ####### 
    delta <- 1
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A4,sigma=sigma0)$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  else if (model == 103) {     ##### DGP M2 δ=1 ####### 
    delta <- 1
    x<-ts(VARMAsim(obs,arlags=c(1),phi=A5,sigma=sigma0)$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  
  
  ######################## second-order moment change ################################ 
  
  ########################## Change in  Variance
  else if (model == 201) {   ##### DGP V1 δ=1     #######   
    delta <- 1
    x <- rmvnorm(obs, mean = mu0, sigma = diag(dataDimension))
    y <-rmvnorm(obs, mean = mu0, sigma = (1+delta)*diag(dataDimension))
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 202) {   ##### DGP V2 δ=1     #######   Variance
    delta <- 1
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A4,sigma=diag(dataDimension) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=A4,sigma=(1+delta)*diag(dataDimension) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }else if (model == 203) {   ##### DGP V2 δ=1     #######   Variance
    delta <- 1
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A5,sigma=diag(dataDimension) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=A5,sigma=(1+delta)*diag(dataDimension) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  ######################## temporal dependence changes ################################ 
  else if (model == 301) { ##### DGP T1 δ=0.8 #######  
    delta <- 0.8
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ### before the change point
    
    xt <- ts(VARMAsim(obs,arlags=c(1),phi=0.1*diag(dataDimension),sigma=sigma0)$series)
    ### after the change point  
    
    yt <- ts(VARMAsim(obs,arlags=c(1),phi=(delta)*diag(dataDimension),sigma=sigma0)$series)
    
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  ######################## Joint distribution
  else if (model == 302) { 
    ##### DGP J1 #######  
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ### before the change point
    xt <- matrix(0,obs,dataDimension)
    rhoXt <-  0.1
    xt[1,] <- rmvnorm(1, mean = mu0, sigma = diag(dataDimension))
    for (i in 2:NROW(xt)){
      xt[i,] <- rhoXt*xt[i-1,]+rmvnorm(1, mean = mu0, sigma =(1-rhoXt^2)*diag(dataDimension)) 
    }
    ### after the change point  
    yt <- matrix(0,obs,dataDimension)
    rhoYt <-  0.7
    yt[1,] <- rmvnorm(1, mean = mu0, sigma = diag(dataDimension))
    for (i in 2:NROW(yt)){
      yt[i,] <- rhoYt*yt[i-1,]+rmvnorm(1, mean = mu0, sigma =(1-rhoYt^2)*diag(dataDimension))
    } 
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  else if (model == 303) { ##### DGP J3 #######  
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ### before the change point
    xt <- reps0(obs,0)
    ### after the change point  
    yt <- matrix(0,obs,dataDimension)
    rhoYt <-  0.7
    yt[1:2,] <- rmvnorm(2, mean = mu0, sigma = diag(dataDimension))
    for (i in 3:NROW(yt)){
      yt[i,] <- rhoYt*yt[i-2,]+rmvnorm(1, mean = mu0, sigma =(1-rhoYt^2)*diag(dataDimension))
    } 
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  x <- x[-(1:burnin),]
}

