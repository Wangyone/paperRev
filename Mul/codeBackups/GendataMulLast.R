#########################################################
########## Data generating process ########
#########################################################
genDataMul <- function(L,T,model) {
  library(mvtnorm)
  library(MTS)
  d <- 2 # Dimension
  ### burn-in(烧链）
  burnin <- ceiling(T/2)
  ## Expected sample size
  N <- (1+L)*T
  
  obs <- burnin + N
  ## parameter
  delta <- 0
  # Different correlation coefficients
  kor0 <- 0
  kor1<-0.2
  kor2<-0.6
  kor3<-0.9
  mu0=rep(0,2)
  mu1=rep(1,2)
  ##  covariance matrix
  sigma0 <- matrix(c(1,kor0,kor0,1),nrow<-2)
  sigma1 <- matrix(c(1,kor1,kor1,1),nrow<-2)
  sigma2 <- matrix(c(1,kor2,kor2,1),nrow<-2)
  sigma3 <- matrix(c(1,kor3,kor3,1),nrow<-2)
  B <- matrix(c(0.2,0.1,0.1,0.2),2,2)/2
  ## VAR(1) coefficients 
  A1 <- matrix(c(0.5,0.2,0.2,0.1),2,2)
  A2 <- matrix(c(0.8,0.1,0.3,0.7),2,2)
  ## BEKK-GARCH(1,1)
  C <- t(matrix(c(0.004,0.005,0,0.003),nrow<-2))
  Amat <- matrix(c(0.254,0.004,-0.004,0.332),nrow<-2)
  Bmat <- matrix(c(0.941,-0.019,0.023,0.864),nrow<-2)
  intercept <- t(C)%*%(C)
  ## If the absolute values of the matrix eigenvalues are all less than 1, 
  #then the VAR (1) model is stationary.
  reps0 <-  function(n,kr=kor) {rmvnorm(n, mean = rep(0,2), sigma = matrix(c(1,kr,kr,1),nrow=2))}
  
  reps1 <-  function(n,kr=kor) {rmvt(n, sigma = matrix(c(1,kr,kr,1),nrow = 2), df = 5)} ## for model 304
  
  reps2 <-  function(n,kr=kor) {rmvt(n, sigma = (3/5)*matrix(c(1,kr,kr,1),nrow = 2), df = 5)} #t5(0,sigma)
  
  x <- matrix(0,obs,d)
  epsN <- reps0(obs,kor0)
  epsT <- reps1(obs,kor2) ## for model 304
  
  if (model == 1) {           ##### DGP N1 #######
    x <- epsN
  } 
  else if (model == 2) {      ##### DGP N2 #######
    x<-reps0(obs,kor2)
  }
  else if (model == 3) {      ##### DGP N3 #######
    x <- rmvt(obs, sigma =(3/5)*sigma2, df = 5)    
    # sigma =(3/5)*sigma2
    # t5(0,Sigma2)
    # Ensure the generation of accurate data
  }
  else if (model == 4) {      
    ##### DGP N4 #######  
    x<-ts(VARMAsim(obs,arlags=c(1),phi=0.1*diag(2),sigma=sigma0)$series)
  }
  else if (model == 5) {      ##### DGP N5 #######  
    x<-ts(VARMAsim(obs,arlags=c(1),phi=0.3*diag(2),sigma=sigma0)$series)
  }
  else if (model == 6) {      ##### DGP N6 #######  
    x<-ts(VARMAsim(obs,arlags=c(1),phi=0.5*diag(2),sigma=sigma0)$series)
  }
  else if (model == 7) {      ##### DGP N7 #######  
    B <- matrix(c(0.2,0.1,0.1,0.2),2,2)/2
    #eigen(B)$values
    x<-ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
  }
  else if (model == 8) {      ##### DGP N8 #######  
    x<-ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=sigma1)$series)
  }
  else if (model == 9) {      ##### DGP N9 #######  
    # BEKK-GARCH(1,1)
    eps1 <- reps0(obs,kor0)
    x <- matrix(0,obs,2)
    H <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    H[[1]] <- 0.01*diag(2)      # initial 
    x[1,] <- t(chol(H[[1]]))%*%eps1[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( x[(i-1),])
      Inn <- Ri%*%t(Ri)  
      H[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%H[[(i-1)]]%*%Bmat
      h12 <-  t(chol(H[[i]]))    
      x[i,] <- h12%*%eps1[i,]
    }
  }
  
  ######################## mean change ################################
  else if (model == 101) {     ##### DGP M1 δ=0.5 ####### 
    delta <- 0.5
    x <- epsN
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
    
  } 
  else if (model == 102) {    #### DGP M1 δ=1 ####### 
    delta <- 1
    x <- epsN
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
    
  } 
  else if (model == 103) {     ##### DGP M2 δ=0.5 ####### 
    delta <- 0.5
    x<-ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  else if (model == 104) {     ##### DGP M2 δ=1 ####### 
    delta <- 1
    x<-ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  else if (model == 105) {     ##### DGP M3 δ=0.5 ####### 
    delta <- 0.5
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=sigma1)$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  else if (model == 106) {     ##### DGP M3 δ=1 ####### 
    delta <- 1
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=sigma1)$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  else if (model == 107) {     ##### DGP M4 δ=0.5 ####### 
    delta <- 0.5
    ### before change
    eps1 <- rmvnorm(obs, mean = c(0,0), sigma = diag(2) )
    x <- matrix(0,obs,2)
    H <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    H[[1]] <- 0.01*diag(2)      # initial 
    
    x[1,] <- t(chol(H[[1]]))%*%eps1[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( x[(i-1),])
      Inn <- Ri%*%t(Ri)  
      H[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%H[[(i-1)]]%*%Bmat
      h12 <-  t(chol(H[[i]]))    
      x[i,] <- h12%*%eps1[i,]
    }
    ### after change
    eps2 <- rmvnorm(obs, mean = c(delta,0), sigma = diag(2) )
    y <- matrix(0,obs,2)
    Hy <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    Hy[[1]] <- 0.01*diag(2)      # initial 
    y[1,] <- t(chol(Hy[[1]]))%*%eps2[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( y[(i-1),])
      Inn <- Ri%*%t(Ri)  
      Hy[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%Hy[[(i-1)]]%*%Bmat
      h12 <-  t(chol(Hy[[i]]))    
      y[i,] <- h12%*%eps2[i,]
    }
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- y[change.pt:obs,] 
  } 
  else if (model == 108) {     ##### DGP M4 δ=1 ####### 
    delta <- 1
    ### before change
    eps1 <- rmvnorm(obs, mean = c(0,0), sigma = diag(2) )
    x <- matrix(0,obs,2)
    H <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    H[[1]] <- 0.01*diag(2)      # initial 
    
    x[1,] <- t(chol(H[[1]]))%*%eps1[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( x[(i-1),])
      Inn <- Ri%*%t(Ri)  
      H[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%H[[(i-1)]]%*%Bmat
      h12 <-  t(chol(H[[i]]))    
      x[i,] <- h12%*%eps1[i,]
    }
    ### after change
    eps2 <- rmvnorm(obs, mean = c(delta,0), sigma = diag(2) )
    y <- matrix(0,obs,2)
    Hy <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    Hy[[1]] <- 0.01*diag(2)      # initial 
    y[1,] <- t(chol(Hy[[1]]))%*%eps2[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( y[(i-1),])
      Inn <- Ri%*%t(Ri)  
      Hy[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%Hy[[(i-1)]]%*%Bmat
      h12 <-  t(chol(Hy[[i]]))    
      y[i,] <- h12%*%eps2[i,]
    }
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- y[change.pt:obs,] 
  } 
  
################# second-order moment change ##################
  
  ########################## Change in  Variance
  else if (model == 201) {   ##### DGP V1 δ=1  Variance    #######   
    delta <- 1
    x <- rmvnorm(obs, mean = rep(0,2), sigma = diag(2))
    y <-rmvnorm(obs, mean = rep(0,2), sigma = (1+delta)*diag(2))
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 202) {  ##### DGP V1 δ=2  Variance #######   
    delta <- 2
    x <- rmvnorm(obs, mean = rep(0,2), sigma = diag(2))
    y <-rmvnorm(obs, mean = rep(0,2), sigma = (1+delta)*diag(2) )
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  
  else if (model == 203) {   ##### DGP V2 δ=1   Variance   #######  
    delta <- 1
    x <- ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=(1+delta)*diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 204) {   ##### DGP V2 δ=2 Variance    #######   
    delta <- 2
    x <- ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=(1+delta)*diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 205) {   ##### DGP V3 δ=1  Variance   #######   
    delta <- 1
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=diag(2) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=(1+delta)*diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 206) {   ##### DGP V3 δ=2  Variance   #######   
    delta <- 2
    x <-  ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=diag(2) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=(1+delta)*diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  
  else if (model == 207) {   ##### DGP V4 δ=1  Variance   #######   
    delta <- 1
    ### before change
    eps1 <- rmvnorm(obs, mean = rep(0,2), sigma =diag(2)  )
    x <- matrix(0,obs,2)
    H <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    H[[1]] <- 0.01*diag(2)      # initial 
    x[1,] <- t(chol(H[[1]]))%*%eps1[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( x[(i-1),])
      Inn <- Ri%*%t(Ri)  
      H[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%H[[(i-1)]]%*%Bmat
      h12 <-  t(chol(H[[i]]))    
      x[i,] <- h12%*%eps1[i,]
    }
    ### after change   matrix(c(1+delta,0,0,1),nrow<-2)
    eps2 <-  rmvnorm(obs, mean = rep(0,2), sigma =matrix(c(1+delta,0,0,1),nrow<-2)  )
    y <- matrix(0,obs,2)
    Hy <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    Hy[[1]] <- 0.01*diag(2)      # initial 
    y[1,] <- t(chol(Hy[[1]]))%*%eps2[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( y[(i-1),])
      Inn <- Ri%*%t(Ri)  
      Hy[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%Hy[[(i-1)]]%*%Bmat
      h12 <-  t(chol(Hy[[i]]))    
      y[i,] <- h12%*%eps2[i,]
    }
    
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 208) {   ##### DGP V5 δ=1   Variance  #######   
    delta <- 1
    ### before change
    eps1 <- rmvnorm(obs, mean = rep(0,2), sigma =matrix(c(1,0,0,1),nrow<-2)  )
    
    x <- matrix(0,obs,2)
    H <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    H[[1]] <- 0.01*diag(2)      # initial 
    x[1,] <- t(chol(H[[1]]))%*%eps1[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( x[(i-1),])
      Inn <- Ri%*%t(Ri)  
      H[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%H[[(i-1)]]%*%Bmat
      h12 <-  t(chol(H[[i]]))    
      x[i,] <- h12%*%eps1[i,]
    }
    ### after change
    eps2 <-  rmvnorm(obs, mean = rep(0,2), sigma =matrix(c(1+delta,0,0,1+delta),nrow<-2)  )
    y <- matrix(0,obs,2)
    Hy <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    Hy[[1]] <- 0.01*diag(2)      # initial 
    y[1,] <- t(chol(Hy[[1]]))%*%eps2[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( y[(i-1),])
      Inn <- Ri%*%t(Ri)  
      Hy[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%Hy[[(i-1)]]%*%Bmat
      h12 <-  t(chol(Hy[[i]]))    
      y[i,] <- h12%*%eps2[i,]
    }
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  
  ########################## Change in  Correlation
  else if (model == 211) {  ##### DGP C1 ρ=0.9 Correlation #######   
    x <-  reps0(obs,kor0)
    y <- reps0(obs,kor3)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  } 
  else if (model == 212) {   ##### DGP C2 ρ=0.9 Correlation    #######   
    x <- ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=sigma3 )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 213) {  ##### DGP C3 ρ=0.9 Correlation #######   
    x<-   ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=sigma0 )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=sigma3 )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  } 
  ######################## temporal dependence changes ################################ 
  else if (model == 301) { ##### DGP T1 δ=0.8 #######  
    delta <- 0.8
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ### before the change point
    
    xt <- ts(VARMAsim(obs,arlags=c(1),phi=0.1*diag(2),sigma=sigma0)$series)
    ### after the change point  
    
    yt <- ts(VARMAsim(obs,arlags=c(1),phi=(delta)*diag(2),sigma=sigma0)$series)
    
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  else if (model == 302) {  ##### DGP T2 VAR(1) #######  
    
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    
    xt <- ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=sigma1)$series) ### before the change point
    
    yt <- ts(VARMAsim(obs,arlags=c(1),phi=A2,sigma=sigma1)$series)  ### after the change point  
    
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  ## reps1 <-  function(n,kr=kor) {rmvt(n, sigma = (5/3)*matrix(c(1,kr,kr,1),nrow = 2), df = 5)}
  
  else if (model == 303) {    ##### DGP T3  VAR(1) ####### 
    #t5(0,sigma)
    #reps2 <-  function(n,kr=kor) {rmvt(n, sigma = (3/5)*matrix(c(1,kr,kr,1),nrow = 2), df = 5)} 
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ## reps1 <-  function(n,kr=kor) {rmvt(n, sigma = matrix(c(1,kr,kr,1),nrow = 2), df = 5)}
    p <- 1
    d <- 2
    ### before the change point
    xt <- matrix(rep(NA,(obs+p)*d),nrow=(obs+p),ncol=d)
    xt[1:p,] <- reps1(p,kor1)
    for(i in (p+1):(obs)) {
      xt[i,]=t(xt[(i-p):(i-1),])%*%t(A1)+reps2(1,kor1)  
    }
    ### after the change point  
    yt <- matrix(rep(NA,(obs+p)*d),nrow=(obs+p),ncol=d)
    yt[1:p,] <- reps1(p,kor1)
    for(i in (p+1):(obs)) {
      yt[i,]=t(yt[(i-p):(i-1),])%*%t(A2)+reps2(1,kor1) 
    }
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  ######################## Joint distribution
  else if (model == 304) { ##### DGP J1 #######  
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ### before the change point
    xt <- matrix(0,obs,2)
    rhoXt <-  0.1
    xt[1,] <- rmvnorm(1, mean = rep(0,2), sigma = diag(2))
    for (i in 2:NROW(xt)){
      xt[i,] <- rhoXt*xt[i-1,]+rmvnorm(1, mean = rep(0,2), sigma =(1-rhoXt^2)*diag(2)) 
    }
    ### after the change point  
    yt <- matrix(0,obs,2)
    rhoYt <-  0.7
    yt[1,] <- rmvnorm(1, mean = rep(0,2), sigma = diag(2))
    for (i in 2:NROW(yt)){
      yt[i,] <- rhoYt*yt[i-1,]+rmvnorm(1, mean = rep(0,2), sigma =(1-rhoYt^2)*diag(2))
    } 
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  else if (model == 305) { ##### DGP J2 #######  
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ### before the change point
    xt <- reps0(obs,0)
    #xt <- ts(VARMAsim(obs,arlags=c(1),phi=0*diag(2),sigma=diag(2))$series)
    ### after the change point 
    yt <- matrix(0,obs,2)
    Phi <- matrix(c(0.8,0.2,-0.2,0.7),2,2) 
    abs( eigen(Phi)$values )
    I4 <- diag(4)
    pp <- kronecker(Phi,Phi) # Kronecker product
    dd <- I4-pp
    sig <- dd%*%matrix(diag(2),4,1) # Obtain sigma-et
    sigma <- matrix(sig,2,2)
    #is.positive.definite(sigma)
    yt <- ts(VARMAsim(obs,arlags=c(1),phi=Phi,sigma=sigma)$series)
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  else if (model == 306) { ##### DGP J3 #######  
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    ### before the change point
    xt <- reps0(obs,0)
    ### after the change point  
    yt <- matrix(0,obs,2)
    rhoYt <-  0.7
    yt[1:2,] <- rmvnorm(2, mean = rep(0,2), sigma = diag(2))
    for (i in 3:NROW(yt)){
      yt[i,] <- rhoYt*yt[i-2,]+rmvnorm(1, mean = rep(0,2), sigma =(1-rhoYt^2)*diag(2))
    } 
    x <- ts(rbind(xt[1:change.pt,], yt[(change.pt+1):obs,]))
  }
  x <- x[-(1:burnin),]
}

