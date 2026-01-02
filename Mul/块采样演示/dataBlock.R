dataBlock <- function(sampleSize,beta,dataDimension) {
  library(mvtnorm)
  library(MTS) 
  # Using VARMAsim fuction
  ### Data length
  burnin <- ceiling(sampleSize/2)
  Numberdata <- burnin + sampleSize
  ##  covariance matrix
  createCormat <- function(d, s) {
    mat <- matrix(NA, nrow = d, ncol = d)
    # ifelse 判断是否在对角线上
    matrix(ifelse(row(mat) == col(mat), 1, s), nrow = d)}
  
  x<-ts(VARMAsim(Numberdata,arlags=c(1),phi=beta*diag(dataDimension),
                 sigma = diag(dataDimension) )$series)
  
  yt <- x[-c(1:burnin),]
  
  return(yt)
}


