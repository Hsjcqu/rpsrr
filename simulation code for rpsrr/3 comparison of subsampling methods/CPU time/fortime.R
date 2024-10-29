##### This file is used to calculate the CPU time
rm(list=ls())
source("main functions.R")

library(mvtnorm)
library(doParallel)
library(foreach)

cores <- 12
cl <- makeCluster(cores)
registerDoParallel(cl)

beta_0 <- c(-4, 0.11, -0.9, -0.5, 1)
N <- 10000
p <- length(beta_0)

H <- 1000
#####time for global estimator
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  # Calculate perturbed differences
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  t1 <- proc.time()
  betafull<- quantreg::rq(y_diff ~ 0 + X_diff, method = "fnb", tau = 0.5)$coefficients
  t2 <- proc.time()
  t <- t2 - t1
  c(t[1])
}
stopCluster(cl)
resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))

#####time for perturbation subsampling 
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  G <- generate(beta_0, N, "Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  # Calculate perturbed differences
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)

  atime1 <- sumtime(N,p,y_diff, X_diff, 50)
  atime2 <- sumtime(N,p,y_diff, X_diff, 150)
  atime3 <- sumtime(N,p,y_diff, X_diff, 250)
  ptime4 <- multime(N,p,y_diff, X_diff, 50)
  ptime5 <- multime(N,p,y_diff, X_diff, 150)
  ptime6 <- multime(N,p,y_diff, X_diff, 250)
  c(atime1,atime2,atime3,ptime4,ptime5,ptime6)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))