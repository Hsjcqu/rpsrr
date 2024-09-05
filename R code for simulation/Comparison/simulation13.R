######Case 3
##### Perturbation Subsampling with product weghts
rm(list=ls())
setwd("E:\\HS\\OSRR\\code")

library(mvtnorm)
library(Rfit)
library(doParallel)
library(foreach)


cores <- 16
cl <- makeCluster(cores)
registerDoParallel(cl)


clusterEvalQ(cl, {
  source("function8.R")
})

beta_0 <- c(-4, 0.11, -0.9, -0.5, 1)
N <- 10000
p <- length(beta_0)

H <- 1000
r <- 50

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  G <- generate(beta_0, N, "Normal", "Mixed") 
  x <- G[[2]]
  y <- G[[1]]
  
  
  beta_U <- optsample(y, x, r, 'U')
  beta_L <- optsample(y, x, r, 'L')
  mse_U <- MSE(beta_0, beta_U)
  mse_L <- MSE(beta_0, beta_L)
  
 
  M <- 10
  betaper_pois<- betaper_linear <- matrix(0, M, p)
  
  for (l in 1:M) {
    q <- r / N
    u <- rbinom(N, 1, q) 
    v_pois<- rpois(N, 1/q)
    betaper_pois[l,] <- persample(y, x, u, v_pois,"R") 
    betaper_linear[l,] <- persample(y, x, u, v_pois,"L") 
  }
  beta_pois <- colMeans(betaper_pois)
  mse_pois <- MSE(beta_0, beta_pois)
  beta_linear <- colMeans(betaper_linear)
  mse_linear <- MSE(beta_0, beta_linear)
  
  c(mse_U, mse_L,mse_linear, mse_pois)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE_U", "MSE_L","MSE_Linear", "MSE_pois")