######this file is the code of selection of m.
##### Perturbation Subsample
rm(list=ls())
setwd("E:\\HS\\OSRR\\code")

#### perturbation subsampling with product weights
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
  
  G <- generate(beta_0, N, "Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  
  beta_1 <- mulwrfit2(N,p,y,x,1,r)
  beta_5 <- mulwrfit2(N,p,y,x,5,r)
  beta_10 <- mulwrfit2(N,p,y,x,10,r)
  beta_15 <- mulwrfit2(N,p,y,x,15,r)
  beta_20 <- mulwrfit2(N,p,y,x,20,r)
  beta_25 <- mulwrfit2(N,p,y,x,25,r)
  beta_30 <- mulwrfit2(N,p,y,x,30,r)

  mse_1 <- MSE(beta_0, beta_1)
  mse_5 <- MSE(beta_0, beta_5)
  mse_10 <- MSE(beta_0, beta_10)
  mse_15 <- MSE(beta_0, beta_15)
  mse_20 <- MSE(beta_0, beta_20)
  mse_25 <- MSE(beta_0, beta_25)
  mse_30 <- MSE(beta_0, beta_30)
  c(mse_1, mse_5,mse_10,mse_15,mse_20,mse_25,mse_30)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE_1", "MSE_5","MSE_10","MSE_15","MSE_20","MSE_25","MSE_30")

##### Perturbation Subsampling with additive weights
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

H <- 200

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  beta_full <- rfit(y ~ x, data = data.frame(y, x))$coefficients[-1]
  mse <- MSE(beta_0, beta_full)
  c(mse)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE_1", "MSE_5","MSE_10","MSE_15","MSE_20","MSE_25","MSE_30")

