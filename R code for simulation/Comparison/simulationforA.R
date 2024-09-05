##### Perturbation Subsampling with additive weights in different cases.
###case 2
rm(list=ls())
setwd("E:\\DSF\\R")

library(mvtnorm)
library(doParallel)
library(foreach)

cores <- 40
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("function8.R")
})

beta_0 <- c(-4, 0.11, -0.9, -0.5, 1)
N <- 10000
p <- length(beta_0)

H <- 1000

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "t3") 
  x <- G[[2]]
  y <- G[[1]]
  # Calculate perturbed differences
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  beta1 <- sumwrfit(N,p,y_diff, X_diff, 10, 50,'pois')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 10, 100,'pois')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, 150,'pois')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 10, 200,'pois')
  mse4 <- MSE(beta_0, beta4)
  beta5 <- sumwrfit(N,p,y_diff, X_diff, 10, 250,'pois')
  mse5 <- MSE(beta_0, beta5)
  beta6 <- sumwrfit(N,p,y_diff, X_diff, 10, 300,'pois')
  mse6 <- MSE(beta_0, beta6)
  c(mse1,mse2,mse3,mse4,mse5,mse6)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE1", "MSE2","MSE3","MSE4","MSE5","MSE6")


##### Case 3
rm(list=ls())
setwd("E:\\DSF\\R")

library(mvtnorm)
library(doParallel)
library(foreach)

cores <- 40
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("function8.R")
})

beta_0 <- c(-4, 0.11, -0.9, -0.5, 1)
N <- 10000
p <- length(beta_0)

H <- 1000

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "Mixed") 
  x <- G[[2]]
  y <- G[[1]]
  # Calculate perturbed differences
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  # 扰动
  beta1 <- sumwrfit(N,p,y_diff, X_diff, 10, 50,'pois')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 10, 100,'pois')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, 150,'pois')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 10, 200,'pois')
  mse4 <- MSE(beta_0, beta4)
  beta5 <- sumwrfit(N,p,y_diff, X_diff, 10, 250,'pois')
  mse5 <- MSE(beta_0, beta5)
  beta6 <- sumwrfit(N,p,y_diff, X_diff, 10, 300,'pois')
  mse6 <- MSE(beta_0, beta6)
  c(mse1,mse2,mse3,mse4,mse5,mse6)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE1", "MSE2","MSE3","MSE4","MSE5","MSE6")