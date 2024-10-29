#####Case 2
rm(list=ls())

library(mvtnorm)
library(Rfit)
library(doParallel)
library(foreach)

cores <- 16
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("main functions.R")
})

#beta_0 <- c(-4, 0.11, -0.9, -0.5) # for 5-dim setting
beta_0 <- c(-4, 0.11, -0.9, -0.5, 1,1,1,1,1,1)  # for 10-dim setting
N <- 10000
p <- length(beta_0)

H <- 1000

###### global estimator
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "t3") 
  x <- G[[2]]
  y <- G[[1]]
  beta_full <- rfit(y ~ x, data = data.frame(y, x))$coefficients[-1]
  mse_full <- MSE(beta_0,beta_full)
  mse_full
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))



##### r \in {50,100,150,200,250,300}
r <- 300 #change r to get the results

#####product weights
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "t3") 
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


######additive weightes
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  # 生成数据
  G <- generate(beta_0, N, "Normal", "t3") 
  x <- G[[2]]
  y <- G[[1]]
  # Calculate perturbed differences
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  # 扰动
  beta1 <- sumwrfit(N,p,y_diff, X_diff, 10, r,'exp')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 10, r,'geom')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, r,'beta')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 10, r,'pois')
  mse4 <- MSE(beta_0, beta4)
  c(mse1,mse2,mse3,mse4)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE1", "MSE2","MSE3","MSE4","MSE5","MSE6")