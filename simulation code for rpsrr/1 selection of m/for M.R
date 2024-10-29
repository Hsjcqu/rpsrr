######this file is the code of selection of m.
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


beta_0 <- c(-4, 0.11, -0.9, -0.5, 1)
N <- 10000
p <- length(beta_0)

H <- 1000
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
write.csv(result_df, file = "r=50g.csv", row.names = FALSE)

##### r=50 / r = 100 
r <- 50
# r <- 100

##### product weights
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

write.csv(result_df, file = "r=50p.csv", row.names = FALSE)
#####

##### additive weights
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  beta1 <- sumwrfit(N,p,y_diff, X_diff, 1, r,'exp')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 5, r,'exp')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, r,'exp')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 15, r,'exp')
  mse4 <- MSE(beta_0, beta4)
  beta5 <- sumwrfit(N,p,y_diff, X_diff, 20, r,'exp')
  mse5 <- MSE(beta_0, beta5)
  beta6 <- sumwrfit(N,p,y_diff, X_diff, 25, r,'exp')
  mse6 <- MSE(beta_0, beta6)
  beta7 <- sumwrfit(N,p,y_diff, X_diff, 30, r,'exp')
  mse7 <- MSE(beta_0, beta7)
  
  c(mse1, mse2,mse3,mse4,mse5,mse6,mse7)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE_1", "MSE_5","MSE_10","MSE_15","MSE_20","MSE_25","MSE_30")

write.csv(result_df, file = "r=50a.csv", row.names = FALSE)
