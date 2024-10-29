rm(list=ls())


library(mvtnorm)
library(doParallel)
library(foreach)
library(Rfit)

cores <- 16
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("main functions.R")
})

data <- read.csv("F.csv")
y <- as.matrix(data[,5])  
x <- as.matrix(data[,-5]) 
N <- nrow(x)
p <- ncol(x)

beta_0 <- rfit(y ~ x, data = data.frame(y, x))$coefficients[-1]

H <- 500
###### product weights
###change the method in function mulwrfit2 can get the result for stochastic weights with different distributions
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  beta_U1 <- optsample(y, x, 100, 'U')
  beta_L1 <- optsample(y, x, 100, 'L')
  mse_U1 <- MSE(beta_0,beta_U1)
  mse_L1 <- MSE(beta_0,beta_L1)
  beta_U2 <- optsample(y, x, 200, 'U')
  beta_L2 <- optsample(y, x, 200, 'L')
  mse_U2 <- MSE(beta_0,beta_U2)
  mse_L2 <- MSE(beta_0,beta_L2)
  beta_U3 <- optsample(y, x, 300, 'U')
  beta_L3 <- optsample(y, x, 300, 'L')
  mse_U3 <- MSE(beta_0,beta_U3)
  mse_L3 <- MSE(beta_0,beta_L3)
  beta_U4 <- optsample(y, x, 400, 'U')
  beta_L4 <- optsample(y, x, 400, 'L')
  mse_U4 <- MSE(beta_0,beta_U4)
  mse_L4 <- MSE(beta_0,beta_L4)
  
  beta1 <- mulwrfit2(N,p,y,x,10,100,'pois')
  mse1 <- MSE(beta_0,beta1)
  beta2 <- mulwrfit2(N,p,y,x,10,200,'pois')
  mse2 <- MSE(beta_0,beta2)
  beta3 <- mulwrfit2(N,p,y,x,10,300,'pois')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- mulwrfit2(N,p,y,x,10,400,'pois')
  mse4 <- MSE(beta_0, beta4)
  c(mse_U1,mse_U2,mse_U3,mse_U4,mse_L1,mse_L2,mse_L3,mse_L4,mse1,mse2,mse3,mse4)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSEU1","MSEU2","MSEU3","MSEU4","MSEL1","MSEL2","MSEL3","MSEL4","MSE1","MSE2","MSE3","MSE4")


###### additive weights
###pois
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  beta1 <- sumwrfit(N,p,y_diff, X_diff, 10, 100,'pois')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 10, 200,'pois')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, 300,'pois')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 10, 400,'pois')
  mse4 <- MSE(beta_0, beta4)
  c(mse1,mse2,mse3,mse4)
}

stopCluster(cl)

# 计算平均 MSE 值
resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
# 创建临时数据框
result_df <- data.frame(t(resultM))

###geom
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)

  beta1 <- sumwrfit(N,p,y_diff, X_diff, 10, 100,'geom')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 10, 200,'geom')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, 300,'geom')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 10, 400,'geom')
  mse4 <- MSE(beta_0, beta4)
  c(mse1,mse2,mse3,mse4)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))

###exp
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  
  beta1 <- sumwrfit(N,p,y_diff, X_diff, 10, 100,'exp')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 10, 200,'exp')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, 300,'exp')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 10, 400,'exp')
  mse4 <- MSE(beta_0, beta4)
  c(mse1,mse2,mse3,mse4)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))

###beta
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  
  beta1 <- sumwrfit(N,p,y_diff, X_diff, 10, 100,'beta')
  mse1 <- MSE(beta_0, beta1)
  beta2 <- sumwrfit(N,p,y_diff, X_diff, 10, 200,'beta')
  mse2 <- MSE(beta_0, beta2)
  beta3 <- sumwrfit(N,p,y_diff, X_diff, 10, 300,'beta')
  mse3 <- MSE(beta_0, beta3)
  beta4 <- sumwrfit(N,p,y_diff, X_diff, 10, 400,'beta')
  mse4 <- MSE(beta_0, beta4)
  c(mse1,mse2,mse3,mse4)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
result_df <- data.frame(t(resultM))