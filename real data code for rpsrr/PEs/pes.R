rm(list=ls())
source("main functions.R")
# 在并行环境中加载必要的包
library(mvtnorm)
library(doParallel)
library(foreach)
library(Rfit)

# 核心数
cores <- 20
cl <- makeCluster(cores)
registerDoParallel(cl)


data <- read.csv("F.csv")
data <- scale(data,center=T,scale=T)
N <- nrow(data)
N.train <- round(N*0.75)
N.test <- N - N.train
traindata <- data[1:N.train,]
testdata <- data[(N.train+1):N,]
y <- as.matrix(traindata[,5])  
x <- as.matrix(traindata[,-5]) 

y.test <- as.matrix(testdata[,5])
x.test <- as.matrix(testdata[,-5])
p <- ncol(x)

H <- 500
###### product weights
###change the method in function mulwrfit2 can get the result for stochastic weights with different distributions
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  beta_U1 <- optsample(y, x, 100, 'U')
  beta_L1 <- optsample(y, x, 100, 'L')
  mse_U1 <- MSE(y.test, x.test%*%beta_U1)
  mse_L1 <- MSE(y.test, x.test%*%beta_L1)
  beta_U2 <- optsample(y, x, 200, 'U')
  beta_L2 <- optsample(y, x, 200, 'L')
  mse_U2 <- MSE(y.test, x.test%*%beta_U2)
  mse_L2 <- MSE(y.test, x.test%*%beta_L2)
  beta_U3 <- optsample(y, x, 300, 'U')
  beta_L3 <- optsample(y, x, 300, 'L')
  mse_U3 <- MSE(y.test, x.test%*%beta_U3)
  mse_L3 <- MSE(y.test, x.test%*%beta_L3)
  beta_U4 <- optsample(y, x, 400, 'U')
  beta_L4 <- optsample(y, x, 400, 'L')
  mse_U4 <- MSE(y.test, x.test%*%beta_U4)
  mse_L4 <- MSE(y.test, x.test%*%beta_L4)
  
  beta1 <- mulwrfit2(N.train,p,y,x,10,100,'pois')
  mse1 <- MSE(y.test, x.test%*%beta1)
  beta2 <- mulwrfit2(N.train,p,y,x,10,200,'pois')
  mse2 <- MSE(y.test, x.test%*%beta2)
  beta3 <- mulwrfit2(N.train,p,y,x,10,300,'pois')
  mse3 <- MSE(y.test, x.test%*%beta3)
  beta4 <- mulwrfit2(N.train,p,y,x,10,400,'pois')
  mse4 <- MSE(y.test, x.test%*%beta4)
  c(mse_U1,mse_U2,mse_U3,mse_U4,mse_L1,mse_L2,mse_L3,mse_L4,mse1,mse2,mse3,mse4)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSEU1","MSEU2","MSEU3","MSEU4","MSEL1","MSEL2","MSEL3","MSEL4","MSE1","MSE2","MSE3","MSE4")



#####additive weights
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  
  beta1 <- sumwrfit(N.train,p,y_diff, X_diff, 10, 100,'pois')
  mse1 <- MSE(y.test, x.test%*%beta1)
  beta2 <- sumwrfit(N.train,p,y_diff, X_diff, 10, 200,'pois')
  mse2 <- MSE(y.test, x.test%*%beta2)
  beta3 <- sumwrfit(N.train,p,y_diff, X_diff, 10, 300,'pois')
  mse3 <- MSE(y.test, x.test%*%beta3)
  beta4 <- sumwrfit(N.train,p,y_diff, X_diff, 10, 400,'pois')
  mse4 <- MSE(y.test, x.test%*%beta4)
  c(mse1,mse2,mse3,mse4)
}

# 停止并行计算
stopCluster(cl)

# 计算平均 MSE 值
resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
# 创建临时数据框
result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE1","MSE2","MSE3","MSE4")

# 将结果写入 CSV 文件
write.csv(result_df, file = "r=realdata12.csv", row.names = FALSE)