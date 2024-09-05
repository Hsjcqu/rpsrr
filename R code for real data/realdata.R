##### Perturbation Subsample
rm(list=ls())
setwd("E:\\DSF\\R")

# 在并行环境中加载必要的包
library(mvtnorm)
library(doParallel)
library(foreach)
library(Rfit)

# 核心数
cores <- 56
cl <- makeCluster(cores)
registerDoParallel(cl)
#在并行环境中加载必要的函数
clusterEvalQ(cl, {
  source("function8.R")
})

data <- read.csv("E:/DSF/R/F.csv")
y <- as.matrix(data[,5])  
x <- as.matrix(data[,-5]) 
N <- nrow(x)
p <- ncol(x)

beta_0 <- rfit(y ~ x, data = data.frame(y, x))$coefficients[-1]
### 重复次数
H <- 200

# 执行并行计算
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  # 在线性回归中使用均匀分布和 L
  
  # Calculate perturbed differences
  y_diff <- compute_yij(y)
  X_diff <- compute_Xij(x)
  # 扰动
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

# 停止并行计算
stopCluster(cl)

# 计算平均 MSE 值
resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))
# 创建临时数据框
result_df <- data.frame(t(resultM))

# 将结果写入 CSV 文件
write.csv(result_df, file = "r=realdata.csv", row.names = FALSE)