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

####### r \in {50,100,150,200,250,300}
r <- 300 #change r to get the results

##### product weights
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  
  M <- 10
  betaper_exp <- betaper_geom <- betaper_beta <- betaper_pois <- matrix(0, M, p)
  
  for (l in 1:M) {
    q <- r / N
    u <- rbinom(N, 1, q) 
    v_exp <- rexp(N, q) 
    v_beta <- runif(N,0,2/q)
    v_geom <- rgeom(N, q)
    v_pois <- rpois(N, 1/q)
    betaper_exp[l,] <- persample(y, x, u, v_exp,"R") 
    betaper_geom[l,] <- persample(y, x, u, v_geom,"R") 
    betaper_beta[l,] <- persample(y, x, u, v_beta,"R") 
    betaper_pois[l,] <- persample(y, x, u, v_pois,"R") 
  }
  
  beta_exp <- colMeans(betaper_exp)
  mse_exp <- MSE(beta_0, beta_exp)
  beta_geom <- colMeans(betaper_geom)
  mse_geom <- MSE(beta_0, beta_geom)
  beta_beta <- colMeans(betaper_beta)
  mse_beta <- MSE(beta_0, beta_beta)
  beta_pois <- colMeans(betaper_pois)
  mse_pois <- MSE(beta_0, beta_pois)
  
  c(mse_exp, mse_geom,mse_beta,mse_pois)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_dfp <- data.frame(t(resultM))
colnames(result_dfp) <- c("MSE_Exp", "MSE_Geom","MSE_beta","MSE_pois")

##### additive weights
result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "Normal", "Normal") 
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

result_dfa <- data.frame(t(resultM))
colnames(result_dfa) <- c("MSE1", "MSE2","MSE3","MSE4")
