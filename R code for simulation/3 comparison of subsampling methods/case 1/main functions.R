library(mvtnorm)
library(ggpubr)
library(quantreg)
library(lava)

###generate the random number
generate <- function(beta,N,Xtype,etype){
  p <- length(beta)
  mean <- rep(0,p)
  sigma <- matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i,j] <- 0.5^(abs(i-j)) 
    }
  }
  if(Xtype == "Normal"){
    x <- mvtnorm::rmvnorm(N,mean,sigma)
  }
  else if(Xtype == "t3"){
    x <- mvtnorm::rmvt(N,sigma,3)
  }
  
  if(etype == "Normal"){
    error <- rnorm(N) 
  }else if(etype == "Mixed"){
    u <- runif(N)
    error <- rep(0, N)
    for (i in 1:N){
      if(u[i]<0.9){
        error[i] <- rnorm(1,0,1)
      }else{
        error[i] <- rchisq(1, 5)
      }
    }
  }else if(etype=="t3"){
    error <- rt(N,4)
  }
  y <- x%*%beta+error
  return(list(y,x,error,sigma))
}


###Calculate the MSE
MSE <- function(x, y) {
  mse <- sum((x - y)^2)
  return(mse)
}


mynorm <- function(x){
  norma <- sqrt(x%*%x)
  return(norma)
}

##
norm_x <- function(x){
  nm <- apply(x, 1, mynorm)
  return(nm)
}

### optimal subsampling for linear model

linearsample <- function(y, x, prob) {
  yp <- diag(1/sqrt(prob)) %*% y
  xp <- diag(1/sqrt(prob)) %*% x
  beta <- lm(yp ~ . - 1, data = data.frame(xp, yp))$coefficients
  return(beta)
}

optsample <- function(y, x, r, method) {
  N <- nrow(x)
  p <- ncol(x) 
  
  if (method == "U") {
    prob <- rep(1/N, N)
  } else {
    r0 <- 100
    idx0 <- sample(1:N, r0, replace = TRUE, prob = rep(1/N, N))
    x0 <- x[idx0, ]
    y0 <- y[idx0]
    beta_0 <- lm(y0 ~ . - 1, data = data.frame(y0, x0))$coefficients
    e <- abs(y - x %*% beta_0)
    
    if (method == "A") {
      C <- crossprod(x) / (N * r)
      prob <- (e * norm_x(t(C %*% t(x)))) / sum(e * norm_x(t(C %*% t(x))))
    } else if (method == "L") {
      prob <- (e * norm_x(x)) / sum(e * norm_x(x))
    }
  }
  
  idx <- sample(1:N, r, replace = TRUE, prob = prob)
  y.opt <- y[idx]
  x.opt <- x[idx, ]
  prob.opt <- prob[idx]
  beta <- linearsample(y.opt, x.opt, prob.opt)
  return(beta)
}

### Weighted rank regression estimation 
##equivalent to LAD
Wrfit <- function(y, x, w) {
  n <- nrow(x)
  p <- ncol(x)
  
  ypairs <- pairup(y)
  yi <- ypairs[, 1]
  yj <- ypairs[, 2]
  xpairs <- pairup(x)
  xi <- xpairs[, 1:p]
  xj <- xpairs[, (p + 1):(2 * p)]
  wpairs <- pairup(w)
  wi <- wpairs[, 1]
  wj <- wpairs[, 2]
  bij <- wi * wj
  yw <- bij * (yi - yj)
  xw <- bij * (xi - xj)
  est <- quantreg::rq(yw ~ 0+xw, method = "fnb", tau = 0.5)$coefficients
  return(est)
}


persample <- function(y, x, u, v, modeltype) {
  N <- length(y)
  p <- ncol(x)
  num <- sum(u)
  idx.per <- which(u == 1)
  
  x.per <- x[idx.per, ]
  y.per <- y[idx.per]
  w.per <- v[idx.per]
  
  if (modeltype == "R") {
    beta.per <- Wrfit(y.per, x.per, w.per) 
  } else if (modeltype == "L") {
    beta.per <- linearsample(y.per, x.per, w.per)
  }
  
  return(beta.per)
}
sumwrfit <- function(n,p,y_diff, X_diff, m, r,method) {
  q <- r/n
  betas <- matrix(0, nrow = p, ncol = m)
  
  for (l in 1:m) {
    # Subset
    u <- rbinom(n, 1, q)
    
    # Stochastic weighting
    if(method == 'exp'){
      v <- rexp(n, q) 
    }else if(method == 'geom'){
      v <- rgeom(n, q)
    }else if(method == 'beta'){
      v <- runif(n,0,2/q)
    }else if(method == 'pois'){
      v <- rpois(n, 1/q)
    }
    
    # Weights
    W <- u * v
    W_diff <- compute_wij1(u, v)
    
    # Filter non-zero weights
    non_zero_idx <- which(W_diff != 0)
    y_per <- y_diff[non_zero_idx]
    X_per <- X_diff[non_zero_idx, ]
    W_per <- W_diff[non_zero_idx]
    
    # Weighted objective function
    yw <- W_per * y_per
    xw <- sweep(X_per, 1, W_per, "*")
    
    # Estimation
    rq_model <- quantreg::rq(yw ~ 0 + xw, method = "fnb", tau = 0.5)
    betas[, l] <- rq_model$coefficients
  }
  
  # Combination
  beta_est <- rowMeans(betas)
  
  return(beta_est)
}

mulwrfit <- function(n,p,y_diff, X_diff, m, r,method) {
  q <- r/n
  betas <- matrix(0, nrow = p, ncol = m)
  
  for (l in 1:m) {
    # Subset
    u <- rbinom(n, 1, q)
    
    # Stochastic weighting
    if(method == 'exp'){
      v <- rexp(n, q) 
    }else if(method == 'geom'){
      v <- rgeom(n, q)
    }else if(method == 'beta'){
      v <- runif(n,0,2/q)
    }else if(method == 'pois'){
      v <- rpois(n, 1/q)
    }
    
    # Weights
    W <- u * v
    W_diff <- compute_wij2(u, v)
    
    # Filter non-zero weights
    non_zero_idx <- which(W_diff != 0)
    y_per <- y_diff[non_zero_idx]
    X_per <- X_diff[non_zero_idx, ]
    W_per <- W_diff[non_zero_idx]
    
    # Weighted objective function
    yw <- W_per * y_per
    xw <- sweep(X_per, 1, W_per, "*")
    
    # Estimation
    rq_model <- quantreg::rq(yw ~ 0 + xw, method = "fnb", tau = 0.5)
    betas[, l] <- rq_model$coefficients
  }
  
  # Combination
  beta_est <- rowMeans(betas)
  
  return(beta_est)
}

mulwrfit2 <- function(n,p,y, x, m, r,method) {
  q <- r/n
  betas <- matrix(0, nrow = p, ncol = m)
  
  
  for (l in 1:m) {
    # Subset
    u <- rbinom(n, 1, q)
    # Stochastic weighting
    if(method == 'exp'){
      v <- rexp(n, q) 
    }else if(method == 'geom'){
      v <- rgeom(n, q)
    }else if(method == 'beta'){
      v <- runif(n,0,2/q)
    }else if(method == 'pois'){
      v <- rpois(n, 1/q)
    }
    num <- sum(u)
    idx.per <- which(u == 1)
    
    x.per <- x[idx.per, ]
    y.per <- y[idx.per]
    w.per <- v[idx.per]
    
    betas[, l]<- Wrfit(y.per, x.per, w.per) 
  }
  
  # Combination
  beta_est <- rowMeans(betas)
  
  return(beta_est)
}

# Helper function to compute differences
compute_yij <- function(y) {
  y <- as.matrix(y)
  N <- length(y)
  yij <- numeric(N * (N - 1) / 2)
  
  idx <- 1
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      yij[idx] <- y[i] - y[j]
      idx <- idx + 1
    }
  }
  
  return(yij)
}

compute_Xij <- function(X) {
  X <- as.matrix(X)
  N <- nrow(X)
  p <- ncol(X)
  
  combs <- combn(N, 2)
  num_combs <- ncol(combs)
  
  Xij <- matrix(0, nrow = num_combs, ncol = p)
  
  for (k in 1:num_combs) {
    i <- combs[1, k]
    j <- combs[2, k]
    Xij[k, ] <- X[i, ] - X[j, ]
  }
  
  return(Xij)
}

compute_wij1 <- function(u, v) {
  w <- u * v
  N <- length(w)
  wij <- numeric(N * (N - 1) / 2)
  
  idx <- 1
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      wij[idx] <- (w[i] + w[j]) / 2
      idx <- idx + 1
    }
  }
  return(wij)
}
compute_wij2 <- function(u, v) {
  w <- u * v
  N <- length(w)
  wij <- numeric(N * (N - 1) / 2)
  
  idx <- 1
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      wij[idx] <- w[i]*w[j]
      idx <- idx + 1
    }
  }
  return(wij)
}
compute_wijm <- function(u, v) {
  w <- u * v
  N <- length(w)
  wij <- numeric(N * (N - 1) / 2)
  
  idx <- 1
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      wij[idx] <- w[i]*w[j]
      idx <- idx + 1
    }
  }
  
  return(wij)
}

sumtime <- function(n, p, y_diff, X_diff, r) {
  q <- r / n
  starttime <- proc.time()
  
  # Subset
  u <- rbinom(n, 1, q)
  v <- rpois(n, 1 / q)
  
  # Weights
  W <- u * v
  W_diff <- compute_wij1(u, v)
  
  # Filter non-zero weights
  non_zero_idx <- which(W_diff != 0)
  y_per <- y_diff[non_zero_idx]
  X_per <- X_diff[non_zero_idx, ]
  W_per <- W_diff[non_zero_idx]
  
  # Weighted objective function
  yw <- W_per * y_per
  xw <- sweep(X_per, 1, W_per, "*")
  
  # Estimation
  betA <- quantreg::rq(yw ~ 0 + xw, method = "fnb", tau = 0.5)$coefficients
  
  endtime <- proc.time()
  timeM <- endtime[1] - starttime[1] 
  
  return(timeM)
}

multime <- function(n,p,y_diff, X_diff, r) {
  q <- r/n
  starttime <- proc.time()
  
  # Subset
  u <- rbinom(n, 1, q)
  v <- rpois(n, 1/q)
  
  # Weights
  W <- u * v
  W_diff <- compute_wij2(u, v)
  
  # Filter non-zero weights
  non_zero_idx <- which(W_diff != 0)
  y_per <- y_diff[non_zero_idx]
  X_per <- X_diff[non_zero_idx, ]
  W_per <- W_diff[non_zero_idx]
  
  # Weighted objective function
  yw <- W_per * y_per
  xw <- sweep(X_per, 1, W_per, "*")
  
  # Estimation
  betA <- quantreg::rq(yw ~ 0 + xw, method = "fnb", tau = 0.5)$coefficients
  
  endtime <- proc.time()
  timeM <- endtime[1] - starttime[1] 
  
  
  return(timeM)
}
