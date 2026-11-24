library(MASS)
library(invgamma)

wd <- 'F:/course/statistic/tvm/S6/s6a'
setwd(wd)

data0 <- read.csv(file='raw_ed_1.txt', header=TRUE, sep='\t')

# ---- 修改1：显式初始化为 NA 的 10000×1 矩阵 ----
sig_eta_x <- sig_eta_m <- sig_eta_y <- sig_eta_xm <- sig_eta_xy <-
  sig_eta_my <- matrix(NA, nrow=10000, ncol=1)

lambda_x   <- as.matrix(data0[,1:3])
psi_m      <- data0[,49]
sigma_umi  <- data0[,50]
lambda_m   <- as.matrix(data0[,4:6])
lambda_y   <- as.matrix(data0[,7:12])
ve_x       <- as.matrix(data0[,13:15])
ve_m       <- as.matrix(data0[,16:18])
ve_y       <- as.matrix(data0[,19:24])
um         <- as.matrix(data0[,25:36])
aa         <- as.matrix(data0[,37:48])
uy         <- as.matrix(data0[,51:62])
bb         <- as.matrix(data0[,63:74])
cc         <- as.matrix(data0[,75:86])
psi_y      <- data0[,87]
sigma_uyi  <- data0[,88]

data1 <- read.csv(file='data0.csv', header=TRUE)
data1 <- na.omit(data1)
data1$survey <- data1$survey/50

# 注意：只 scale 一次即可
data1[,6:17] <- scale(data1[,6:17])

nkot <- 9
n    <- nrow(data1)
t    <- data1[,2]
t0   <- rep(1,n); t1 <- t; t2 <- t^2
pre.k  <- seq(0,1,length.out=nkot+2)
myknot <- pre.k[2:(length(pre.k)-1)]

cf <- function(t, knot){
  mycf <- NULL
  for (i in 1:length(knot)){
    mycf[i] <- (ifelse(t >= knot[i], (t-knot[i])^2, 0))
  }
  return(mycf)
}

z_pi <- t(sapply(t, cf, knot=myknot))
z_t  <- cbind(t0, t1, t2)
eta_x <- eta_m <- eta_y <- matrix(0, n, 1)

xx <- as.matrix(data1[,6:8])
mm <- as.matrix(data1[,9:11])
yy <- as.matrix(data1[,12:17])

# ---- 修改2：固定随机种子，保证 mvrnorm 可重复 ----
set.seed(12345)

for (k in 1:10000) {
  theta_y <- diag(ve_y[k,])
  theta_x <- diag(ve_x[k,])
  theta_m <- diag(ve_m[k,])
  
  eta_x <- (xx - mvrnorm(n, rep(0,3), theta_x)) %*% t(ginv(lambda_x[k,]))
  eta_m <- (mm - mvrnorm(n, rep(0,3), theta_m)) %*% t(ginv(lambda_m[k,]))
  eta_y <- (yy - mvrnorm(n, rep(0,6), theta_y)) %*% t(ginv(lambda_y[k,]))
  
  sig_eta_x[k]  <- cov(eta_x)
  sig_eta_m[k]  <- cov(eta_m)
  sig_eta_y[k]  <- cov(eta_y)
  sig_eta_xm[k] <- cov(eta_x, eta_m)
  sig_eta_xy[k] <- cov(eta_x, eta_y)
  sig_eta_my[k] <- cov(eta_m, eta_y)
}

ll2 <- array(dim=c(10000,1))
SS    <- cov(cbind(xx, mm, yy))
# ---- 修改3：使用 ML 协方差 SS_ML，与 n*log|Sigma| 的似然更匹配 ----
SS_ML <- SS * (n-1) / n
p_obs <- ncol(SS_ML)   # 12

for (k in 1:10000) {
  theta_model <- matrix(c(sig_eta_x[k],  sig_eta_xm[k], sig_eta_xy[k],
                          sig_eta_xm[k], sig_eta_m[k],  sig_eta_my[k],
                          sig_eta_xy[k], sig_eta_my[k], sig_eta_y[k]),
                        nrow=3, ncol=3)
  
  lambda_model <- cbind(as.matrix(c(lambda_x[k,], rep(0,9))),
                        as.matrix(c(rep(0,3), lambda_m[k,], rep(0,6))),
                        as.matrix(c(rep(0,6), lambda_y[k,])))
  
  model_psi <- lambda_model %*% theta_model %*% t(lambda_model) +
    diag(c(ve_x[k,], ve_m[k,], ve_y[k,]))
  
  ll2[k] <- n * (log(det(model_psi)) +
                   sum(diag(SS_ML %*% solve(model_psi))) -
                   p_obs - log(det(SS_ML)))
}

# ---- 修改4：Dbar 与 Dhat 统一使用后 5000 个样本 ----
burn <- 5001:10000

# ---- 修改5：后验均值参数直接由后验样本均值构造，不再重新 mvrnorm 生成 eta ----
m_sig_x  <- mean(sig_eta_x[burn])
m_sig_m  <- mean(sig_eta_m[burn])
m_sig_y  <- mean(sig_eta_y[burn])
m_sig_xm <- mean(sig_eta_xm[burn])
m_sig_xy <- mean(sig_eta_xy[burn])
m_sig_my <- mean(sig_eta_my[burn])

theta_model_m <- matrix(c(m_sig_x,  m_sig_xm, m_sig_xy,
                          m_sig_xm, m_sig_m,  m_sig_my,
                          m_sig_xy, m_sig_my, m_sig_y),
                        nrow=3, ncol=3)

m_lambda_x <- colMeans(lambda_x[burn,])
m_lambda_m <- colMeans(lambda_m[burn,])
m_lambda_y <- colMeans(lambda_y[burn,])
m_ve_x     <- colMeans(ve_x[burn,])
m_ve_m     <- colMeans(ve_m[burn,])
m_ve_y     <- colMeans(ve_y[burn,])

lambda_model_m <- cbind(as.matrix(c(m_lambda_x, rep(0,9))),
                        as.matrix(c(rep(0,3), m_lambda_m, rep(0,6))),
                        as.matrix(c(rep(0,6), m_lambda_y)))

model_psi_m <- lambda_model_m %*% theta_model_m %*% t(lambda_model_m) +
  diag(c(m_ve_x, m_ve_m, m_ve_y))

ll2_m <- n * (log(det(model_psi_m)) +
                sum(diag(SS_ML %*% solve(model_psi_m))) -
                p_obs - log(det(SS_ML)))

# ---- 修改6：DIC / pD 标准公式 ----
Dbar <- mean(ll2[burn])
Dhat <- ll2_m

p_d <- Dbar - Dhat
DIC <- Dbar + p_d

p_d
DIC
