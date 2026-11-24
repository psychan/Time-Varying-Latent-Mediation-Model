gen_mediation_flat <- function(n=200,m=10){
  library(MASS)
  #n个人，m个时间点
  nm <- n*m
  #题目数
  i_x <- 3
  i_y <- 3
  i_m <- 3
  #从n(0,1)中生成x，观测和潜变量的方差都要是1
  t <- x <- nid <- um <- uy <- matrix(0,m,n)
  ui <- mvrnorm(n,numeric(2),diag(2)*0.1)
  for (i in 1:n) {
    t[,i] <- sort(runif(m,0,1))
    x[,i] <- as.matrix(rnorm(m,0,1))
    nid[,i] <- rep(i,m)
    um[,i] <- rep(0.4,m) + ui[i,1]
    uy[,i] <- rep(0.6,m) + ui[i,2]
  }
  time <- matrix(t,ncol=1)
  xx <- matrix(x,ncol=1)
  id <- matrix(nid,ncol=1)
  umi <- matrix(um,ncol=1)
  uyi <- matrix(uy,ncol=1)
  #x、m、y的载荷
  lam_x <- as.matrix(c(1,0.8,0.8))
  lam_y <- as.matrix(c(1,0.8,0.8))
  lam_m <- as.matrix(c(1,0.8,0.8))
  #各题目的残差
  sigma <- mvrnorm(nm,numeric(i_x+i_y+i_m),diag(i_x+i_y+i_m)*0.1)
  #x的项目
  item_x <- xx %*% t(lam_x) + sigma[,1:i_x]
  #回归各系数
  alpha <- 0.4
  beta <- 0.4
  gamma <- 0.4
  #回归的残差
  epsilon <- mvrnorm(nm,numeric(2),diag(2)*0.2)
  #生成m和y
  mm <- umi + alpha * xx + epsilon[,1]
  yy <- uyi + beta * mm + gamma * xx + epsilon[,2]
  #根据生成的m和y生成对应的项目
  item_m <- mm %*% t(lam_m) + sigma[,(i_x+1):(i_x+i_m)]
  item_y <- yy %*% t(lam_y) + sigma[,(i_x+i_m+1):(i_x+i_y+i_m)]
  #生成变量名
  x_name <- array('',dim = i_x)
  m_name <- array('',dim = i_m)
  y_name <- array('',dim = i_y)
  for (i in 1:i_x) {
    x_name[i] <- paste0('x',as.character(i))
  }
  for (i in 1:i_m) {
    m_name[i] <- paste0('m',as.character(i))
  }
  for (i in 1:i_y) {
    y_name[i] <- paste0('y',as.character(i))
  }
  data0 <- cbind(id,item_x,item_m,item_y,time)
  data_name <- c('id',x_name,m_name,y_name,'time')
  colnames(data0) <- data_name
  #把所有项目合并、存储
  #data0 <- list(x1=item_x[,1],x2=item_x[,2],x3=item_x[,3],
  #              m1=item_m[,1],m2=item_m[,2],m3=item_m[,3],
  #              y1=item_y[,1],y2=item_y[,2],y3=item_y[,3])
  #这个自创函数的返回值
  return(data0)
}

gen_mediation_valley <- function(n=200,m=10){
  library(MASS)
  #n个人，m个时间点
  nm <- n*m
  #题目数
  i_x <- 3
  i_y <- 3
  i_m <- 3
  #从n(0,1)中生成x，观测和潜变量的方差都要是1
  t <- x <- nid <- um <- uy <- matrix(0,m,n)
  ui <- mvrnorm(n,numeric(2),diag(2)*0.1)
  for (i in 1:n) {
    t[,i] <- sort(runif(m,0,1))
    x[,i] <- as.matrix(rnorm(m,0,1))
    nid[,i] <- rep(i,m)
    um[,i] <- rep(0.4,m) + ui[i,1]
    uy[,i] <- rep(0.6,m) + ui[i,2]
  }
  time <- matrix(t,ncol=1)
  xx <- matrix(x,ncol=1)
  id <- matrix(nid,ncol=1)
  umi <- matrix(um,ncol=1)
  uyi <- matrix(uy,ncol=1)
  #x、m、y的载荷
  lam_x <- as.matrix(c(1,0.8,0.8))
  lam_y <- as.matrix(c(1,0.8,0.8))
  lam_m <- as.matrix(c(1,0.8,0.8))
  #各题目的残差
  sigma <- mvrnorm(nm,numeric(i_x+i_y+i_m),diag(i_x+i_y+i_m)*0.1)
  #x的项目
  item_x <- xx %*% t(lam_x) + sigma[,1:i_x]
  #回归各系数
  alpha <- 1/(1+exp((-7)*(time-0.5)))-0.5
  beta <- 1/(1+exp((-7)*(time-0.5)))-0.5
  gamma <- 1/(1+exp((-7)*(time-0.5)))-0.5
  #回归的残差
  epsilon <- mvrnorm(nm,numeric(2),diag(2)*0.2)
  #生成m和y
  mm <- umi + alpha * xx + epsilon[,1]
  yy <- uyi + beta * mm + gamma * xx + epsilon[,2]
  #根据生成的m和y生成对应的项目
  item_m <- mm %*% t(lam_m) + sigma[,(i_x+1):(i_x+i_m)]
  item_y <- yy %*% t(lam_y) + sigma[,(i_x+i_m+1):(i_x+i_y+i_m)]
  #生成变量名
  x_name <- array('',dim = i_x)
  m_name <- array('',dim = i_m)
  y_name <- array('',dim = i_y)
  for (i in 1:i_x) {
    x_name[i] <- paste0('x',as.character(i))
  }
  for (i in 1:i_m) {
    m_name[i] <- paste0('m',as.character(i))
  }
  for (i in 1:i_y) {
    y_name[i] <- paste0('y',as.character(i))
  }
  data0 <- cbind(id,item_x,item_m,item_y,time)
  data_name <- c('id',x_name,m_name,y_name,'time')
  colnames(data0) <- data_name
  #把所有项目合并、存储
  #data0 <- list(x1=item_x[,1],x2=item_x[,2],x3=item_x[,3],
  #              m1=item_m[,1],m2=item_m[,2],m3=item_m[,3],
  #              y1=item_y[,1],y2=item_y[,2],y3=item_y[,3])
  #这个自创函数的返回值
  return(data0)
}

gen_mediation_hump <- function(n=200,m=10){
  library(MASS)
  #n个人，m个时间点
  nm <- n*m
  #题目数
  i_x <- 3
  i_y <- 3
  i_m <- 3
  #从n(0,1)中生成x，观测和潜变量的方差都要是1
  t <- x <- nid <- um <- uy <- matrix(0,m,n)
  ui <- mvrnorm(n,numeric(2),diag(2)*0.1)
  for (i in 1:n) {
    t[,i] <- sort(runif(m,0,1))
    x[,i] <- as.matrix(rnorm(m,0,1))
    nid[,i] <- rep(i,m)
    um[,i] <- rep(0.4,m) + ui[i,1]
    uy[,i] <- rep(0.6,m) + ui[i,2]
  }
  time <- matrix(t,ncol=1)
  xx <- matrix(x,ncol=1)
  id <- matrix(nid,ncol=1)
  umi <- matrix(um,ncol=1)
  uyi <- matrix(uy,ncol=1)
  #x、m、y的载荷
  lam_x <- as.matrix(c(1,0.8,0.8))
  lam_y <- as.matrix(c(1,0.8,0.8))
  lam_m <- as.matrix(c(1,0.8,0.8))
  #各题目的残差
  sigma <- mvrnorm(nm,numeric(i_x+i_y+i_m),diag(i_x+i_y+i_m)*0.1)
  #x的项目
  item_x <- xx %*% t(lam_x) + sigma[,1:i_x]
  #回归各系数
  alpha <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
  beta <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
  gamma <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
  #回归的残差
  epsilon <- mvrnorm(nm,numeric(2),diag(2)*0.2)
  #生成m和y
  mm <- umi + alpha * xx + epsilon[,1]
  yy <- uyi + beta * mm + gamma * xx + epsilon[,2]
  #根据生成的m和y生成对应的项目
  item_m <- mm %*% t(lam_m) + sigma[,(i_x+1):(i_x+i_m)]
  item_y <- yy %*% t(lam_y) + sigma[,(i_x+i_m+1):(i_x+i_y+i_m)]
  #生成变量名
  x_name <- array('',dim = i_x)
  m_name <- array('',dim = i_m)
  y_name <- array('',dim = i_y)
  for (i in 1:i_x) {
    x_name[i] <- paste0('x',as.character(i))
  }
  for (i in 1:i_m) {
    m_name[i] <- paste0('m',as.character(i))
  }
  for (i in 1:i_y) {
    y_name[i] <- paste0('y',as.character(i))
  }
  data0 <- cbind(id,item_x,item_m,item_y,time)
  data_name <- c('id',x_name,m_name,y_name,'time')
  colnames(data0) <- data_name
  #把所有项目合并、存储
  #data0 <- list(x1=item_x[,1],x2=item_x[,2],x3=item_x[,3],
  #              m1=item_m[,1],m2=item_m[,2],m3=item_m[,3],
  #              y1=item_y[,1],y2=item_y[,2],y3=item_y[,3])
  #这个自创函数的返回值
  return(data0)
}

gen_mediation_decay <- function(n=200,m=10){
  library(MASS)
  #n个人，m个时间点
  nm <- n*m
  #题目数
  i_x <- 3
  i_y <- 3
  i_m <- 3
  #从n(0,1)中生成x，观测和潜变量的方差都要是1
  t <- x <- nid <- um <- uy <- matrix(0,m,n)
  ui <- mvrnorm(n,numeric(2),diag(2)*0.1)
  for (i in 1:n) {
    t[,i] <- sort(runif(m,0,1))
    x[,i] <- as.matrix(rnorm(m,0,1))
    nid[,i] <- rep(i,m)
    um[,i] <- rep(0.4,m) + ui[i,1]
    uy[,i] <- rep(0.6,m) + ui[i,2]
  }
  time <- matrix(t,ncol=1)
  xx <- matrix(x,ncol=1)
  id <- matrix(nid,ncol=1)
  umi <- matrix(um,ncol=1)
  uyi <- matrix(uy,ncol=1)
  #x、m、y的载荷
  lam_x <- as.matrix(c(1,0.8,0.8))
  lam_y <- as.matrix(c(1,0.8,0.8))
  lam_m <- as.matrix(c(1,0.8,0.8))
  #各题目的残差
  sigma <- mvrnorm(nm,numeric(i_x+i_y+i_m),diag(i_x+i_y+i_m)*0.1)
  #x的项目
  item_x <- xx %*% t(lam_x) + sigma[,1:i_x]
  #回归各系数
  alpha <- -0.1+0.5*exp(-time^0.5)
  beta <- -0.1+0.5*exp(-time^0.5)
  gamma <- -0.1+0.5*exp(-time^0.5)
  #回归的残差
  epsilon <- mvrnorm(nm,numeric(2),diag(2)*0.2)
  #生成m和y
  mm <- umi + alpha * xx + epsilon[,1]
  yy <- uyi + beta * mm + gamma * xx + epsilon[,2]
  #根据生成的m和y生成对应的项目
  item_m <- mm %*% t(lam_m) + sigma[,(i_x+1):(i_x+i_m)]
  item_y <- yy %*% t(lam_y) + sigma[,(i_x+i_m+1):(i_x+i_y+i_m)]
  #生成变量名
  x_name <- array('',dim = i_x)
  m_name <- array('',dim = i_m)
  y_name <- array('',dim = i_y)
  for (i in 1:i_x) {
    x_name[i] <- paste0('x',as.character(i))
  }
  for (i in 1:i_m) {
    m_name[i] <- paste0('m',as.character(i))
  }
  for (i in 1:i_y) {
    y_name[i] <- paste0('y',as.character(i))
  }
  data0 <- cbind(id,item_x,item_m,item_y,time)
  data_name <- c('id',x_name,m_name,y_name,'time')
  colnames(data0) <- data_name
  #把所有项目合并、存储
  #data0 <- list(x1=item_x[,1],x2=item_x[,2],x3=item_x[,3],
  #              m1=item_m[,1],m2=item_m[,2],m3=item_m[,3],
  #              y1=item_y[,1],y2=item_y[,2],y3=item_y[,3])
  #这个自创函数的返回值
  return(data0)
}