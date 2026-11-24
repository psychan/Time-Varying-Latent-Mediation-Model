#m=u_m(t)+alpha(t)*x+epsilon_m, epsilon_m~N(0,psi_m)
#y=u_y(t)+beta(t)*m+gamma(t)*x+epsilon_y, epsilon_y~N(0,psi_y)
##X=v_x+lambda_x*x+ve_x, ve_x~N(0,theta_x)
##M=v_m+lambda_m*m+ve_m, ve_m~N(0,theta_m)
##Y=v_y+lambda_y*y+ve_y, ve_y~N(0,theta_y)
###alpha(t)=phi_a0+phi_a1*t+phi_a2*t^2+ak*(t-pik)^2_+, ak~N(0,zeta_a)
###beta(t)=phi_b0+phi_b1*t+phi_b2*t^2+bk*(t-pik)^2_+, bk~N(0,zeta_b)
###gamma(t)=phi_c0+phi_c1*t+phi_c2*t^2+ck*(t-pik)^2_+, ck~N(0,zeta_c)
###u_m(t)=phi_d0+phi_d1*t+phi_d2*t^2+dk*(t-pik)^2_+ +umi, dk~N(0,zeta_d) umi~(0,sigma_um)
###u_y(t)=phi_e0+phi_e1*t+phi_e2*t^2+ek*(t-pik)^2_+ +uyi, ek~N(0,zeta_e) uyi~(0,sigma_uy)
tvmediation <- function(data=NULL,x=NULL,m=NULL,y=NULL,cluster=NULL,t=Null,
                        nkot=9,iter=10000,burnin=5000){
  library(MASS)
  library(invgamma)
  library(LaplacesDemon)
  
  timepoint <- (1:nkot)*0.1
  ind.est <- alpha.est <- beta.est <- gamma.est <- um.est <- uy.est <- 
    matrix(0,iter,length(timepoint));
  
  data <- as.data.frame(data)
  n <- nrow(data) #全部观测值个数
  nid <- unique(cluster)
  n_sample <- length(nid) #样本量
  sub <- as.numeric(match(cluster, nid))
  n_obs <- as.numeric(as.vector(table(cluster)))
  
  t0 <- rep(1,n); t1 <- t; t2 <- t^2;
  pre.k <- seq(0,1,length.out=nkot+2) #生成节点序列
  myknot <- pre.k[2:(length(pre.k)-1)] #去掉头尾的节点序列
  cf<-function(t,knot){
    #截断指数基的截断值
    mycf <- NULL
    for (i in 1:length(knot)){
      mycf[i]<-(ifelse(t>=knot[i],(t-knot[i])^2,0))
    }
    return(mycf)
  }
  z_pi <- t(sapply(t,cf,knot=myknot)) # n*9
  z_t <- cbind(t0,t1,t2) # n*3的矩阵
  
  my.hat <- function(t,alpha0,alpha1,alpha2,a1,a2,a3,a4,a5,a6,a7,a8,a9){
    #拟合函数
    f1 <- (ifelse(t>=0.1,(t-0.1)^2,0))
    f2 <- (ifelse(t>=0.2,(t-0.2)^2,0))
    f3 <- (ifelse(t>=0.3,(t-0.3)^2,0))
    f4 <- (ifelse(t>=0.4,(t-0.4)^2,0))
    f5 <- (ifelse(t>=0.5,(t-0.5)^2,0))
    f6 <- (ifelse(t>=0.6,(t-0.6)^2,0))
    f7 <- (ifelse(t>=0.7,(t-0.7)^2,0))
    f8 <- (ifelse(t>=0.8,(t-0.8)^2,0))
    f9 <- (ifelse(t>=0.9,(t-0.9)^2,0))
    rr <- alpha0+alpha1*t+alpha2*t^2+a1*f1+a2*f2+a3*f3+a4*f4+a5*f5+a6*f6+a7*f7+a8*f8+a9*f9;
    return(rr)
  }
  
  xx <- as.matrix(x); mm <- as.matrix(m); yy <- as.matrix(y);
  n_xitem <- ncol(xx)
  n_mitem <- ncol(mm)
  n_yitem <- ncol(yy)
  
  #需要估计的参数及初始化
  ##载荷
  lambda_x <- matrix(0,iter+1,n_xitem) #x的载荷
  lambda_m <- matrix(0,iter+1,n_mitem) #m的载荷
  lambda_y <- matrix(0,iter+1,n_yitem) #y的载荷
  lambda_x[,1] <- lambda_m[,1] <- lambda_y[,1] <- 1
  lambda_x[1,2:n_xitem] <- rnorm(n_xitem-1,0.5,0.1)
  lambda_m[1,2:n_mitem] <- rnorm(n_mitem-1,0.5,0.1)
  lambda_y[1,2:n_yitem] <- rnorm(n_yitem-1,0.5,0.1)
  ##项目残差
  ve_x <- matrix(0,iter+1,n_xitem)
  ve_m <- matrix(0,iter+1,n_mitem)
  ve_y <- matrix(0,iter+1,n_yitem)
  ve_x[1,] <- diag(cov(xx))
  ve_m[1,] <- diag(cov(mm))
  ve_y[1,] <- diag(cov(yy))
  theta_x <- diag(ve_x[1,])
  theta_m <- diag(ve_m[1,])
  theta_y <- diag(ve_y[1,])
  ##项目截距
  v_x <- matrix(0,iter+1,n_xitem)
  v_m <- matrix(0,iter+1,n_mitem)
  v_y <- matrix(0,iter+1,n_yitem)
  ##回归系数
  alpha <- beta <- gamma <- matrix(0,iter+1,3)
  alpha[1,] <- beta[1,] <- gamma[1,] <- rnorm(1,1,0.1)
  ##回归截距
  u_m <- u_y <- matrix(0,iter+1,3)
  ###截距残差
  sigma_umi <- sigma_uyi <- matrix(0,iter+1,1)
  sigma_umi[1,] <- sigma_uyi[1,] <- rnorm(1,1,0.1)
  u_m[1,1] <- rnorm(1,0,sigma_umi[1,])
  u_y[1,1] <- rnorm(1,0,sigma_uyi[1,])
  ##回归残差
  psi_m <- psi_y <- matrix(0,iter+1,1)
  psi_m[1,] <- psi_y[1,] <- rnorm(1,1,0.1)
  ##截断基
  ak <- bk <- ck <- dk <- ek <- matrix(0,iter+1,9)
  zeta_a <- zeta_b <- zeta_c <- zeta_d <- zeta_e <- matrix(0,iter+1,1)
  zeta_a[1,] <- zeta_b[1,] <- zeta_c[1,] <- zeta_d[1,] <- zeta_e[1,] <- 
    rnorm(1,1,0.1); #相当于截断基的先验
  ak[1,] <- rnorm(9,0,zeta_a[1,])
  bk[1,] <- rnorm(9,0,zeta_b[1,])
  ck[1,] <- rnorm(9,0,zeta_c[1,])
  dk[1,] <- rnorm(9,0,zeta_d[1,])
  ek[1,] <- rnorm(9,0,zeta_e[1,])
  ###整合成大矩阵
  abc1 <- cbind(u_m,dk,alpha,ak) #m的回归(3+9+3+9)*iter
  abc2 <- cbind(u_y,ek,beta,bk,gamma,ck) #y的回归(3+9+3+9+3+9)*iter
  u1 <- u2 <- matrix(0,n_sample,iter) #每个被试自身的均值，被试数*iter
  
  #先验
  ##截距斜率载荷
  mu_lx <- mu_vx <- rep(0,n_xitem)
  mu_lm <- mu_vm <- rep(0,n_mitem)
  mu_ly <- mu_vy <- rep(0,n_yitem)
  mu_alpha <- mu_beta <- mu_gamma <- rep(0,3)
  mu_um <- mu_uy <- rep(0,3)
  sigma_lx <- sigma_vx <- diag(10^4,n_xitem,n_xitem)
  sigma_lm <- sigma_vm <- diag(10^4,n_mitem,n_mitem)
  sigma_ly <- sigma_vy <- diag(10^4,n_yitem,n_yitem)
  sigma_alpha <- sigma_beta <- sigma_gamma <- diag(10^4,3,3)
  sigma_um <- sigma_uy <- diag(10^4,3,3)
  ##方差
  a_em <- a_ey <- 0.001
  a_sa <- a_sb <- a_sc <- a_sd <- a_se <- 0.001
  a_um <- a_uy <- 0.001
  b_em <- b_ey <- 0.001 
  b_sa <- b_sb <- b_sc <- b_sd <- b_se <- 0.001
  b_um <- b_uy <- 0.001
  
  eta_x <- eta_m <- eta_y <- matrix(0,n,1)
  sig_eta_x <- 10^4
  
  for (k in 1:iter) {
    #B1:eta_x lambda_x v_x
    a_y <- ginv(t(lambda_y[k,])%*%ginv(theta_y)%*%lambda_y[k,]+
                  ginv(psi_y[k,]+sigma_uyi[k,]));
    for (i in 1:n) {
      at <- abc1[k,(3+1+nkot):(3+nkot+3)]%*%z_t[i,]+
        abc1[k,(3+1+nkot+3):(2*(3+nkot))]%*%z_pi[i,];
      bt <- abc2[k,(3+1+nkot):(3+nkot+3)]%*%z_t[i,]+
        abc2[k,(3+1+nkot+3):(2*(nkot+3))]%*%z_pi[i,];
      gt <- abc2[k,(2*(nkot+3)+1):(2*(nkot+3)+3)]%*%z_t[i,]+
        abc2[k,(2*(nkot+3)+1+3):(3*(nkot+3))]%*%z_pi[i,];
      
      m0 <- abc1[k,1]+abc1[k,2:3]%*%z_t[i,2:3]+abc1[k,(3+1):(3+nkot)]%*%z_pi[i,];
      y0 <- abc2[k,1]+abc2[k,2:3]%*%z_t[i,2:3]+abc2[k,(3+1):(3+nkot)]%*%z_pi[i,];
      
      psi_mx <- (at^(-1))*(psi_m[k,]+sigma_umi[k,])*(at^(-1))
      psi_yx <- (gt^(-1))*(psi_y[k,]+sigma_uyi[k,])*(gt^(-1))
      psi_ym <- (bt^(-1))*(psi_y[k,]+sigma_uyi[k,])*(bt^(-1))
      #eta_x
      a_x <- ginv(t(lambda_x[k,])%*%ginv(theta_x)%*%lambda_x[k,]+ginv(psi_mx)+
                    ginv(psi_yx)+sig_eta_x^(-1));
      b_x <- t(lambda_x[k,])%*%ginv(theta_x)%*%(xx[i,]-v_x[k,])+
        ginv(psi_mx)*(eta_m[i]-m0)/at+
        ginv(psi_yx)*(eta_y[i]-y0-bt*eta_m[i])/gt;
      eta_x[i] <- rnorm(1,a_x*b_x,sqrt(a_x));
      
      #eta_m
      a_m <- ginv(t(lambda_m[k,])%*%ginv(theta_m)%*%lambda_m[k,]+
                    ginv(psi_m[k,]+sigma_umi[k,])+ginv(psi_ym));
      b_m <- t(lambda_m[k,])%*%ginv(theta_m)%*%(mm[i,]-v_m[k,])+
        ginv(psi_m[k,]+sigma_umi[k,])*(m0+at*eta_x[i])+
        ginv(psi_ym)*(eta_y[i]-y0-gt*eta_x[i])/bt;
      eta_m[i] <- rnorm(1,a_m*b_m,sqrt(a_m));
      
      #eta_y
      b_y <- t(lambda_y[k,])%*%ginv(theta_y)%*%(yy[i,]-v_y[k,])+
        ginv(psi_y[k,]+sigma_uyi[k,])*(y0+bt*eta_m[i]+gt*eta_x[i]);
      eta_y[i] <- rnorm(1,a_y*b_y,sqrt(a_y));
    }
    
    sig_eta_x <- rinvgamma(1,0.001+n/2,0.001+sum((eta_x)^2)/2);
    
    a_lx <- ginv(as.numeric(t(eta_x)%*%diag(rep(1,n))%*%eta_x)*ginv(theta_x)+
                   ginv(sigma_lx));
    b_lx <- ginv(theta_x)%*%t(t(eta_x)%*%diag(rep(1,n))%*%(xx-v_x[k,]))+
      ginv(sigma_lx)%*%mu_lx;
    lamx <- mvrnorm(1,a_lx%*%b_lx,a_lx)
    lambda_x[k+1,2:n_xitem] <- lamx[2:n_xitem]
    
    #a_vx <- ginv(n*ginv(theta_x)+ginv(sigma_vx))
    #b_vx <- ginv(theta_x)%*%colSums(xx-eta_x%*%lambda_x[k+1,])+ginv(sigma_vx)%*%mu_vx
    #v_x[k+1,] <- mvrnorm(1,a_vx%*%b_vx,a_vx)
    v_x[k+1,] <- 0
    
    a_lm <- ginv(as.numeric(t(eta_m)%*%diag(rep(1,n))%*%eta_m)*ginv(theta_m)+
                   ginv(sigma_lm));
    b_lm <- ginv(theta_m)%*%t(t(eta_m)%*%diag(rep(1,n))%*%(mm-v_m[k,]))+
      ginv(sigma_lm)%*%mu_lm;
    lamm <- mvrnorm(1,a_lm%*%b_lm,a_lm)
    lambda_m[k+1,2:n_mitem] <- lamm[2:n_mitem]
    
    #a_vm <- ginv(n*ginv(theta_m)+ginv(sigma_vm))
    #b_vm <- ginv(theta_m)%*%colSums(mm-eta_m%*%lambda_m[k+1,])+ginv(sigma_vm)%*%mu_vm
    #v_m[k+1,] <- mvrnorm(1,a_vm%*%b_vm,a_vm)
    v_m[k+1,] <- 0
    
    a_ly <- ginv(as.numeric(t(eta_y)%*%diag(rep(1,n))%*%eta_y)*ginv(theta_y)+
                   ginv(sigma_ly))
    b_ly <- ginv(theta_y)%*%t(t(eta_y)%*%diag(rep(1,n))%*%(yy-v_y[k,]))+
      ginv(sigma_ly)%*%mu_ly
    lamy <- mvrnorm(1,a_ly%*%b_ly,a_ly)
    lambda_y[k+1,2:n_yitem] <- lamy[2:n_yitem]
    
    #a_vy <- ginv(n*ginv(theta_y)+ginv(sigma_vy))
    #b_vy <- ginv(theta_y)%*%colSums(yy-eta_y%*%lambda_y[k+1,])+ginv(sigma_vy)%*%mu_vy
    #v_y[k+1,] <- mvrnorm(1,a_vy%*%b_vy,a_vy)
    v_y[k+1,] <- 0
    
    #B4:m的回归
    wx1 <- cbind(z_t[,2:3],z_pi,z_t*matrix(rep(eta_x,3),ncol=3),
                 z_pi*matrix(rep(eta_x,nkot),ncol=nkot));
    wx2 <- cbind(z_t[,2:3],z_pi,z_t*matrix(rep(eta_m,3),ncol=3),
                 z_pi*matrix(rep(eta_m,nkot),ncol=nkot),
                 z_t*matrix(rep(eta_x,3),ncol=3),
                 z_pi*matrix(rep(eta_x,nkot),ncol=nkot));
    for (i in 1:n_sample) {
      sig_u1 <- (n_obs[i]*psi_m[k,]^(-1)+sigma_umi[k,]^(-1))^(-1) #m的回归
      sig_u2 <- (n_obs[i]*psi_y[k,]^(-1)+sigma_uyi[k,]^(-1))^(-1)
      mu_u1 <- sig_u1*(psi_m[k,]^(-1)*sum(eta_m[sub==i]-wx1[sub==i,]%*%abc1[k,2:(2*(3+nkot))])+
                         sigma_umi[k,]^(-1)*abc1[k,1]);
      mu_u2 <- sig_u2*(psi_y[k,]^(-1)*sum(eta_y[sub==i]-wx2[sub==i,]%*%abc2[k,2:(3*(3+nkot))])+
                         sigma_uyi[k,]^(-1)*abc2[k,1]);
      u1[i,k] <- rnorm(1,mu_u1,sqrt(sig_u1)) #每个被试一个截距u
      u2[i,k] <- rnorm(1,mu_u2,sqrt(sig_u2))
    }
    
    sig_a0 <- (n_sample*sigma_umi[k,]^(-1)+10^(-4))^(-1)
    sig_b0 <- (n_sample*sigma_uyi[k,]^(-1)+10^(-4))^(-1)
    mu_a0 <- sig_a0*sigma_umi[k,]^(-1)*sum(u1[,k])
    mu_b0 <- sig_b0*sigma_uyi[k,]^(-1)*sum(u2[,k])
    
    abc1[k+1,1]= rnorm(1,mu_a0,sqrt(sig_a0)) # 每个被试的截距=截距均值+残差
    abc2[k+1,1]= rnorm(1,mu_b0,sqrt(sig_b0))
    
    sig_m <- diag(c(rep(1,2)*10^4,rep(1,nkot)*zeta_d[k],
                    rep(1,3)*10^4,rep(1,nkot)*zeta_a[k])); #方程M的残差方差的先验
    temp_mu_a <- array(dim=c(2*(3+nkot)-1,n_sample)); # 23（截距和x的值）*被试数，存均值
    temp_v_a <- array(dim=c(2*(3+nkot)-1,2*(3+nkot)-1,n_sample));  # 23（截距和x的值）*23（截距和x的值）*被试数，存方差协方差矩阵
    
    for(i in 1:n_sample){ 
      if(n_obs[i]==1){ #如果只测量一次
        temp_v_a[,,i]=(wx1[sub==i,])%*%t(wx1[sub==i,])*psi_m[k]^(-1);  
        temp_mu_a[,i]=(wx1[sub==i,])%*%t(eta_m[sub==i]-u1[i,k])*psi_m[k]^(-1);    
      }else{
        temp_v_a[,,i]=t(wx1[sub==i,])%*%diag(rep(1,n_obs[i]))%*%(wx1[sub==i,])*psi_m[k]^(-1);  
        temp_mu_a[,i]=t(wx1[sub==i,])%*%diag(rep(1,n_obs[i]))%*%(eta_m[sub==i]-u1[i,k])*psi_m[k]^(-1);
      }
    }
    
    sig_abc1 <- ginv(apply(temp_v_a,c(1,2),sum)+ginv(sig_m))
    mu_abc1 <- sig_abc1%*%(apply(temp_mu_a,c(1),sum))
    
    abc1[k+1,2:(2*(3+nkot))] <- mvrnorm(1,mu_abc1,sig_abc1)
    
    #B5：y的回归
    sig_y <- diag(c(rep(1,2)*10^4,rep(1,nkot)*zeta_e[k],rep(1,3)*10^4,
                    rep(1,nkot)*zeta_b[k],rep(1,3)*10^4,rep(1,nkot)*zeta_c[k]));#方程y的残差方差的先验
    
    temp_mu_b <- array(dim=c(3*(3+nkot)-1,n_sample));
    temp_v_b <- array(dim=c(3*(3+nkot)-1,3*(3+nkot)-1,n_sample));
    
    for(i in 1:n_sample){ 
      if(n_obs[i]==1){
        temp_v_b[,,i]=(wx2[sub==i,])%*%t(wx2[sub==i,])*psi_y[k]^(-1);  
        temp_mu_b[,i]=(wx2[sub==i,])%*%t(eta_y[sub==i]-u2[i,k])*psi_y[k]^(-1);
      }else{
        temp_v_b[,,i]=t(wx2[sub==i,])%*%diag(rep(1,n_obs[i]))%*%wx2[sub==i,]*psi_y[k]^(-1);  
        temp_mu_b[,i]=t(wx2[sub==i,])%*%diag(rep(1,n_obs[i]))%*%(eta_y[sub==i]-u2[i,k])*psi_y[k]^(-1);    
      }
    }
    sig_abc2 <- ginv(apply(temp_v_b,c(1,2),sum)+ginv(sig_y))
    mu_abc2 <- sig_abc2%*%(apply(temp_mu_b,c(1),sum))
    
    abc2[k+1,2:(3*(3+nkot))] <- mvrnorm(1,mu_abc2,sig_abc2)
    
    #B6
    for (j in 1:n_xitem) {
      ve_x[k+1,j] <- rinvgamma(1,0.001+n/2,0.001+sum((xx[,j]-v_x[k+1,j]-eta_x%*%lambda_x[k+1,j])^2)/2);
    }
    theta_x <- diag(ve_x[k+1,])
    for (j in 1:n_mitem) {
      ve_m[k+1,j] <- rinvgamma(1,0.001+n/2,0.001+sum((mm[,j]-v_m[k+1,j]-eta_m%*%lambda_m[k+1,j])^2)/2);
    }
    theta_m <- diag(ve_m[k+1,])
    for (j in 1:n_yitem) {
      ve_y[k+1,j] <- rinvgamma(1,0.001+n/2,0.001+sum((yy[,j]-v_y[k+1,j]-eta_y%*%lambda_y[k+1,j])^2)/2);
    }
    theta_y <- diag(ve_y[k+1,]);
    
    zeta_a[k+1,] <- rinvgamma(1,a_sa+nkot/2,b_sa+
                                sum(abc1[k+1,(3+nkot+3+1):(2*(3+nkot))]^2)/2);
    zeta_b[k+1,] <- rinvgamma(1,a_sb+nkot/2,b_sb+
                                sum(abc2[k+1,(3+nkot+3+1):(2*(3+nkot))]^2)/2);
    zeta_c[k+1,] <- rinvgamma(1,a_sc+nkot/2,b_sc+
                                sum(abc2[k+1,(2*(3+nkot)+3+1):(3*(3+nkot))]^2)/2);
    zeta_d[k+1,] <- rinvgamma(1,a_sd+nkot/2,b_sd+
                                sum(abc1[k+1,(3+1):(3+nkot)]^2)/2);
    zeta_e[k+1,] <- rinvgamma(1,a_se+nkot/2,b_se+
                                sum(abc2[k+1,(3+1):(3+nkot)]^2)/2);
    
    tem_E <- tem_F <- matrix(0,n_sample,1);
    for(i in 1:n_sample){
      tem_E[i] <- sum((eta_m[sub==i]-u1[i,k]-wx1[sub==i,]%*%abc1[k+1,2:(2*(3+nkot))])^2); #m回归的残差
      tem_F[i] <- sum((eta_y[sub==i]-u2[i,k]-wx2[sub==i,]%*%abc2[k+1,2:(3*(3+nkot))])^2); #y回归的残差
    }
    psi_m[k+1] <- rinvgamma(1,a_em+n/2,b_em+sum(tem_E)/2); #m回归的残差
    psi_y[k+1] <- rinvgamma(1,a_ey+n/2,b_ey+sum(tem_F)/2); #y回归的残差
    
    sigma_umi[k+1,] <- rinvgamma(1,a_um+n_sample/2,b_um+sum((u1[,k]-abc1[k+1,1])^2)/2); #截距的残差
    sigma_uyi[k+1,] <- rinvgamma(1,a_uy+n_sample/2,b_uy+sum((u2[,k]-abc2[k+1,1])^2)/2); #截距的残差
  }
  for (k in 1:iter) {
    #计算估计值
    aa <- abc1[k,(3+nkot+1):(2*(3+nkot))]; # 3+9 x对m的回归
    bb <- abc2[k,(3+nkot+1):(2*(3+nkot))]; # 3+9 m对y的回归
    cc <- abc2[k,(2*(3+nkot)+1):(3*(3+nkot))]; # 3+9 x对y的回归
    dd <- abc1[k,1:(3+nkot)]; # m的截距
    ee <- abc2[k,1:(3+nkot)]; # y的截距
    for(j in 1:length(timepoint)){
      alpha.est[k,j] <- my.hat(timepoint[j], aa[1],aa[2],aa[3],aa[4],aa[5],aa[6],aa[7],aa[8],aa[9],aa[10],aa[11],aa[12]);
      beta.est[k,j] <- my.hat(timepoint[j], bb[1],bb[2],bb[3],bb[4],bb[5],bb[6],bb[7],bb[8],bb[9],bb[10],bb[11],bb[12]);
      gamma.est[k,j] <- my.hat(timepoint[j], cc[1],cc[2],cc[3],cc[4],cc[5],cc[6],cc[7],cc[8],cc[9],cc[10],cc[11],cc[12]);
      ind.est[k,j] <- alpha.est[k,j]*beta.est[k,j];
      um.est[k,j] <- my.hat(timepoint[j], dd[1],dd[2],dd[3],dd[4],dd[5],dd[6],dd[7],dd[8],dd[9],dd[10],dd[11],dd[12]);
      uy.est[k,j] <- my.hat(timepoint[j], ee[1],ee[2],ee[3],ee[4],ee[5],ee[6],ee[7],ee[8],ee[9],ee[10],ee[11],ee[12]);
    }
  }
  outputResults <- function(para){
    if(is.array(para)==TRUE){
      output <- rbind(colMeans(para),apply(para,2,sd,na.rm = TRUE),
                      apply(para,2,quantile,probs = c(0.025),na.rm = TRUE),
                      apply(para,2,quantile,probs = c(0.975),na.rm = TRUE));
    }else{
      output <- rbind(mean(para),sd(para),quantile(para,0.025),quantile(para,0.975));
    }
    rownames(output) <- c("mean", "sd", "95% lower CI", "95% upper CI");
    return(output)
  }
  
  lam_x <- outputResults(lambda_x[(burnin+2):(iter+1),])
  lam_m <- outputResults(lambda_m[(burnin+2):(iter+1),])
  lam_y <- outputResults(lambda_y[(burnin+2):(iter+1),])
  #vx <- outputResults(v_x[(burnin+2):(iter+1),])
  #vm <- outputResults(v_m[(burnin+2):(iter+1),])
  #vy <- outputResults(v_y[(burnin+2):(iter+1),])
  vex <- outputResults(ve_x[(burnin+2):(iter+1),])
  vem <- outputResults(ve_m[(burnin+2):(iter+1),])
  vey <- outputResults(ve_y[(burnin+2):(iter+1),])
  um <- outputResults(um.est[(burnin+1):iter,])
  aa <- outputResults(alpha.est[(burnin+1):iter,])
  uy <- outputResults(uy.est[(burnin+1):iter,])
  bb <- outputResults(beta.est[(burnin+1):iter,])
  cc <- outputResults(gamma.est[(burnin+1):iter,])
  ind <- outputResults(ind.est[(burnin+1):iter,])
  psim <- outputResults(psi_m[(burnin+2):(iter+1),])
  psiy <- outputResults(psi_y[(burnin+2):(iter+1),])
  sigmaumi <- outputResults(sigma_umi[(burnin+2):(iter+1),])
  sigmauyi <- outputResults(sigma_uyi[(burnin+2):(iter+1),])
  
  results.est <- cbind(lam_x,lam_m,lam_y,vex,vem,vey,
                       um,aa,psim,sigmaumi,uy,bb,cc,psiy,sigmauyi,ind);
  
  lx_name <- vex_name <- array('',dim = n_xitem)
  lm_name <- vem_name <- array('',dim = n_mitem)
  ly_name <- vey_name <- array('',dim = n_yitem)
  for (i in 1:n_xitem) {
    lx_name[i] <- paste0('lam_x',as.character(i))
    #vx_name[i] <- paste0('v_x',as.character(i))
    vex_name[i] <- paste0('ve_x',as.character(i))
  }
  for (i in 1:n_mitem) {
    lm_name[i] <- paste0('lam_m',as.character(i))
    #vm_name[i] <- paste0('v_m',as.character(i))
    vem_name[i] <- paste0('ve_m',as.character(i))
  }
  for (i in 1:n_yitem) {
    ly_name[i] <- paste0('lam_y',as.character(i))
    #vy_name[i] <- paste0('v_y',as.character(i))
    vey_name[i] <- paste0('ve_y',as.character(i))
  }
  um_name <- alpha_name <- uy_name <- beta_name <- gamma_name <- ind_name <- 
    array('',dim = 9);
  for (i in 1:9) {
    um_name[i] <- paste0('u_m0.',as.character(i))
    alpha_name[i] <- paste0('alpha0.',as.character(i))
    uy_name[i] <- paste0('u_y0.',as.character(i))
    beta_name[i] <- paste0('beta0.',as.character(i))
    gamma_name[i] <- paste0('gamma0.',as.character(i))
    ind_name[i] <- paste0('ind0.',as.character(i))
  }
  #  um_name1 <- alpha_name1 <- uy_name1 <- beta_name1 <- gamma_name1 <- 
  #    array('',dim = 3);
  #  um_name2 <- alpha_name2 <- uy_name2 <- beta_name2 <- gamma_name2 <- 
  #    array('',dim = nkot);
  #  for (i in 1:3) {
  #    um_name1[i] <- paste0('u_m',as.character(i-1))
  #    alpha_name1[i] <- paste0('alpha',as.character(i-1))
  #    uy_name1[i] <- paste0('u_y',as.character(i-1))
  #    beta_name1[i] <- paste0('beta',as.character(i-1))
  #    gamma_name1[i] <- paste0('gamma',as.character(i-1))
  #  }
  #  for (i in 1:nkot) {
  #    um_name2[i] <- paste0('d',as.character(i))
  #    alpha_name2[i] <- paste0('a',as.character(i))
  #    uy_name2[i] <- paste0('e',as.character(i))
  #    beta_name2[i] <- paste0('b',as.character(i))
  #    gamma_name2[i] <- paste0('c',as.character(i))
  #  }
  colnames(results.est) <- c(lx_name,lm_name,ly_name,vex_name,vem_name,vey_name,
                             um_name,alpha_name,"psi_m","sigma_umi",
                             uy_name,beta_name,gamma_name,"psi_y","sigma_uyi",ind_name);
  results.estimation <- t(results.est)
  sig <- array(dim = (2*n_xitem+2*n_mitem+2*n_yitem+9*6+4))
  output <- cbind(results.estimation,sig)
  output[results.estimation[,3]>0|results.estimation[,4]<0,5] <- 1
  output[results.estimation[,3]<0&results.estimation[,4]>0,5] <- 0
  
  results.burnin <- data.frame(lambda_x=lambda_x[(burnin+2):(iter+1),],
                               lambda_m=lambda_m[(burnin+2):(iter+1),],
                               lambda_y=lambda_y[(burnin+2):(iter+1),],
                               ve_x=ve_x[(burnin+2):(iter+1),],
                               ve_m=ve_m[(burnin+2):(iter+1),],
                               ve_y=ve_y[(burnin+2):(iter+1),],
                               u_m=abc1[(burnin+2):(iter+1),1:3],
                               dk=abc1[(burnin+2):(iter+1),4:(3+nkot)],
                               alpha=abc1[(burnin+2):(iter+1),(4+nkot):(3+nkot+3)],
                               ak=abc1[(burnin+2):(iter+1),(3+nkot+4):(2*(3+nkot))],
                               psi_m=psi_m[(burnin+2):(iter+1),],
                               sigma_umi=sigma_umi[(burnin+2):(iter+1),],
                               u_y=abc2[(burnin+2):(iter+1),1:3],
                               ek=abc2[(burnin+2):(iter+1),4:(3+nkot)],
                               beta=abc2[(burnin+2):(iter+1),(4+nkot):(3+nkot+3)],
                               bk=abc2[(burnin+2):(iter+1),(3+nkot+4):(2*(3+nkot))],
                               gamma=abc2[(burnin+2):(iter+1),(2*(3+nkot)+1):(2*(3+nkot)+3)],
                               ck=abc2[(burnin+2):(iter+1),(2*(3+nkot)+4):(3*(3+nkot))],
                               psi_y=psi_y[(burnin+2):(iter+1),],
                               sigma_uyi=sigma_uyi[(burnin+2):(iter+1),]);
  
  results.raw <- data.frame(lambda_x=lambda_x[2:(iter+1),],
                            lambda_m=lambda_m[2:(iter+1),],
                            lambda_y=lambda_y[2:(iter+1),],
                            ve_x=ve_x[2:(iter+1),],
                            ve_m=ve_m[2:(iter+1),],
                            ve_y=ve_y[2:(iter+1),],
                            u_m=abc1[2:(iter+1),1:3],
                            dk=abc1[2:(iter+1),4:(3+nkot)],
                            alpha=abc1[2:(iter+1),(4+nkot):(3+nkot+3)],
                            ak=abc1[2:(iter+1),(3+nkot+4):(2*(3+nkot))],
                            psi_m=psi_m[2:(iter+1),],
                            sigma_umi=sigma_umi[2:(iter+1),],
                            u_y=abc2[2:(iter+1),1:3],
                            ek=abc2[2:(iter+1),4:(3+nkot)],
                            beta=abc2[2:(iter+1),(4+nkot):(3+nkot+3)],
                            bk=abc2[2:(iter+1),(3+nkot+4):(2*(3+nkot))],
                            gamma=abc2[2:(iter+1),(2*(3+nkot)+1):(2*(3+nkot)+3)],
                            ck=abc2[2:(iter+1),(2*(3+nkot)+4):(3*(3+nkot))],
                            psi_y=psi_y[2:(iter+1),],
                            sigma_uyi=sigma_uyi[2:(iter+1),]);
  
  results <- list(est=output,raw=results.raw,remain=results.burnin)
  return(results)
}
