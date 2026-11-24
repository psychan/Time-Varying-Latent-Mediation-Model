n_sample = c(20,50,100,200)
m_measure = c(10,20)

time <- c(1:10)*0.1
path.true <- 1/(1+exp((-7)*(time-0.5)))-0.5
ind.true <- path.true*path.true
lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),10)
me.true <- 0.1
se.true <- 0.2

for (n in 2:4) {
  ##m=10
  wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n_sample[n]),
              'm10/sp')
  setwd(wd)
  
  lam <- me <- array(dim=c(90,3,200))
  aa <- bb <- cc <- ind <- array(dim=c(10,3,200))
  se <- array(dim=c(20,3,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      lam[,,i] <- as.matrix(data0[c(1:90),c(3,6:7)])
      me[,,i] <- as.matrix(data0[c(311:400),c(3,6:7)])
      se[,,i] <- as.matrix(data0[c(401:420),c(3,6:7)])
      for (j in 1:10) {
        aa[j,,i] <- as.matrix(data0[(88+j*3),c(3,6:7)])
        bb[j,,i] <- as.matrix(data0[(90+j*3),c(3,6:7)])
        cc[j,,i] <- as.matrix(data0[(89+j*3),c(3,6:7)])
        ind[j,,i] <- as.matrix(data0[(420+j),c(3,6:7)])
      }
    }
  }
  
  cover_rate <- mse <- power <- cvg <- array(dim=c(5,1))
  #cover_rate
  cr.lam <- cr.me <- array(0,dim=c(90,200))
  cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(10,200))
  cr.se <- array(0,dim=c(20,200))
  cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
  cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
  cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
  cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
  cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
  cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
  cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
  cr.lam <- cr.lam[,!is.na(lam[1,1,])]
  cr.me <- cr.me[,!is.na(lam[1,1,])]
  cr.a <- cr.a[,!is.na(lam[1,1,])]
  cr.b <- cr.b[,!is.na(lam[1,1,])]
  cr.c <- cr.c[,!is.na(lam[1,1,])]
  cr.se <- cr.se[,!is.na(lam[1,1,])]
  cr.ind <- cr.ind[,!is.na(lam[1,1,])]
  cover_rate[1] <- mean(cr.lam)
  cover_rate[2] <- mean(cr.me)
  cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
  cover_rate[4] <- mean(cr.se)
  cover_rate[5] <- mean(cr.ind)
  #mse
  mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
  mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
  mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                 mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                 mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
  mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
  mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
  #power
  tpr.lam <- tpr.me <- array(0,dim=c(90,200))
  tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(10,200))
  tpr.se <- array(0,dim=c(20,200))
  tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
  tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
  tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
  tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
  tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
  tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
  tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
  tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
  tpr.me <- tpr.me[,!is.na(lam[1,1,])]
  tpr.a <- tpr.a[,!is.na(lam[1,1,])]
  tpr.b <- tpr.b[,!is.na(lam[1,1,])]
  tpr.c <- tpr.c[,!is.na(lam[1,1,])]
  tpr.se <- tpr.se[,!is.na(lam[1,1,])]
  tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
  power[1] <- mean(tpr.lam)
  power[2] <- mean(tpr.me)
  power[3] <- mean(mean(tpr.a[c(1:4,6:10),]),mean(tpr.b[c(1:4,6:10),]),
                   mean(tpr.c[c(1:4,6:10),]))
  power[4] <- mean(tpr.se)
  power[5] <- mean(tpr.ind[c(1:4,6:10),])
  
  cvg <- mean(!is.na(lam[1,1,]))
  results <- cbind(cover_rate,mse,power,cvg)
  rownames(results) <- c('loadings','error','regression','residual','indirect')
  colnames(results) <- c('cover_rate','mse','power','convergence')
  
  out_name <- paste0('F:/course/statistic/tvm/S1/valley/sp_n',as.character(n_sample[n]),
                     'm10.txt')
  write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
}

time <- c(1:20)*0.05
path.true <- 1/(1+exp((-7)*(time-0.5)))-0.5
ind.true <- path.true*path.true
lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),20)
for (n in 3:4) {
  ##m=20
  wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n_sample[n]),
              'm20/sp')
  setwd(wd)
  
  lam <- me <- array(dim=c(180,3,200))
  aa <- bb <- cc <- ind <- array(dim=c(20,3,200))
  se <- array(dim=c(40,3,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      lam[,,i] <- as.matrix(data0[c(1:180),c(3,6:7)])
      me[,,i] <- as.matrix(data0[c(821:1000),c(3,6:7)])
      se[,,i] <- as.matrix(data0[c(1001:1040),c(3,6:7)])
      for (j in 1:20) {
        aa[j,,i] <- as.matrix(data0[(178+j*3),c(3,6:7)])
        bb[j,,i] <- as.matrix(data0[(180+j*3),c(3,6:7)])
        cc[j,,i] <- as.matrix(data0[(179+j*3),c(3,6:7)])
        ind[j,,i] <- as.matrix(data0[(1040+j),c(3,6:7)])
      }
    }
  }
  
  cover_rate <- mse <- power <- array(dim=c(5,1))
  #cover_rate
  cr.lam <- cr.me <- array(0,dim=c(180,200))
  cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(20,200))
  cr.se <- array(0,dim=c(40,200))
  cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
  cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
  cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
  cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
  cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
  cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
  cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
  cr.lam <- cr.lam[,!is.na(lam[1,1,])]
  cr.me <- cr.me[,!is.na(lam[1,1,])]
  cr.a <- cr.a[,!is.na(lam[1,1,])]
  cr.b <- cr.b[,!is.na(lam[1,1,])]
  cr.c <- cr.c[,!is.na(lam[1,1,])]
  cr.se <- cr.se[,!is.na(lam[1,1,])]
  cr.ind <- cr.ind[,!is.na(lam[1,1,])]
  cover_rate[1] <- mean(cr.lam)
  cover_rate[2] <- mean(cr.me)
  cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
  cover_rate[4] <- mean(cr.se)
  cover_rate[5] <- mean(cr.ind)
  #mse
  mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
  mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
  mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                 mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                 mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
  mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
  mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
  #power
  tpr.lam <- tpr.me <- array(0,dim=c(180,200))
  tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(20,200))
  tpr.se <- array(0,dim=c(40,200))
  tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
  tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
  tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
  tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
  tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
  tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
  tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
  tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
  tpr.me <- tpr.me[,!is.na(lam[1,1,])]
  tpr.a <- tpr.a[,!is.na(lam[1,1,])]
  tpr.b <- tpr.b[,!is.na(lam[1,1,])]
  tpr.c <- tpr.c[,!is.na(lam[1,1,])]
  tpr.se <- tpr.se[,!is.na(lam[1,1,])]
  tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
  power[1] <- mean(tpr.lam)
  power[2] <- mean(tpr.me)
  power[3] <- mean(mean(tpr.a[c(1:4,6:10),]),mean(tpr.b[c(1:4,6:10),]),
                   mean(tpr.c[c(1:4,6:10),]))
  power[4] <- mean(tpr.se)
  power[5] <- mean(tpr.ind[c(1:4,6:10),])
  
  cvg <- mean(!is.na(lam[1,1,]))
  results <- cbind(cover_rate,mse,power,cvg)
  rownames(results) <- c('loadings','error','regression','residual','indirect')
  colnames(results) <- c('cover_rate','mse','power','convergence')
  
  out_name <- paste0('F:/course/statistic/tvm/S1/valley/sp_n',as.character(n_sample[n]),
                     'm20.txt')
  write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
}

for (it in 1:1) {
  n_sample = c(20,50,100,200)
  m_measure = c(10,20)
  
  time <- c(1:10)*0.1
  path.true <- -0.1+0.5*exp(-time^0.5)
  ind.true <- path.true*path.true
  lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),10)
  me.true <- 0.1
  se.true <- 0.2
  
  for (n in 2:4) {
    ##m=10
    wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n_sample[n]),
                'm10/sp')
    setwd(wd)
    
    lam <- me <- array(dim=c(90,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(10,3,200))
    se <- array(dim=c(20,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      if(ncol(data0)==8){
        lam[,,i] <- as.matrix(data0[c(1:90),c(3,6:7)])
        me[,,i] <- as.matrix(data0[c(311:400),c(3,6:7)])
        se[,,i] <- as.matrix(data0[c(401:420),c(3,6:7)])
        for (j in 1:10) {
          aa[j,,i] <- as.matrix(data0[(88+j*3),c(3,6:7)])
          bb[j,,i] <- as.matrix(data0[(90+j*3),c(3,6:7)])
          cc[j,,i] <- as.matrix(data0[(89+j*3),c(3,6:7)])
          ind[j,,i] <- as.matrix(data0[(420+j),c(3,6:7)])
        }
      }
    }
    
    cover_rate <- mse <- power <- cvg <- array(dim=c(5,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(90,200))
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(10,200))
    cr.se <- array(0,dim=c(20,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cr.lam <- cr.lam[,!is.na(lam[1,1,])]
    cr.me <- cr.me[,!is.na(lam[1,1,])]
    cr.a <- cr.a[,!is.na(lam[1,1,])]
    cr.b <- cr.b[,!is.na(lam[1,1,])]
    cr.c <- cr.c[,!is.na(lam[1,1,])]
    cr.se <- cr.se[,!is.na(lam[1,1,])]
    cr.ind <- cr.ind[,!is.na(lam[1,1,])]
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[4] <- mean(cr.se)
    cover_rate[5] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
    mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
    mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
    mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
    mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(90,200))
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(10,200))
    tpr.se <- array(0,dim=c(20,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
    tpr.me <- tpr.me[,!is.na(lam[1,1,])]
    tpr.a <- tpr.a[,!is.na(lam[1,1,])]
    tpr.b <- tpr.b[,!is.na(lam[1,1,])]
    tpr.c <- tpr.c[,!is.na(lam[1,1,])]
    tpr.se <- tpr.se[,!is.na(lam[1,1,])]
    tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(mean(tpr.a),mean(tpr.b),mean(tpr.c))
    power[4] <- mean(tpr.se)
    power[5] <- mean(tpr.ind)
    
    cvg <- mean(!is.na(lam[1,1,]))
    results <- cbind(cover_rate,mse,power,cvg)
    rownames(results) <- c('loadings','error','regression','residual','indirect')
    colnames(results) <- c('cover_rate','mse','power','convergence')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/decay/sp_n',as.character(n_sample[n]),
                       'm10.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
  
  time <- c(1:20)*0.05
  path.true <- -0.1+0.5*exp(-time^0.5)
  ind.true <- path.true*path.true
  lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),20)
  for (n in 3:4) {
    ##m=20
    wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n_sample[n]),
                'm20/sp')
    setwd(wd)
    
    lam <- me <- array(dim=c(180,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(20,3,200))
    se <- array(dim=c(40,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      if(ncol(data0)==8){
        lam[,,i] <- as.matrix(data0[c(1:180),c(3,6:7)])
        me[,,i] <- as.matrix(data0[c(821:1000),c(3,6:7)])
        se[,,i] <- as.matrix(data0[c(1001:1040),c(3,6:7)])
        for (j in 1:20) {
          aa[j,,i] <- as.matrix(data0[(178+j*3),c(3,6:7)])
          bb[j,,i] <- as.matrix(data0[(180+j*3),c(3,6:7)])
          cc[j,,i] <- as.matrix(data0[(179+j*3),c(3,6:7)])
          ind[j,,i] <- as.matrix(data0[(1040+j),c(3,6:7)])
        }
      }
    }
    
    cover_rate <- mse <- power <- array(dim=c(5,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(180,200))
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(20,200))
    cr.se <- array(0,dim=c(40,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cr.lam <- cr.lam[,!is.na(lam[1,1,])]
    cr.me <- cr.me[,!is.na(lam[1,1,])]
    cr.a <- cr.a[,!is.na(lam[1,1,])]
    cr.b <- cr.b[,!is.na(lam[1,1,])]
    cr.c <- cr.c[,!is.na(lam[1,1,])]
    cr.se <- cr.se[,!is.na(lam[1,1,])]
    cr.ind <- cr.ind[,!is.na(lam[1,1,])]
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[4] <- mean(cr.se)
    cover_rate[5] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
    mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
    mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
    mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
    mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(180,200))
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(20,200))
    tpr.se <- array(0,dim=c(40,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
    tpr.me <- tpr.me[,!is.na(lam[1,1,])]
    tpr.a <- tpr.a[,!is.na(lam[1,1,])]
    tpr.b <- tpr.b[,!is.na(lam[1,1,])]
    tpr.c <- tpr.c[,!is.na(lam[1,1,])]
    tpr.se <- tpr.se[,!is.na(lam[1,1,])]
    tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(mean(tpr.a),mean(tpr.b),mean(tpr.c))
    power[4] <- mean(tpr.se)
    power[5] <- mean(tpr.ind[c(1:4,6:10),])
    
    cvg <- mean(!is.na(lam[1,1,]))
    results <- cbind(cover_rate,mse,power,cvg)
    rownames(results) <- c('loadings','error','regression','residual','indirect')
    colnames(results) <- c('cover_rate','mse','power','convergence')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/decay/sp_n',as.character(n_sample[n]),
                       'm20.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
  
}

for (it in 1:1) {
  n_sample = c(20,50,100,200)
  m_measure = c(10,20)
  
  time <- c(1:10)*0.1
  path.true <- 0.4
  ind.true <- path.true*path.true
  lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),10)
  me.true <- 0.1
  se.true <- 0.2
  
  for (n in 2:4) {
    ##m=10
    wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n_sample[n]),
                'm10/sp')
    setwd(wd)
    
    lam <- me <- array(dim=c(90,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(10,3,200))
    se <- array(dim=c(20,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      if(ncol(data0)==8){
        lam[,,i] <- as.matrix(data0[c(1:90),c(3,6:7)])
        me[,,i] <- as.matrix(data0[c(311:400),c(3,6:7)])
        se[,,i] <- as.matrix(data0[c(401:420),c(3,6:7)])
        for (j in 1:10) {
          aa[j,,i] <- as.matrix(data0[(88+j*3),c(3,6:7)])
          bb[j,,i] <- as.matrix(data0[(90+j*3),c(3,6:7)])
          cc[j,,i] <- as.matrix(data0[(89+j*3),c(3,6:7)])
          ind[j,,i] <- as.matrix(data0[(420+j),c(3,6:7)])
        }
      }
    }
    
    cover_rate <- mse <- power <- cvg <- array(dim=c(5,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(90,200))
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(10,200))
    cr.se <- array(0,dim=c(20,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cr.lam <- cr.lam[,!is.na(lam[1,1,])]
    cr.me <- cr.me[,!is.na(lam[1,1,])]
    cr.a <- cr.a[,!is.na(lam[1,1,])]
    cr.b <- cr.b[,!is.na(lam[1,1,])]
    cr.c <- cr.c[,!is.na(lam[1,1,])]
    cr.se <- cr.se[,!is.na(lam[1,1,])]
    cr.ind <- cr.ind[,!is.na(lam[1,1,])]
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[4] <- mean(cr.se)
    cover_rate[5] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
    mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
    mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
    mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
    mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(90,200))
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(10,200))
    tpr.se <- array(0,dim=c(20,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
    tpr.me <- tpr.me[,!is.na(lam[1,1,])]
    tpr.a <- tpr.a[,!is.na(lam[1,1,])]
    tpr.b <- tpr.b[,!is.na(lam[1,1,])]
    tpr.c <- tpr.c[,!is.na(lam[1,1,])]
    tpr.se <- tpr.se[,!is.na(lam[1,1,])]
    tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(mean(tpr.a),mean(tpr.b),mean(tpr.c))
    power[4] <- mean(tpr.se)
    power[5] <- mean(tpr.ind[c(1:4,6:10),])
    
    cvg <- mean(!is.na(lam[1,1,]))
    results <- cbind(cover_rate,mse,power,cvg)
    rownames(results) <- c('loadings','error','regression','residual','indirect')
    colnames(results) <- c('cover_rate','mse','power','convergence')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/flat/sp_n',as.character(n_sample[n]),
                       'm10.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
  
  time <- c(1:20)*0.05
  path.true <- 0.4
  ind.true <- path.true*path.true
  lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),20)
  for (n in 3:4) {
    ##m=20
    wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n_sample[n]),
                'm20/sp')
    setwd(wd)
    
    lam <- me <- array(dim=c(180,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(20,3,200))
    se <- array(dim=c(40,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      if(ncol(data0)==8){
        lam[,,i] <- as.matrix(data0[c(1:180),c(3,6:7)])
        me[,,i] <- as.matrix(data0[c(821:1000),c(3,6:7)])
        se[,,i] <- as.matrix(data0[c(1001:1040),c(3,6:7)])
        for (j in 1:20) {
          aa[j,,i] <- as.matrix(data0[(178+j*3),c(3,6:7)])
          bb[j,,i] <- as.matrix(data0[(180+j*3),c(3,6:7)])
          cc[j,,i] <- as.matrix(data0[(179+j*3),c(3,6:7)])
          ind[j,,i] <- as.matrix(data0[(1040+j),c(3,6:7)])
        }
      }
    }
    
    cover_rate <- mse <- power <- array(dim=c(5,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(180,200))
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(20,200))
    cr.se <- array(0,dim=c(40,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cr.lam <- cr.lam[,!is.na(lam[1,1,])]
    cr.me <- cr.me[,!is.na(lam[1,1,])]
    cr.a <- cr.a[,!is.na(lam[1,1,])]
    cr.b <- cr.b[,!is.na(lam[1,1,])]
    cr.c <- cr.c[,!is.na(lam[1,1,])]
    cr.se <- cr.se[,!is.na(lam[1,1,])]
    cr.ind <- cr.ind[,!is.na(lam[1,1,])]
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[4] <- mean(cr.se)
    cover_rate[5] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
    mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
    mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
    mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
    mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(180,200))
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(20,200))
    tpr.se <- array(0,dim=c(40,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
    tpr.me <- tpr.me[,!is.na(lam[1,1,])]
    tpr.a <- tpr.a[,!is.na(lam[1,1,])]
    tpr.b <- tpr.b[,!is.na(lam[1,1,])]
    tpr.c <- tpr.c[,!is.na(lam[1,1,])]
    tpr.se <- tpr.se[,!is.na(lam[1,1,])]
    tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(mean(tpr.a),mean(tpr.b),mean(tpr.c))
    power[4] <- mean(tpr.se)
    power[5] <- mean(tpr.ind[c(1:4,6:10),])
    
    cvg <- mean(!is.na(lam[1,1,]))
    results <- cbind(cover_rate,mse,power,cvg)
    rownames(results) <- c('loadings','error','regression','residual','indirect')
    colnames(results) <- c('cover_rate','mse','power','convergence')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/flat/sp_n',as.character(n_sample[n]),
                       'm20.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
  
}

for (it in 1:1) {
  n_sample = c(20,50,100,200)
  m_measure = c(10,20)
  
  time <- c(1:10)*0.1
  path.true <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
  ind.true <- path.true*path.true
  lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),10)
  me.true <- 0.1
  se.true <- 0.2
  
  for (n in 2:4) {
    ##m=10
    wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n_sample[n]),
                'm10/sp')
    setwd(wd)
    
    lam <- me <- array(dim=c(90,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(10,3,200))
    se <- array(dim=c(20,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      if(ncol(data0)==8){
        lam[,,i] <- as.matrix(data0[c(1:90),c(3,6:7)])
        me[,,i] <- as.matrix(data0[c(311:400),c(3,6:7)])
        se[,,i] <- as.matrix(data0[c(401:420),c(3,6:7)])
        for (j in 1:10) {
          aa[j,,i] <- as.matrix(data0[(88+j*3),c(3,6:7)])
          bb[j,,i] <- as.matrix(data0[(90+j*3),c(3,6:7)])
          cc[j,,i] <- as.matrix(data0[(89+j*3),c(3,6:7)])
          ind[j,,i] <- as.matrix(data0[(420+j),c(3,6:7)])
        }
      }
    }
    
    cover_rate <- mse <- power <- cvg <- array(dim=c(5,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(90,200))
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(10,200))
    cr.se <- array(0,dim=c(20,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cr.lam <- cr.lam[,!is.na(lam[1,1,])]
    cr.me <- cr.me[,!is.na(lam[1,1,])]
    cr.a <- cr.a[,!is.na(lam[1,1,])]
    cr.b <- cr.b[,!is.na(lam[1,1,])]
    cr.c <- cr.c[,!is.na(lam[1,1,])]
    cr.se <- cr.se[,!is.na(lam[1,1,])]
    cr.ind <- cr.ind[,!is.na(lam[1,1,])]
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[4] <- mean(cr.se)
    cover_rate[5] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
    mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
    mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
    mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
    mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(90,200))
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(10,200))
    tpr.se <- array(0,dim=c(20,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
    tpr.me <- tpr.me[,!is.na(lam[1,1,])]
    tpr.a <- tpr.a[,!is.na(lam[1,1,])]
    tpr.b <- tpr.b[,!is.na(lam[1,1,])]
    tpr.c <- tpr.c[,!is.na(lam[1,1,])]
    tpr.se <- tpr.se[,!is.na(lam[1,1,])]
    tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(mean(tpr.a[c(1:4,6:10),]),mean(tpr.b[c(1:4,6:10),]),
                     mean(tpr.c[c(1:4,6:10),]))
    power[4] <- mean(tpr.se)
    power[5] <- mean(tpr.ind[c(1:4,6:10),])
    
    cvg <- mean(!is.na(lam[1,1,]))
    results <- cbind(cover_rate,mse,power,cvg)
    rownames(results) <- c('loadings','error','regression','residual','indirect')
    colnames(results) <- c('cover_rate','mse','power','convergence')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/hump/sp_n',as.character(n_sample[n]),
                       'm10.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
  
  time <- c(1:20)*0.05
  path.true <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
  ind.true <- path.true*path.true
  lam.true <- rep(c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8),20)
  for (n in 3:4) {
    ##m=20
    wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n_sample[n]),
                'm20/sp')
    setwd(wd)
    
    lam <- me <- array(dim=c(180,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(20,3,200))
    se <- array(dim=c(40,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      if(ncol(data0)==8){
        lam[,,i] <- as.matrix(data0[c(1:180),c(3,6:7)])
        me[,,i] <- as.matrix(data0[c(821:1000),c(3,6:7)])
        se[,,i] <- as.matrix(data0[c(1001:1040),c(3,6:7)])
        for (j in 1:20) {
          aa[j,,i] <- as.matrix(data0[(178+j*3),c(3,6:7)])
          bb[j,,i] <- as.matrix(data0[(180+j*3),c(3,6:7)])
          cc[j,,i] <- as.matrix(data0[(179+j*3),c(3,6:7)])
          ind[j,,i] <- as.matrix(data0[(1040+j),c(3,6:7)])
        }
      }
    }
    
    cover_rate <- mse <- power <- array(dim=c(5,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(180,200))
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(20,200))
    cr.se <- array(0,dim=c(40,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se[se[,2,]<=se.true&se[,3,]>=se.true] <- 1
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cr.lam <- cr.lam[,!is.na(lam[1,1,])]
    cr.me <- cr.me[,!is.na(lam[1,1,])]
    cr.a <- cr.a[,!is.na(lam[1,1,])]
    cr.b <- cr.b[,!is.na(lam[1,1,])]
    cr.c <- cr.c[,!is.na(lam[1,1,])]
    cr.se <- cr.se[,!is.na(lam[1,1,])]
    cr.ind <- cr.ind[,!is.na(lam[1,1,])]
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[4] <- mean(cr.se)
    cover_rate[5] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,!is.na(lam[1,1,])]-lam.true)^2)
    mse[2] <- mean((me[,1,!is.na(lam[1,1,])]-me.true)^2)
    mse[3] <- mean(mean((aa[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((bb[,1,!is.na(lam[1,1,])]-path.true)^2),
                   mean((cc[,1,!is.na(lam[1,1,])]-path.true)^2))
    mse[4] <- mean((se[,1,!is.na(lam[1,1,])]-se.true)^2)
    mse[5] <- mean((ind[,1,!is.na(lam[1,1,])]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(180,200))
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(20,200))
    tpr.se <- array(0,dim=c(40,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se[se[,2,]>=0|se[,3,]<=0] <- 1
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    tpr.lam <- tpr.lam[,!is.na(lam[1,1,])]
    tpr.me <- tpr.me[,!is.na(lam[1,1,])]
    tpr.a <- tpr.a[,!is.na(lam[1,1,])]
    tpr.b <- tpr.b[,!is.na(lam[1,1,])]
    tpr.c <- tpr.c[,!is.na(lam[1,1,])]
    tpr.se <- tpr.se[,!is.na(lam[1,1,])]
    tpr.ind <- tpr.ind[,!is.na(lam[1,1,])]
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(mean(tpr.a[c(1:4,6:10),]),mean(tpr.b[c(1:4,6:10),]),
                     mean(tpr.c[c(1:4,6:10),]))
    power[4] <- mean(tpr.se)
    power[5] <- mean(tpr.ind[c(1:4,6:10),])
    
    cvg <- mean(!is.na(lam[1,1,]))
    results <- cbind(cover_rate,mse,power,cvg)
    rownames(results) <- c('loadings','error','regression','residual','indirect')
    colnames(results) <- c('cover_rate','mse','power','convergence')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/hump/sp_n',as.character(n_sample[n]),
                       'm20.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
  
}
