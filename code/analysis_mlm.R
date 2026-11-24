n_sample = c(20,50,100,200)
m_measure = c(10,20)

time <- c(1:9)*0.1
path.true <- 1/(1+exp((-7)*(time-0.5)))-0.5
ind.true <- path.true*path.true
lam.true <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8)
me.true <- rep(0.1,9)
int.true <- c(0.4,0.6)
se.b.true <- c(0.2,0.2)
se.w.true <- c(0.1,0.1)

for (n in 1:4) {
  for (m in 1:2) {
    wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n_sample[n]),
                'm',as.character(m_measure[m]),'/multilevel')
    setwd(wd)
    
    lam <- me <- array(dim=c(9,3,200))
    se.b <- se.w <- int <- array(dim=c(2,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(9,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      lam[,,i] <- as.matrix(data0[c(1:9),c(3,6:7)])
      me[,,i] <- as.matrix(data0[c(16:24),c(3,6:7)])
      int[,,i] <- as.matrix(data0[13:14,c(3,6:7)])
      se.b[,,i] <- as.matrix(data0[c(25:26),c(3,6:7)])
      se.w[,,i] <- as.matrix(data0[c(60:61),c(3,6:7)])
      for (j in 1:9) {
        aa[j,,i] <- as.matrix(data0[10,c(3,6:7)])
        bb[j,,i] <- as.matrix(data0[12,c(3,6:7)])
        cc[j,,i] <- as.matrix(data0[11,c(3,6:7)])
        ind[j,,i] <- as.matrix(data0[62,c(3,6:7)])
      }
    }
    cover_rate <- mse <- power <- array(dim=c(6,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(9,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se.b <- cr.se.w <- cr.int <- array(0,dim=c(2,200))
    cr.se.b[se.b[,2,]<=se.b.true&se.b[,3,]>=se.b.true] <- 1
    cr.se.w[se.w[,2,]<=se.w.true&se.w[,3,]>=se.w.true] <- 1
    cr.int[int[,2,]<=int.true&int[,3,]>=int.true] <- 1
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(9,200))
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(cr.int)
    cover_rate[4] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[5] <- mean(mean(cr.se.b),mean(cr.se.w))
    cover_rate[6] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,]-lam.true)^2)
    mse[2] <- mean((me[,1,]-me.true)^2)
    mse[3] <- mean((int[,1,]-int.true)^2)
    mse[4] <- mean(mean((aa[,1,]-path.true)^2),mean((bb[,1,]-path.true)^2),
                   mean((cc[,1,]-path.true)^2))
    mse[5] <- mean(mean((se.b[,1,]-se.b.true)^2),mean((se.w[,1,]-se.w.true)^2))
    mse[6] <- mean((ind[,1,]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(9,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se.b <- tpr.se.w <- tpr.int <- array(0,dim=c(2,200))
    tpr.se.b[se.b[,2,]>=0|se.b[,3,]<=0] <- 1
    tpr.se.w[se.w[,2,]>=0|se.w[,3,]<=0] <- 1
    tpr.int[int[,2,]>=0|int[,3,]<=0] <- 1
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(9,200))
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(tpr.int)
    power[4] <- mean(mean(tpr.a[c(1:4,6:9),]),mean(tpr.b[c(1:4,6:9),]),
                     mean(tpr.c[c(1:4,6:9),]))
    power[5] <- mean(mean(tpr.se.b),mean(tpr.se.w))
    power[6] <- mean(tpr.ind[c(1:4,6:9),])
    
    results <- cbind(cover_rate,mse,power)
    rownames(results) <- c('loadings','error','intercept','regression',
                           'residual','indirect')
    colnames(results) <- c('cover_rate','mse','power')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/valley/multilevel_n',as.character(n_sample[n]),
                       'm',as.character(m_measure[m]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}

n_sample = c(20,50,100,200)
m_measure = c(10,20)

time <- c(1:9)*0.1
path.true <- 1/(1+exp((-7)*(time-0.5)))-0.5
ind.true <- path.true*path.true
lam.true <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8)
me.true <- rep(0.1,9)
int.true <- c(0.4,0.6)
se.b.true <- c(0.2,0.2)
se.w.true <- c(0.1,0.1)

for (n in 1:4) {
  for (m in 1:2) {
    wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n_sample[n]),
                'm',as.character(m_measure[m]),'/multilevel')
    setwd(wd)
    
    lam <- me <- array(dim=c(9,3,200))
    se.b <- se.w <- int <- array(dim=c(2,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(9,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      lam[,,i] <- as.matrix(data0[c(1:9),c(3,6:7)])
      me[,,i] <- as.matrix(data0[c(16:24),c(3,6:7)])
      int[,,i] <- as.matrix(data0[13:14,c(3,6:7)])
      se.b[,,i] <- as.matrix(data0[c(25:26),c(3,6:7)])
      se.w[,,i] <- as.matrix(data0[c(60:61),c(3,6:7)])
      for (j in 1:9) {
        aa[j,,i] <- as.matrix(data0[10,c(3,6:7)])
        bb[j,,i] <- as.matrix(data0[12,c(3,6:7)])
        cc[j,,i] <- as.matrix(data0[11,c(3,6:7)])
        ind[j,,i] <- as.matrix(data0[62,c(3,6:7)])
      }
    }
    cover_rate <- mse <- power <- array(dim=c(6,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(9,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se.b <- cr.se.w <- cr.int <- array(0,dim=c(2,200))
    cr.se.b[se.b[,2,]<=se.b.true&se.b[,3,]>=se.b.true] <- 1
    cr.se.w[se.w[,2,]<=se.w.true&se.w[,3,]>=se.w.true] <- 1
    cr.int[int[,2,]<=int.true&int[,3,]>=int.true] <- 1
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(9,200))
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(cr.int)
    cover_rate[4] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[5] <- mean(mean(cr.se.b),mean(cr.se.w))
    cover_rate[6] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,]-lam.true)^2)
    mse[2] <- mean((me[,1,]-me.true)^2)
    mse[3] <- mean((int[,1,]-int.true)^2)
    mse[4] <- mean(mean((aa[,1,]-path.true)^2),mean((bb[,1,]-path.true)^2),
                   mean((cc[,1,]-path.true)^2))
    mse[5] <- mean(mean((se.b[,1,]-se.b.true)^2),mean((se.w[,1,]-se.w.true)^2))
    mse[6] <- mean((ind[,1,]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(9,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se.b <- tpr.se.w <- tpr.int <- array(0,dim=c(2,200))
    tpr.se.b[se.b[,2,]>=0|se.b[,3,]<=0] <- 1
    tpr.se.w[se.w[,2,]>=0|se.w[,3,]<=0] <- 1
    tpr.int[int[,2,]>=0|int[,3,]<=0] <- 1
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(9,200))
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(tpr.int)
    power[4] <- mean(mean(tpr.a[c(1:4,6:9),]),mean(tpr.b[c(1:4,6:9),]),
                     mean(tpr.c[c(1:4,6:9),]))
    power[5] <- mean(mean(tpr.se.b),mean(tpr.se.w))
    power[6] <- mean(tpr.ind[c(1:4,6:9),])
    
    results <- cbind(cover_rate,mse,power)
    rownames(results) <- c('loadings','error','intercept','regression',
                           'residual','indirect')
    colnames(results) <- c('cover_rate','mse','power')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/valley/multilevel_n',as.character(n_sample[n]),
                       'm',as.character(m_measure[m]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}

n_sample = c(20,50,100,200)
m_measure = c(10,20)

time <- c(1:9)*0.1
path.true <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
ind.true <- path.true*path.true
lam.true <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8)
me.true <- rep(0.1,9)
int.true <- c(0.4,0.6)
se.b.true <- c(0.2,0.2)
se.w.true <- c(0.1,0.1)

for (n in 1:4) {
  for (m in 1:2) {
    wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n_sample[n]),
                'm',as.character(m_measure[m]),'/multilevel')
    setwd(wd)
    
    lam <- me <- array(dim=c(9,3,200))
    se.b <- se.w <- int <- array(dim=c(2,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(9,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      lam[,,i] <- as.matrix(data0[c(1:9),c(3,6:7)])
      me[,,i] <- as.matrix(data0[c(16:24),c(3,6:7)])
      int[,,i] <- as.matrix(data0[13:14,c(3,6:7)])
      se.b[,,i] <- as.matrix(data0[c(25:26),c(3,6:7)])
      se.w[,,i] <- as.matrix(data0[c(60:61),c(3,6:7)])
      for (j in 1:9) {
        aa[j,,i] <- as.matrix(data0[10,c(3,6:7)])
        bb[j,,i] <- as.matrix(data0[12,c(3,6:7)])
        cc[j,,i] <- as.matrix(data0[11,c(3,6:7)])
        ind[j,,i] <- as.matrix(data0[62,c(3,6:7)])
      }
    }
    cover_rate <- mse <- power <- array(dim=c(6,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(9,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se.b <- cr.se.w <- cr.int <- array(0,dim=c(2,200))
    cr.se.b[se.b[,2,]<=se.b.true&se.b[,3,]>=se.b.true] <- 1
    cr.se.w[se.w[,2,]<=se.w.true&se.w[,3,]>=se.w.true] <- 1
    cr.int[int[,2,]<=int.true&int[,3,]>=int.true] <- 1
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(9,200))
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(cr.int)
    cover_rate[4] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[5] <- mean(mean(cr.se.b),mean(cr.se.w))
    cover_rate[6] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,]-lam.true)^2)
    mse[2] <- mean((me[,1,]-me.true)^2)
    mse[3] <- mean((int[,1,]-int.true)^2)
    mse[4] <- mean(mean((aa[,1,]-path.true)^2),mean((bb[,1,]-path.true)^2),
                   mean((cc[,1,]-path.true)^2))
    mse[5] <- mean(mean((se.b[,1,]-se.b.true)^2),mean((se.w[,1,]-se.w.true)^2))
    mse[6] <- mean((ind[,1,]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(9,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se.b <- tpr.se.w <- tpr.int <- array(0,dim=c(2,200))
    tpr.se.b[se.b[,2,]>=0|se.b[,3,]<=0] <- 1
    tpr.se.w[se.w[,2,]>=0|se.w[,3,]<=0] <- 1
    tpr.int[int[,2,]>=0|int[,3,]<=0] <- 1
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(9,200))
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(tpr.int)
    power[4] <- mean(mean(tpr.a),mean(tpr.b),mean(tpr.c))
    power[5] <- mean(mean(tpr.se.b),mean(tpr.se.w))
    power[6] <- mean(tpr.ind)
    
    results <- cbind(cover_rate,mse,power)
    rownames(results) <- c('loadings','error','intercept','regression',
                           'residual','indirect')
    colnames(results) <- c('cover_rate','mse','power')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/hump/multilevel_n',as.character(n_sample[n]),
                       'm',as.character(m_measure[m]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}

n_sample = c(20,50,100,200)
m_measure = c(10,20)

time <- c(1:9)*0.1
path.true <- rep(0.4,9)
ind.true <- path.true*path.true
lam.true <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8)
me.true <- rep(0.1,9)
int.true <- c(0.4,0.6)
se.b.true <- c(0.2,0.2)
se.w.true <- c(0.1,0.1)

for (n in 1:4) {
  for (m in 1:2) {
    wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n_sample[n]),
                'm',as.character(m_measure[m]),'/multilevel')
    setwd(wd)
    
    lam <- me <- array(dim=c(9,3,200))
    se.b <- se.w <- int <- array(dim=c(2,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(9,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      lam[,,i] <- as.matrix(data0[c(1:9),c(3,6:7)])
      me[,,i] <- as.matrix(data0[c(16:24),c(3,6:7)])
      int[,,i] <- as.matrix(data0[13:14,c(3,6:7)])
      se.b[,,i] <- as.matrix(data0[c(25:26),c(3,6:7)])
      se.w[,,i] <- as.matrix(data0[c(60:61),c(3,6:7)])
      for (j in 1:9) {
        aa[j,,i] <- as.matrix(data0[10,c(3,6:7)])
        bb[j,,i] <- as.matrix(data0[12,c(3,6:7)])
        cc[j,,i] <- as.matrix(data0[11,c(3,6:7)])
        ind[j,,i] <- as.matrix(data0[62,c(3,6:7)])
      }
    }
    cover_rate <- mse <- power <- array(dim=c(6,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(9,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se.b <- cr.se.w <- cr.int <- array(0,dim=c(2,200))
    cr.se.b[se.b[,2,]<=se.b.true&se.b[,3,]>=se.b.true] <- 1
    cr.se.w[se.w[,2,]<=se.w.true&se.w[,3,]>=se.w.true] <- 1
    cr.int[int[,2,]<=int.true&int[,3,]>=int.true] <- 1
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(9,200))
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(cr.int)
    cover_rate[4] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[5] <- mean(mean(cr.se.b),mean(cr.se.w))
    cover_rate[6] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,]-lam.true)^2)
    mse[2] <- mean((me[,1,]-me.true)^2)
    mse[3] <- mean((int[,1,]-int.true)^2)
    mse[4] <- mean(mean((aa[,1,]-path.true)^2),mean((bb[,1,]-path.true)^2),
                   mean((cc[,1,]-path.true)^2))
    mse[5] <- mean(mean((se.b[,1,]-se.b.true)^2),mean((se.w[,1,]-se.w.true)^2))
    mse[6] <- mean((ind[,1,]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(9,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se.b <- tpr.se.w <- tpr.int <- array(0,dim=c(2,200))
    tpr.se.b[se.b[,2,]>=0|se.b[,3,]<=0] <- 1
    tpr.se.w[se.w[,2,]>=0|se.w[,3,]<=0] <- 1
    tpr.int[int[,2,]>=0|int[,3,]<=0] <- 1
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(9,200))
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(tpr.int)
    power[4] <- mean(mean(tpr.a),mean(tpr.b),mean(tpr.c))
    power[5] <- mean(mean(tpr.se.b),mean(tpr.se.w))
    power[6] <- mean(tpr.ind)
    
    results <- cbind(cover_rate,mse,power)
    rownames(results) <- c('loadings','error','intercept','regression',
                           'residual','indirect')
    colnames(results) <- c('cover_rate','mse','power')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/flat/multilevel_n',as.character(n_sample[n]),
                       'm',as.character(m_measure[m]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}

n_sample = c(20,50,100,200)
m_measure = c(10,20)

time <- c(1:9)*0.1
path.true <- -0.1+0.5*exp(-time^0.5)
ind.true <- path.true*path.true
lam.true <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8)
me.true <- rep(0.1,9)
int.true <- c(0.4,0.6)
se.b.true <- c(0.2,0.2)
se.w.true <- c(0.1,0.1)

for (n in 1:4) {
  for (m in 1:2) {
    wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n_sample[n]),
                'm',as.character(m_measure[m]),'/multilevel')
    setwd(wd)
    
    lam <- me <- array(dim=c(9,3,200))
    se.b <- se.w <- int <- array(dim=c(2,3,200))
    aa <- bb <- cc <- ind <- array(dim=c(9,3,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      lam[,,i] <- as.matrix(data0[c(1:9),c(3,6:7)])
      me[,,i] <- as.matrix(data0[c(16:24),c(3,6:7)])
      int[,,i] <- as.matrix(data0[13:14,c(3,6:7)])
      se.b[,,i] <- as.matrix(data0[c(25:26),c(3,6:7)])
      se.w[,,i] <- as.matrix(data0[c(60:61),c(3,6:7)])
      for (j in 1:9) {
        aa[j,,i] <- as.matrix(data0[10,c(3,6:7)])
        bb[j,,i] <- as.matrix(data0[12,c(3,6:7)])
        cc[j,,i] <- as.matrix(data0[11,c(3,6:7)])
        ind[j,,i] <- as.matrix(data0[62,c(3,6:7)])
      }
    }
    cover_rate <- mse <- power <- array(dim=c(6,1))
    #cover_rate
    cr.lam <- cr.me <- array(0,dim=c(9,200))
    cr.lam[lam[,2,]<=lam.true&lam[,3,]>=lam.true] <- 1
    cr.me[me[,2,]<=me.true&me[,3,]>=me.true] <- 1
    cr.se.b <- cr.se.w <- cr.int <- array(0,dim=c(2,200))
    cr.se.b[se.b[,2,]<=se.b.true&se.b[,3,]>=se.b.true] <- 1
    cr.se.w[se.w[,2,]<=se.w.true&se.w[,3,]>=se.w.true] <- 1
    cr.int[int[,2,]<=int.true&int[,3,]>=int.true] <- 1
    cr.a <- cr.b <- cr.c <- cr.ind <- array(0,dim=c(9,200))
    cr.a[aa[,2,]<=path.true&aa[,2,]>=path.true] <- 1
    cr.b[bb[,2,]<=path.true&bb[,2,]>=path.true] <- 1
    cr.c[cc[,2,]<=path.true&cc[,2,]>=path.true] <- 1
    cr.ind[ind[,2,]<=ind.true&ind[,2,]>=ind.true] <- 1
    cover_rate[1] <- mean(cr.lam)
    cover_rate[2] <- mean(cr.me)
    cover_rate[3] <- mean(cr.int)
    cover_rate[4] <- mean(mean(cr.a),mean(cr.b),mean(cr.c))
    cover_rate[5] <- mean(mean(cr.se.b),mean(cr.se.w))
    cover_rate[6] <- mean(cr.ind)
    #mse
    mse[1] <- mean((lam[,1,]-lam.true)^2)
    mse[2] <- mean((me[,1,]-me.true)^2)
    mse[3] <- mean((int[,1,]-int.true)^2)
    mse[4] <- mean(mean((aa[,1,]-path.true)^2),mean((bb[,1,]-path.true)^2),
                   mean((cc[,1,]-path.true)^2))
    mse[5] <- mean(mean((se.b[,1,]-se.b.true)^2),mean((se.w[,1,]-se.w.true)^2))
    mse[6] <- mean((ind[,1,]-ind.true)^2)
    #power
    tpr.lam <- tpr.me <- array(0,dim=c(9,200))
    tpr.lam[lam[,2,]>=0|lam[,3,]<=0] <- 1
    tpr.me[me[,2,]>=0|me[,3,]<=0] <- 1
    tpr.se.b <- tpr.se.w <- tpr.int <- array(0,dim=c(2,200))
    tpr.se.b[se.b[,2,]>=0|se.b[,3,]<=0] <- 1
    tpr.se.w[se.w[,2,]>=0|se.w[,3,]<=0] <- 1
    tpr.int[int[,2,]>=0|int[,3,]<=0] <- 1
    tpr.a <- tpr.b <- tpr.c <- tpr.ind <- array(0,dim=c(9,200))
    tpr.a[aa[,2,]>=0|aa[,3,]<=0] <- 1
    tpr.b[bb[,2,]>=0|bb[,3,]<=0] <- 1
    tpr.c[cc[,2,]>=0|cc[,3,]<=0] <- 1
    tpr.ind[ind[,2,]>=0|ind[,3,]<=0] <- 1
    power[1] <- mean(tpr.lam)
    power[2] <- mean(tpr.me)
    power[3] <- mean(tpr.int)
    power[4] <- mean(mean(tpr.a),mean(tpr.b),mean(tpr.c))
    power[5] <- mean(mean(tpr.se.b),mean(tpr.se.w))
    power[6] <- mean(tpr.ind)
    
    results <- cbind(cover_rate,mse,power)
    rownames(results) <- c('loadings','error','intercept','regression',
                           'residual','indirect')
    colnames(results) <- c('cover_rate','mse','power')
    
    out_name <- paste0('F:/course/statistic/tvm/S1/decay/multilevel_n',as.character(n_sample[n]),
                       'm',as.character(m_measure[m]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}
