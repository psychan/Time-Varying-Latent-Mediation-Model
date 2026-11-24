library(foreach)
library(iterators)
n <- c(20,50,100,200)
m <- c(10,20)
foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    data<-array(dim=c(76,5,200))
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data[,,i] <- data1
    }
    rn <- rownames(data1)
    cn <- colnames(data1)
    dimnames(data)[[1]] <- rn
    dimnames(data)[[2]] <- cn
    #真值
    ##回归系数真值
    time <- c(1:9)*0.1
    beta <- 1/(1+exp((-7)*(time-0.5)))-0.5
    ##间接效应真值
    ind <- beta*beta
    true_value <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8,
                    rep(0.1,9),rep(0.4,9),beta,0.2,0.1,
                    rep(0.6,9),beta,beta,0.2,0.1,
                    ind)
    #区间覆盖率
    cr <- array(dim=c(76,200))
    cr[data[,3,]<=true_value&data[,4,]>=true_value] <- 1
    cr[is.na(cr)] <- 0
    rownames(cr) <- rn
    cover_rate <- rowMeans(cr)
    cr_se <- sqrt(cover_rate*(1-cover_rate)/200)
    #MSE
    mse <- rowMeans((true_value-data[,1,])^2)
    mse_se <- sqrt(rowMeans(((true_value-data[,1,])^2-mse)^2)/(200*199))
    #relative bias
    relative_bias <- rowMeans((data[,1,]-true_value)/true_value)
    rb_se <- sqrt(rowSums(((data[,1,]-true_value)/true_value-relative_bias)^2)/(200*199))
    #power
    tpr <- tnr <- array(dim=c(76,200))
    true_positive <- true_value!=0
    tpr[data[,5,]==true_positive] <- 1
    tpr[is.na(tpr)] <- 0
    tpr[!true_positive,] <- NA
    tnr[data[,5,]==true_positive] <- 1
    tnr[is.na(tnr)] <- 0
    tnr[true_positive,] <- NA
    rownames(tpr) <- rn
    true_positive_rate <- rowMeans(tpr)
    true_negative_rate <- rowMeans(tnr)
    tpr_se <- sqrt(true_positive_rate*(1-true_positive_rate)/200)
    tnr_se <- sqrt(true_negative_rate*(1-true_negative_rate)/200)
    #总结果
    results <- cbind(cover_rate,mse,relative_bias,true_positive_rate,true_negative_rate)
    results_se <- cbind(cr_se,mse_se,rb_se,tpr_se,tnr_se)
    #输出
    out_name <- paste0('F:/course/statistic/tvm/S1/valley/est_results_n',
                       as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
    f1 <- paste0('F:/course/statistic/tvm/S1/valley/est_se_n',
                 as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results_se,file=f1,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}

foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    data<-array(dim=c(76,5,200))
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data[,,i] <- data1
    }
    rn <- rownames(data1)
    cn <- colnames(data1)
    dimnames(data)[[1]] <- rn
    dimnames(data)[[2]] <- cn
    #真值
    ##回归系数真值
    time <- c(1:9)*0.1
    beta <- -0.1+0.5*exp(-time^0.5)
    ##间接效应真值
    ind <- beta*beta
    true_value <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8,
                    rep(0.1,9),rep(0.4,9),beta,0.2,0.1,
                    rep(0.6,9),beta,beta,0.2,0.1,
                    ind)
    #区间覆盖率
    cr <- array(dim=c(76,200))
    cr[data[,3,]<=true_value&data[,4,]>=true_value] <- 1
    cr[is.na(cr)] <- 0
    rownames(cr) <- rn
    cover_rate <- rowMeans(cr)
    cr_se <- sqrt(cover_rate*(1-cover_rate)/200)
    #MSE
    mse <- rowMeans((true_value-data[,1,])^2)
    mse_se <- sqrt(rowMeans(((true_value-data[,1,])^2-mse)^2)/(200*199))
    #relative bias
    relative_bias <- rowMeans((data[,1,]-true_value)/true_value)
    rb_se <- sqrt(rowSums(((data[,1,]-true_value)/true_value-relative_bias)^2)/(200*199))
    #power
    tpr <- tnr <- array(dim=c(76,200))
    true_positive <- true_value!=0
    tpr[data[,5,]==true_positive] <- 1
    tpr[is.na(tpr)] <- 0
    tpr[!true_positive,] <- NA
    tnr[data[,5,]==true_positive] <- 1
    tnr[is.na(tnr)] <- 0
    tnr[true_positive,] <- NA
    rownames(tpr) <- rn
    true_positive_rate <- rowMeans(tpr)
    true_negative_rate <- rowMeans(tnr)
    tpr_se <- sqrt(true_positive_rate*(1-true_positive_rate)/200)
    tnr_se <- sqrt(true_negative_rate*(1-true_negative_rate)/200)
    #总结果
    results <- cbind(cover_rate,mse,relative_bias,true_positive_rate,true_negative_rate)
    results_se <- cbind(cr_se,mse_se,rb_se,tpr_se,tnr_se)
    #输出
    out_name <- paste0('F:/course/statistic/tvm/S1/decay/est_results_n',
                       as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
    f1 <- paste0('F:/course/statistic/tvm/S1/decay/est_se_n',
                 as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results_se,file=f1,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}

foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    data<-array(dim=c(76,5,200))
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data[,,i] <- data1
    }
    rn <- rownames(data1)
    cn <- colnames(data1)
    dimnames(data)[[1]] <- rn
    dimnames(data)[[2]] <- cn
    #真值
    ##回归系数真值
    time <- c(1:9)*0.1
    beta <- rep(0.4,9)
    ##间接效应真值
    ind <- beta*beta
    true_value <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8,
                    rep(0.1,9),rep(0.4,9),beta,0.2,0.1,
                    rep(0.6,9),beta,beta,0.2,0.1,
                    ind)
    #区间覆盖率
    cr <- array(dim=c(76,200))
    cr[data[,3,]<=true_value&data[,4,]>=true_value] <- 1
    cr[is.na(cr)] <- 0
    rownames(cr) <- rn
    cover_rate <- rowMeans(cr)
    cr_se <- sqrt(cover_rate*(1-cover_rate)/200)
    #MSE
    mse <- rowMeans((true_value-data[,1,])^2)
    mse_se <- sqrt(rowMeans(((true_value-data[,1,])^2-mse)^2)/(200*199))
    #relative bias
    relative_bias <- rowMeans((data[,1,]-true_value)/true_value)
    rb_se <- sqrt(rowSums(((data[,1,]-true_value)/true_value-relative_bias)^2)/(200*199))
    #power
    tpr <- tnr <- array(dim=c(76,200))
    true_positive <- true_value!=0
    tpr[data[,5,]==true_positive] <- 1
    tpr[is.na(tpr)] <- 0
    tpr[!true_positive,] <- NA
    tnr[data[,5,]==true_positive] <- 1
    tnr[is.na(tnr)] <- 0
    tnr[true_positive,] <- NA
    rownames(tpr) <- rn
    true_positive_rate <- rowMeans(tpr)
    true_negative_rate <- rowMeans(tnr)
    tpr_se <- sqrt(true_positive_rate*(1-true_positive_rate)/200)
    tnr_se <- sqrt(true_negative_rate*(1-true_negative_rate)/200)
    #总结果
    results <- cbind(cover_rate,mse,relative_bias,true_positive_rate,true_negative_rate)
    results_se <- cbind(cr_se,mse_se,rb_se,tpr_se,tnr_se)
    #输出
    out_name <- paste0('F:/course/statistic/tvm/S1/flat/est_results_n',
                       as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
    f1 <- paste0('F:/course/statistic/tvm/S1/flat/est_se_n',
                 as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results_se,file=f1,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}
foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    data<-array(dim=c(76,5,200))
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data[,,i] <- data1
    }
    rn <- rownames(data1)
    cn <- colnames(data1)
    dimnames(data)[[1]] <- rn
    dimnames(data)[[2]] <- cn
    #真值
    ##回归系数真值
    time <- c(1:9)*0.1
    beta <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
    ##间接效应真值
    ind <- beta*beta
    true_value <- c(1,0.8,0.8,1,0.8,0.8,1,0.8,0.8,
                    rep(0.1,9),rep(0.4,9),beta,0.2,0.1,
                    rep(0.6,9),beta,beta,0.2,0.1,
                    ind)
    #区间覆盖率
    cr <- array(dim=c(76,200))
    cr[data[,3,]<=true_value&data[,4,]>=true_value] <- 1
    cr[is.na(cr)] <- 0
    rownames(cr) <- rn
    cover_rate <- rowMeans(cr)
    cr_se <- sqrt(cover_rate*(1-cover_rate)/200)
    #MSE
    mse <- rowMeans((true_value-data[,1,])^2)
    mse_se <- sqrt(rowMeans(((true_value-data[,1,])^2-mse)^2)/(200*199))
    #relative bias
    relative_bias <- rowMeans((data[,1,]-true_value)/true_value)
    rb_se <- sqrt(rowSums(((data[,1,]-true_value)/true_value-relative_bias)^2)/(200*199))
    #power
    tpr <- tnr <- array(dim=c(76,200))
    true_positive <- true_value!=0
    tpr[data[,5,]==true_positive] <- 1
    tpr[is.na(tpr)] <- 0
    tpr[!true_positive,] <- NA
    tnr[data[,5,]==true_positive] <- 1
    tnr[is.na(tnr)] <- 0
    tnr[true_positive,] <- NA
    rownames(tpr) <- rn
    true_positive_rate <- rowMeans(tpr)
    true_negative_rate <- rowMeans(tnr)
    tpr_se <- sqrt(true_positive_rate*(1-true_positive_rate)/200)
    tnr_se <- sqrt(true_negative_rate*(1-true_negative_rate)/200)
    #总结果
    results <- cbind(cover_rate,mse,relative_bias,true_positive_rate,true_negative_rate)
    results_se <- cbind(cr_se,mse_se,rb_se,tpr_se,tnr_se)
    #输出
    out_name <- paste0('F:/course/statistic/tvm/S1/hump/est_results_n',
                       as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results,file=out_name,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
    f1 <- paste0('F:/course/statistic/tvm/S1/hump/est_se_n',
                 as.character(n[j]),'m',as.character(m[k]),'.txt')
    write.table(results_se,file=f1,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
  }
}


type = c('flat','decay','hump','valley')
para = c('loadings','measurement error','intercept','path coefficients',
         'regression residual','indirect effect')
data2 <- data3 <- array(dim=c(192,9))
colnames(data2) <- colnames(data3) <- c('type','sample size','measurement occasion',
                                        'coefficient','cover rate','MSE',
                                        'relative bias','true positive rate',
                                        'true negative rate')
all_cr <- all_rb <- all_mse <- all_tpr <- all_tnr <- all_cover_rate <- 
  all_MeanSE <- all_relative_bias <- all_TruePR <- all_TrueNR <- array(dim=c(4,4,2,6))
data2[,1] <- data3[,1] <- rep(type,each=48)
data2[,2] <- data3[,2] <- rep(rep(n,each=12),4)
data2[,3] <- data3[,3] <- rep(rep(m,each=6),16)
data2[,4] <- data3[,4] <- rep(para,32)
foreach (l=1:4,.errorhandling="remove") %do% {
  foreach (j=1:4,.errorhandling="remove") %do% {
    foreach (k=1:2,.errorhandling="remove") %do% {
      f1 <- paste0('F:/course/statistic/tvm/S1/',type[l],'/est_results_n',
                   as.character(n[j]),'m',as.character(m[k]),'.txt')
      data0 <- read.csv(file=f1,header=TRUE,sep='\t')
      tf <- data2[,1]==type[l]&data2[,2]==n[j]&data2[,3]==m[k]
      cr_loadings <- mean(data0$cover_rate[c(2:3,5:6,8:9)])
      cr_me <- mean(data0$cover_rate[c(10:18)])
      cr_int <- mean(data0$cover_rate[c(19:27,39:47)])
      cr_path <- mean(data0$cover_rate[c(28:36,48:65)])
      cr_re <- mean(data0$cover_rate[c(37:38,66:67)])
      cr_ind <- mean(data0$cover_rate[c(68:76)])
      data2[tf,5] <- all_cr[l,j,k,] <- c(cr_loadings,cr_me,cr_int,cr_path,cr_re,cr_ind)
      mse_loadings <- mean(data0$mse[c(2:3,5:6,8:9)])
      mse_me <- mean(data0$mse[c(10:18)])
      mse_int <- mean(data0$mse[c(19:27,39:47)])
      mse_path <- mean(data0$mse[c(28:36,48:65)])
      mse_re <- mean(data0$mse[c(37:38,66:67)])
      mse_ind <- mean(data0$mse[c(68:76)])
      data2[tf,6] <- all_mse[l,j,k,] <- c(mse_loadings,mse_me,mse_int,mse_path,mse_re,mse_ind)
      rb_loadings <- mean(data0$relative_bias[c(2:3,5:6,8:9)])
      rb_me <- mean(data0$relative_bias[c(10:18)])
      rb_int <- mean(data0$relative_bias[c(19:27,39:47)])
      rb_path <- mean(na.omit(data0$relative_bias[c(28:36,48:65)]))
      rb_re <- mean(data0$relative_bias[c(37:38,66:67)])
      rb_ind <- mean(na.omit(data0$relative_bias[c(68:76)]))
      data2[tf,7] <- all_rb[l,j,k,] <- c(rb_loadings,rb_me,rb_int,rb_path,rb_re,rb_ind)
      tpr_loadings <- mean(data0$true_positive_rate[c(2:3,5:6,8:9)])
      tpr_me <- mean(data0$true_positive_rate[c(10:18)])
      tpr_int <- mean(data0$true_positive_rate[c(19:27,39:47)])
      tpr_path <- mean(na.omit(data0$true_positive_rate[c(28:36,48:65)]))
      tpr_re <- mean(data0$true_positive_rate[c(37:38,66:67)])
      tpr_ind <- mean(na.omit(data0$true_positive_rate[c(68:76)]))
      data2[tf,8] <- all_tpr[l,j,k,] <- c(tpr_loadings,tpr_me,tpr_int,tpr_path,tpr_re,tpr_ind)
      tnr_loadings <- mean(na.omit(data0$true_negative_rate[c(2:3,5:6,8:9)]))
      tnr_me <- mean(na.omit(data0$true_negative_rate[c(10:18)]))
      tnr_int <- mean(na.omit(data0$true_negative_rate[c(19:27,39:47)]))
      tnr_path <- mean(na.omit(data0$true_negative_rate[c(28:36,48:65)]))
      tnr_re <- mean(na.omit(data0$true_negative_rate[c(37:38,66:67)]))
      tnr_ind <- mean(na.omit(data0$true_negative_rate[c(68:76)]))
      data2[tf,9] <- all_tnr[l,j,k,] <- c(tnr_loadings,tnr_me,tnr_int,tnr_path,tnr_re,tnr_ind)
      
      f2 <- paste0('F:/course/statistic/tvm/S1/',type[l],'/est_se_n',
                   as.character(n[j]),'m',as.character(m[k]),'.txt')
      data1 <- read.csv(file=f2,header=TRUE,sep='\t')
      cr.se_loadings <- mean(data1$cr_se[c(2:3,5:6,8:9)])
      cr.se_me <- mean(data1$cr_se[c(10:18)])
      cr.se_int <- mean(data1$cr_se[c(19:27,39:47)])
      cr.se_path <- mean(data1$cr_se[c(28:36,48:65)])
      cr.se_re <- mean(data1$cr_se[c(37:38,66:67)])
      cr.se_ind <- mean(data1$cr_se[c(68:76)])
      cr.se <- c(cr.se_loadings,cr.se_me,cr.se_int,cr.se_path,cr.se_re,cr.se_ind)
      for (i in 1:6) {
        all_cover_rate[l,j,k,i] <- paste0(sprintf('%0.3f',all_cr[l,j,k,i]),'(',
                                          sprintf('%0.3f',cr.se[i]),')')
      }
      data3[tf,5] <- all_cover_rate[l,j,k,]
      mse.se_loadings <- mean(data1$mse_se[c(2:3,5:6,8:9)])
      mse.se_me <- mean(data1$mse_se[c(10:18)])
      mse.se_int <- mean(data1$mse_se[c(19:27,39:47)])
      mse.se_path <- mean(data1$mse_se[c(28:36,48:65)])
      mse.se_re <- mean(data1$mse_se[c(37:38,66:67)])
      mse.se_ind <- mean(data1$mse_se[c(68:76)])
      mse.se <- c(mse.se_loadings,mse.se_me,mse.se_int,mse.se_path,mse.se_re,mse.se_ind)
      for (i in 1:6) {
        all_MeanSE[l,j,k,i] <- paste0(sprintf('%0.8f',all_mse[l,j,k,i]),'(',
                                      sprintf('%0.8f',mse.se[i]),')')
      }
      data3[tf,6] <- all_MeanSE[l,j,k,]
      rb.se_loadings <- mean(data1$rb_se[c(2:3,5:6,8:9)])
      rb.se_me <- mean(data1$rb_se[c(10:18)])
      rb.se_int <- mean(data1$rb_se[c(19:27,39:47)])
      rb.se_path <- mean(na.omit(data1$rb_se[c(28:36,48:65)]))
      rb.se_re <- mean(data1$rb_se[c(37:38,66:67)])
      rb.se_ind <- mean(na.omit(data1$rb_se[c(68:76)]))
      rb.se <- c(rb.se_loadings,rb.se_me,rb.se_int,rb.se_path,rb.se_re,rb.se_ind)
      for (i in 1:6) {
        all_relative_bias[l,j,k,i] <- paste0(sprintf('%0.3f',all_rb[l,j,k,i]),'(',
                                             sprintf('%0.3f',rb.se[i]),')')
      }
      data3[tf,7] <- all_relative_bias[l,j,k,]
      tpr.se_loadings <- mean(data1$tpr_se[c(2:3,5:6,8:9)])
      tpr.se_me <- mean(data1$tpr_se[c(10:18)])
      tpr.se_int <- mean(data1$tpr_se[c(19:27,39:47)])
      tpr.se_path <- mean(na.omit(data1$tpr_se[c(28:36,48:65)]))
      tpr.se_re <- mean(data1$tpr_se[c(37:38,66:67)])
      tpr.se_ind <- mean(na.omit(data1$tpr_se[c(68:76)]))
      tpr.se <- c(tpr.se_loadings,tpr.se_me,tpr.se_int,tpr.se_path,tpr.se_re,tpr.se_ind)
      for (i in 1:6) {
        all_TruePR[l,j,k,i] <- paste0(sprintf('%0.3f',all_tpr[l,j,k,i]),'(',
                                      sprintf('%0.3f',tpr.se[i]),')')
      }
      data3[tf,8] <- all_TruePR[l,j,k,]
      tnr.se_loadings <- mean(na.omit(data0$tnr_se[c(2:3,5:6,8:9)]))
      tnr.se_me <- mean(na.omit(data0$tnr_se[c(10:18)]))
      tnr.se_int <- mean(na.omit(data0$tnr_se[c(19:27,39:47)]))
      tnr.se_path <- mean(na.omit(data0$tnr_se[c(28:36,48:65)]))
      tnr.se_re <- mean(na.omit(data0$tnr_se[c(37:38,66:67)]))
      tnr.se_ind <- mean(na.omit(data0$tnr_se[c(68:76)]))
      all_tnr.se[l,j,k,] <- c(tnr.se_loadings,tnr.se_me,tnr.se_int,tnr.se_path,tnr.se_re,tnr.se_ind)
      for (i in 1:6) {
        all_TrueNR[l,j,k,i] <- paste0(sprintf('%0.3f',all_tnr[l,j,k,i]),'(',
                                      sprintf('%0.3f',tnr.se[i]),')')
      }
      data3[tf,9] <- all_TrueNR[l,j,k,]
    }
  }
}
write.table(data2,file='F:/course/statistic/tvm/S1/results.txt',sep='\t',quote=FALSE,
            col.names = TRUE,row.names = FALSE)
write.table(data3,file='F:/course/statistic/tvm/S1/results_se.txt',sep='\t',quote=FALSE,
            col.names = TRUE,row.names = FALSE)
