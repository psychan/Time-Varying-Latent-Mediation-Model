wd='F:/course/statistic/tvm/S6/s6a'
setwd(wd)
source('F:/course/statistic/tvm/tvmediation_noi.R')


data0 <- read.csv(file='data0.csv',header=TRUE)
data1 <- na.omit(data0)

data1$survey <- data1$survey/50
for (i in 6:17) {
  data1[,i] <- scale(data1[,i])
}

write.table(data1[,c(1:2,5:17)],file='data1.txt',sep='\t',quote=FALSE,row.names = FALSE,col.names = FALSE)


results <- tvmediation(data=data1,x=data1[,6:8],m=data1[,9:11],
                       y=data1[,12:17],cluster = data1[,1],t=data1[,2],
                       nkot=9,iter=10000,burnin=5000);
f1 <- paste0('raw_ed.txt')
f2 <- paste0('remain_ed.txt')
f3 <- paste0('est_ed.txt')
write.table(results$raw,file=f1,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
write.table(results$remain,file=f2,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
write.table(results$est,file=f3,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)

data1[,6:17] <- scale(data1[,6:17])
for (chain in 1:2) {
  results <- tvmediation(data=data1,x=cbind(data1[,7:8],data1[,6]),
                         m=cbind(data1[,10:11],data1[,9]),
                         y=data1[,12:17],
                         cluster = data1[,1],t=data1[,2],
                         nkot=9,iter=10000,burnin=5000);
  f1 <- paste0('raw_ed_',chain,'.txt')
  f2 <- paste0('remain_ed_',chain,'.txt')
  f3 <- paste0('est_ed_',chain,'.txt')
  write.table(results$raw,file=f1,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
  write.table(results$remain,file=f2,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
  write.table(results$est,file=f3,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
}


est_ed<- read.csv(file='F:/course/statistic/tvm/S6/s6a/est_ed_1.txt',header=TRUE,sep='\t')
t <- c(1:9)*0.1
ind_m <- est_ed[74:82,1]
ind_l <- est_ed[74:82,3]
ind_u <- est_ed[74:82,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ed.png')
png(f_name)
plot(t,ind_m,xlim=c(0,1),ylim=c(0,0.5),"l", xlab="time",
     ylab = "indirect effect",col='#00bfff',lwd=2)
lines(t,ind_l,col='#00bfff',lty=2,lwd=2)
lines(t,ind_u,col='#00bfff',lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

a_m <- est_ed[34:42,1]
a_l <- est_ed[34:42,3]
a_u <- est_ed[34:42,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ed_alpha.png')
png(f_name)
plot(t,a_m,xlim=c(0,1),ylim=c(0,0.6),"l", xlab="time",
     ylab = "religiousness → spirituality",col="#FF5733",lwd=2)
lines(t,a_l,col="#FF5733",lty=2,lwd=2)
lines(t,a_u,col="#FF5733",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

b_m <- est_ed[54:62,1]
b_l <- est_ed[54:62,3]
b_u <- est_ed[54:62,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ed_beta.png')
png(f_name)
plot(t,b_m,xlim=c(0,1),ylim=c(0,0.6),"l", xlab="time",
     ylab = "spirituality → moral affect",col="#E41A1C",lwd=2)
lines(t,b_l,col="#E41A1C",lty=2,lwd=2)
lines(t,b_u,col="#E41A1C",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

g_m <- est_ed[63:71,1]
g_l <- est_ed[63:71,3]
g_u <- est_ed[63:71,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ed_gamma.png')
png(f_name)
plot(t,g_m,xlim=c(0,1),ylim=c(0,0.6),"l", xlab="time",
     ylab = "religiousness → moral affect",col="#8f001a",lwd=2)
lines(t,g_l,col="#8f001a",lty=2,lwd=2)
lines(t,g_u,col="#8f001a",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

nid <- unique(data1[,1])
sub <- as.numeric(match(data1[,1], nid))
n_obs <- as.numeric(as.vector(table(data1[,1])))
age <- gender <- matrix(0,length(nid),1)
ai <- 1
for (i in nid) {
  age[ai] <- unique(data1[data1[,1]==i,5])
  gender[ai] <- unique(data1[data1[,1]==i,22])
  ai <- ai+1
}
mean(age)
sd(age)
n_sex <- as.numeric(as.vector(table(gender)))


source('F:/course/statistic/tvm/DynMed.R')
results <- DynMed(data=data1,x=rowMeans(data1[,6:8]),
                  m=rowMeans(data1[,9:11]),y=rowMeans(data1[,12:17]),
                  cluster = data1[,1],t=data1[,2])
f1 <- paste0(wd,'raw_nocfa.txt')
f2 <- paste0(wd,'remain_nocfa.txt')
f3 <- paste0(wd,'est_nocfa.txt')
write.table(results$raw,file=f1,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
write.table(results$remain,file=f2,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
write.table(results$est,file=f3,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)

est_ncfa<- read.csv(file='F:/course/statistic/tvm/S6/s6a/est_nocfa.txt',header=TRUE,sep='\t')
t <- c(1:9)*0.1
ind_m <- est_ncfa[50:58,1]
ind_l <- est_ncfa[50:58,3]
ind_u <- est_ncfa[50:58,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ncfa.png')
png(f_name)
plot(t,ind_m,xlim=c(0,1),ylim=c(0,0.5),"l", xlab="time",
     ylab = "indirect effect",col='#00bfff',lwd=2)
lines(t,ind_l,col='#00bfff',lty=2,lwd=2)
lines(t,ind_u,col='#00bfff',lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

a_m <- est_ncfa[10:18,1]
a_l <- est_ncfa[10:18,3]
a_u <- est_ncfa[10:18,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ncfa_alpha.png')
png(f_name)
plot(t,a_m,xlim=c(0,1),ylim=c(0,0.6),"l", xlab="time",
     ylab = "religiousness → spirituality",col="#FF5733",lwd=2)
lines(t,a_l,col="#FF5733",lty=2,lwd=2)
lines(t,a_u,col="#FF5733",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

b_m <- est_ncfa[30:38,1]
b_l <- est_ncfa[30:38,3]
b_u <- est_ncfa[30:38,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ncfa_beta.png')
png(f_name)
plot(t,b_m,xlim=c(0,1),ylim=c(0,0.6),"l", xlab="time",
     ylab = "spirituality → moral affect",col="#E41A1C",lwd=2)
lines(t,b_l,col="#E41A1C",lty=2,lwd=2)
lines(t,b_u,col="#E41A1C",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

g_m <- est_ncfa[39:47,1]
g_l <- est_ncfa[39:47,3]
g_u <- est_ncfa[39:47,4]
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_ncfa_gamma.png')
png(f_name)
plot(t,g_m,xlim=c(0,1),ylim=c(-0.1,0.6),"l", xlab="time",
     ylab = "religiousness → moral affect",col="#8f001a",lwd=2)
lines(t,g_l,col="#8f001a",lty=2,lwd=2)
lines(t,g_u,col="#8f001a",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

#mlm
t <- c(1:9)*0.1
ind_m <- rep(0.717,9)
ind_l <- rep(0.523,9)
ind_u <- rep(0.956,9)
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_mlm.png')
png(f_name)
plot(t,ind_m,xlim=c(0,1),ylim=c(0,1),"l", xlab="time",
     ylab = "indirect effect",col='#00bfff',lwd=2)
lines(t,ind_l,col='#00bfff',lty=2,lwd=2)
lines(t,ind_u,col='#00bfff',lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

a_m <- rep(0.729,9)
a_l <- rep(0.624,9)
a_u <- rep(0.846,9)
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_mlm_alpha.png')
png(f_name)
plot(t,a_m,xlim=c(0,1),ylim=c(0,1),"l", xlab="time",
     ylab = "religiousness → spirituality",col="#FF5733",lwd=2)
lines(t,a_l,col="#FF5733",lty=2,lwd=2)
lines(t,a_u,col="#FF5733",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

b_m <- rep(0.982,9)
b_l <- rep(0.769,9)
b_u <- rep(1.216,9)
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_mlm_beta.png')
png(f_name)
plot(t,b_m,xlim=c(0,1),ylim=c(0.5,1.5),"l", xlab="time",
     ylab = "spirituality → moral affect",col="#E41A1C",lwd=2)
lines(t,b_l,col="#E41A1C",lty=2,lwd=2)
lines(t,b_u,col="#E41A1C",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()

g_m <- rep(-0.063,9)
g_l <- rep(-0.296,9)
g_u <- rep(0.139,9)
f_name <- paste0('F:/course/statistic/tvm/S6/s6a/est_mlm_gamma.png')
png(f_name)
plot(t,g_m,xlim=c(0,1),ylim=c(-0.5,0.5),"l", xlab="time",
     ylab = "religiousness → moral affect",col="#8f001a",lwd=2)
lines(t,g_l,col="#8f001a",lty=2,lwd=2)
lines(t,g_u,col="#8f001a",lty=2,lwd=2)
abline(h=0,col='#000000')
dev.off()
