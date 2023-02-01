tidy.DMN <- function(x, ...) {
  ret <- as.data.frame(x@fit)
  tbl_df(fix_data_frame(ret, c("conf.low", "estimate", "conf.high")))}

ddiric <- function(x,attr2=attr,sim2=sim) {
  
  cer2 <- x
  final <- matrix(0,nrow(cer2),ncol(cer2))
  colnames(final) <- colnames(cer2)
  row.names(final) <- row.names(cer2)
  ddist <- as.matrix(dist(cbind(attr2$EASTING,attr2$NORTHING)))
  ddist <- max(ddist)-ddist
  row.names(ddist) <- row.names(cer2)
  
  ddist <- (ddist-min(ddist))/(max(ddist)-min(ddist))
  mult1 <- matrix(NA,nrow(ddist),ncol(ddist))
  for (i in 1:length(sim2)) {
    mult1[i] <- weighted.mean(x=c(sim2[i],ddist[i]),w=c(0.5,0.5))}
  diag(mult1) <- 0
  
  for (m in 1:nrow(cer2)) {
    cer.a <- as.matrix(cer2)
    cer.a <- cer2[which(mult1[m,]>quantile(mult1[m,],probs=c(0.5))),]
    cer.ret <- which(colSums(cer.a)>5)
    cer.a <- round(cer.a[,cer.ret],0)
    cer.ret2 <- which(rowSums(cer.a)>0)
    cer.a <- cer.a[cer.ret2,]
    cer.b <- as.matrix(cer.a)
    row.names(cer.b) <- row.names(cer.a)
    colnames(cer.b) <- colnames(cer.a)
    cer.a <- cer.b
    
    dm_fit <- DirichletMultinomial::dmn(cer.a,1)
    
    
    
    dm_params <- tidy(dm_fit)
    
    par.e <- dm_params$estimate
    par_total <- sum(par.e)
    
    item_mean <- dm_params$estimate / sum(dm_params$estimate)
    
    cer.p <- prop.table(cer.a,1)
    out.type <- cer.a
    for (i in 1:ncol(cer.a)) {
      out.type[,i] <- (cer.a[,i]+par.e[i])/(rowSums(cer.a)+par_total)}
    
    site <- t(as.matrix(as.matrix(cer)[m,cer.ret]))
    out.site <- site
    out.site <- prop.table(as.matrix(out.site),1)
    for (i in 1:ncol(site)) {
      out.site[1,i] <- (site[1,i]+par.e[i])/(rowSums(site)+par_total)}
    row.names(out.site) <- v1 <- row.names(cer)[m]
    colnames(out.site) <- v2 <- colnames(cer)[cer.ret]
    
    for(i in 1:length(cer.ret)) {
      final[m,cer.ret[i]] <- out.site[1,i]}
  }
  
  logcer <- log(rowSums(cer))
  logab <- log(rowSums(ab))
  
  fit <- lm(logcer~logab)
  
  outlier <- names(outlierTest(fit)$rstudent)
  temp <- NULL
  for (i in 1:length(outlier)) {
    temp[i] <- which(row.names(cer)==outlier[i])}
  
  outlier <- temp
  
  
  final2 <- as.data.frame(final)
  for (i in outlier) {
    final2[i,] <- cer[i,]/rowSums(cer[i,])}
  
  return(final2)
}


simulation2 <- function(orig100,frac) {
  out.list <- list()
  nsim <- 250 # set number of simulations
  
  for (j in 1:nsim) {
    data.sim <- NULL
    # this for loop creates 1 random replicate of the underlying ceramic data 
    for (m in 1:nrow(orig100)) {
      data.sim <- rbind(data.sim,t(rmultinom(1,rowSums(orig100)[m]*frac,prob=orig100[m,])))}
    out.list[[j]] <- data.sim}
  return(out.list)}


simulation.d <- function(diric,orig100) {
  out.list <- list()
  nsim <- 250 # set number of simulations
  
  for (j in 1:nsim) {
    data.sim <- NULL
    # this for loop creates 1 random replicate of the underlying ceramic data 
    for (m in 1:nrow(diric)) {
      data.sim <- rbind(data.sim,t(rmultinom(1,rowSums(orig100)[m],prob=diric[m,])))}
    out.list[[j]] <- data.sim}
  return(out.list)}


simulation.o <- function(orig100) {
  out.list <- list()
  nsim <- 250 # set number of simulations
  
  for (j in 1:nsim) {
    data.sim <- NULL
    # this for loop creates 1 random replicate of the underlying ceramic data 
    for (m in 1:nrow(orig100)) {
      data.sim <- rbind(data.sim,t(rmultinom(1,rowSums(orig100)[m]*sample(seq(0.25,0.75,by=0.01),1),prob=orig100[m,])))}
    out.list[[j]] <- data.sim}
  return(out.list)}


for(i in 1:250) {vv[[i]] <- sim.mat(test2[[i]])}
vva <- matrix(NA,222,250)
for(i in 1:250) {vva[,i] <- rowSums(vv[[i]])}

zz <- matrix(NA,nsim,4)
for (i in 1:nsim) {
  zz[i,1] <- mean((sim.mat(test2[[i]])[which(z==1),]-sim[which(z==1),])^2)
  zz[i,2] <- mean((sim.mat(test2[[i]])[which(z==10),]-sim[which(z==10),])^2)
  zz[i,3] <- mean((sim.mat(test2[[i]])[which(z==100),]-sim[which(z==100),])^2)
  zz[i,4] <- mean((sim.mat(test2[[i]])[which(z==1000),]-sim[which(z==1000),])^2)}

zz <- matrix(NA,nsim,4)
for (i in 1:nsim) {
  zz[i,1] <- mean((sim.mat(test2[[i]])[which(z==1),]-sim[which(z==1),])^2)
  zz[i,2] <- mean((sim.mat(test2[[i]])[which(z==10),]-sim[which(z==10),])^2)
  zz[i,3] <- mean((sim.mat(test2[[i]])[which(z==100),]-sim[which(z==100),])^2)
  zz[i,4] <- mean((sim.mat(test2[[i]])[which(z==1000),]-sim[which(z==1000),])^2)}


bins2 <- c(1,10,100,1000)
plot(log10(bins2),(zz2[2,]),type='l',col='blue',lwd=2,ylim=c(0,0.03),xlab='Log Sample Size',ylab='Mean Squared Error')
points(log10(bins2),(zz2.orig[2,]),type='l',col='red',lwd=2)
legend('topright',c('Bootstrap Resampling','Empirical Bayesian Estimation'),cex=0.75,bty='n',lwd=c(2,2),col=c('red','blue'))

mat2 <- matrix(NA,222,222)
row.names(mat2) <- row.names(cer)
colnames(mat2) <- row.names(cer)
for(i in 1:222) {
  for(j in 1:222) {
    mat2[i,j] <- t.test(vva[i,],vva[j,])$p.value
  }
}

mat3 <- matrix(NA,222,222)
row.names(mat3) <- row.names(cer)
colnames(mat3) <- row.names(cer)
for(i in 1:222) {
  for(j in 1:222) {
    mat3[i,j] <- mean(abs(vva[i,]-vva[j,]))
  }
}


polygon(c(log10(bins2),log10(rev(bins2))),c((zz2[1,]),rev((zz2[3,]))),border=0,col=rgb(0,0,1,0.25))
polygon(c(log10(bins2),log10(rev(bins2))),c((zz2.orig[1,]),rev((zz2.orig[3,]))),border=0,col=rgb(1,0,0,0.25))


orig.sample <- simulation.o(orig100)

orig.list90 <- simulation.o(orig100,frac=0.9)
orig.list80 <- simulation.o(orig100,frac=0.8)
orig.list70 <- simulation.o(orig100,frac=0.7)
orig.list60 <- simulation.o(orig100,frac=0.6)
orig.list50 <- simulation.o(orig100,frac=0.5)
orig.list40 <- simulation.o(orig100,frac=0.4)
orig.list30 <- simulation.o(orig100,frac=0.3)
orig.list20 <- simulation.o(orig100,frac=0.2)

out.list90 <- simulation2(orig100,frac=90)
out.list80 <- simulation2(orig100,frac=80)
out.list70 <- simulation2(orig100,frac=70)
out.list60 <- simulation2(orig100,frac=60)
out.list50 <- simulation2(orig100,frac=50)
out.list40 <- simulation2(orig100,frac=40)
out.list30 <- simulation2(orig100,frac=30)
out.list20 <- simulation2(orig100,frac=20)

diric100 <- ddiric(orig100,attr,sim)
diric90 <- ddiric(orig90,attr,sim)
diric80 <- ddiric(orig80,attr,sim)
diric70 <- ddiric(orig70,attr,sim)
diric60 <- ddiric(orig60,attr,sim)
diric50 <- ddiric(orig50,attr,sim)
diric40 <- ddiric(orig40,attr,sim)
diric30 <- ddiric(orig30,attr,sim)
diric20 <- ddiric(orig20,attr,sim)

diric.sim <- simulation.d(diric100,orig100)
diric.rs90 <- simulation.d(diric90,orig100)
diric.rs80 <- simulation.d(diric80,orig100)
diric.rs70 <- simulation.d(diric70,orig100)
diric.rs60 <- simulation.d(diric60,orig100)
diric.rs50 <- simulation.d(diric50,orig100)
diric.rs40 <- simulation.d(diric40,orig100)
diric.rs30 <- simulation.d(diric30,orig100)
diric.rs20 <- simulation.d(diric20,orig100)

diric2.rs.samp <- simulation.d(diric100,orig100)

diric.2rs90 <- simulation.d(diric90,orig90)
diric.2rs80 <- simulation.d(diric80,orig80)
diric.2rs70 <- simulation.d(diric70,orig70)
diric.2rs60 <- simulation.d(diric60,orig60)
diric.2rs50 <- simulation.d(diric50,orig50)
diric.2rs40 <- simulation.d(diric40,orig40)
diric.2rs30 <- simulation.d(diric30,orig30)
diric.2rs20 <- simulation.d(diric20,orig20)

diric.list100b <- list()
for(i in 1:nsim) {diric.list100b[[i]] <- ddiric(orig3[[i]],attr,sim)}


for(i in 1:nsim) {diric.list90[[i]] <- ddiric(out.list90[[i]],attr,sim)}
for(i in 1:nsim) {diric.list80[[i]] <- ddiric(out.list80[[i]],attr,sim)}
for(i in 1:nsim) {diric.list70[[i]] <- ddiric(out.list70[[i]],attr,sim)}
for(i in 1:nsim) {diric.list60[[i]] <- ddiric(out.list60[[i]],attr,sim)}
for(i in 1:nsim) {diric.list50[[i]] <- ddiric(out.list50[[i]],attr,sim)}
for(i in 1:nsim) {diric.list40[[i]] <- ddiric(out.list40[[i]],attr,sim)}
for(i in 1:nsim) {diric.list30[[i]] <- ddiric(out.list30[[i]],attr,sim)}
for(i in 1:nsim) {diric.list20[[i]] <- ddiric(out.list20[[i]],attr,sim)}

diric.sample <- list()
for(i in 1:nsim) {diric.sample[[i]] <- ddiric(orig.sample[[i]],attr,sim)}


bayes2 <- function (x,cer.prop) {
  bayes <- NULL
  bayes.st <- matrix(NA,nrow(cer.prop),nsim)
  for(i in 1:nsim) {
    bayes[i] <- mean((prop.table(as.matrix(x[[i]]),1)-cer.prop)^2)
    bayes.st[,i] <- apply(((prop.table(as.matrix(x[[i]]),1)-cer.prop)^2),1,mean)}
  v.out <- list()
  v.out[[1]] <- bayes
  v.out[[2]] <- bayes.st
  return(v.out)}

bayes.rs90 <- bayes2(diric.rs90,cer.prop)
bayes.rs80 <- bayes2(diric.rs80,cer.prop)
bayes.rs70 <- bayes2(diric.rs70,cer.prop)
bayes.rs60 <- bayes2(diric.rs60,cer.prop)
bayes.rs50 <- bayes2(diric.rs50,cer.prop)
bayes.rs40 <- bayes2(diric.rs40,cer.prop)
bayes.rs30 <- bayes2(diric.rs30,cer.prop)
bayes.rs20 <- bayes2(diric.rs20,cer.prop)


bayes.2rs90 <- bayes2(diric.2rs90,cer.prop)
bayes.2rs80 <- bayes2(diric.2rs80,cer.prop)
bayes.2rs70 <- bayes2(diric.2rs70,cer.prop)
bayes.2rs60 <- bayes2(diric.2rs60,cer.prop)
bayes.2rs50 <- bayes2(diric.2rs50,cer.prop)
bayes.2rs40 <- bayes2(diric.2rs40,cer.prop)
bayes.2rs30 <- bayes2(diric.2rs30,cer.prop)
bayes.2rs20 <- bayes2(diric.2rs20,cer.prop)

test <- bayes2(diric.sample,cer.prop)

boot2 <- function(x,cer.prop) {
  boot <- NULL
  boot.st <- matrix(NA,nrow(cer),nsim)
  for(i in 1:nsim) {
    boot[i] <- mean((prop.table(as.matrix(x[[i]]),1)-cer.prop)^2)
    boot.st[,i] <- apply(((prop.table(as.matrix(x[[i]]),1)-cer.prop)^2),1,mean)}
  vv.out <- list()
  vv.out[[1]] <- boot
  vv.out[[2]] <- boot.st
  return(vv.out)}

rmse <- function(x,cer.prop) {
  boot <- NULL
  boot.st <- matrix(NA,nrow(cer),nsim)
  for(i in 1:nsim) {
    boot[i] <- mean((x[[i]]-cer.prop)^2)
    boot.st[,i] <- apply(((x[[i]]-cer.prop)^2),1,mean)}
  vv.out <- list()
  vv.out[[1]] <- sqrt(boot)
  vv.out[[2]] <- sqrt(boot.st)
  return(vv.out)}


rmse2 <- function(x,cer.prop) {
  boot <- NULL
  boot.st <- matrix(NA,nrow(cer),nsim)
  for(i in 1:nsim) {
    boot[i] <- mean((prop.table(as.matrix(x[[i]]),1)-cer.prop)^2)
    boot.st[,i] <- apply(((prop.table(as.matrix(x[[i]]),1)-cer.prop)^2),1,mean)}
  vv.out <- list()
  vv.out[[1]] <- sqrt(boot)
  vv.out[[2]] <- sqrt(boot.st)
  return(vv.out)}



boot.sim20 <- rmse(sim.b.20,sim)
boot.sim30 <- rmse(sim.b.30,sim)
boot.sim40 <- rmse(sim.b.40,sim)
boot.sim50 <- rmse(sim.b.50,sim)
boot.sim60 <- rmse(sim.b.60,sim)
boot.sim70 <- rmse(sim.b.70,sim)
boot.sim80 <- rmse(sim.b.80,sim)
boot.sim90 <- rmse(sim.b.90,sim)

bayes.sim20 <- rmse(sim.d2.20,sim)
bayes.sim30 <- rmse(sim.d2.30,sim)
bayes.sim40 <- rmse(sim.d2.40,sim)
bayes.sim50 <- rmse(sim.d2.50,sim)
bayes.sim60 <- rmse(sim.d2.60,sim)
bayes.sim70 <- rmse(sim.d2.70,sim)
bayes.sim80 <- rmse(sim.d2.80,sim)
bayes.sim90 <- rmse(sim.d2.90,sim)

bayes.cer20 <- rmse(diric.list20,cer.prop)
bayes.cer30 <- rmse(diric.list30,cer.prop)
bayes.cer40 <- rmse(diric.list40,cer.prop)
bayes.cer50 <- rmse(diric.list50,cer.prop)
bayes.cer60 <- rmse(diric.list60,cer.prop)
bayes.cer70 <- rmse(diric.list70,cer.prop)
bayes.cer80 <- rmse(diric.list80,cer.prop)
bayes.cer90 <- rmse(diric.list90,cer.prop)


boot.rs90 <- rmse2(out.list90,cer.prop)
boot.rs80 <- rmse2(out.list80,cer.prop)
boot.rs70 <- rmse2(out.list70,cer.prop)
boot.rs60 <- rmse2(out.list60,cer.prop)
boot.rs50 <- rmse2(out.list50,cer.prop)
boot.rs40 <- rmse2(out.list40,cer.prop)
boot.rs30 <- rmse2(out.list30,cer.prop)
boot.rs20 <- rmse2(out.list20,cer.prop)

orig.rs90 <- boot2(orig.list90,cer.prop)
orig.rs80 <- boot2(orig.list80,cer.prop)
orig.rs70 <- boot2(orig.list70,cer.prop)
orig.rs60 <- boot2(orig.list60,cer.prop)
orig.rs50 <- boot2(orig.list50,cer.prop)
orig.rs40 <- boot2(orig.list40,cer.prop)
orig.rs30 <- boot2(orig.list30,cer.prop)
orig.rs20 <- boot2(orig.list20,cer.prop)

brs20 <- quantile(bayes.rs20[[1]],probs=c(0.05,0.5,0.95))
brs30 <- quantile(bayes.rs30[[1]],probs=c(0.05,0.5,0.95))
brs40 <- quantile(bayes.rs40[[1]],probs=c(0.05,0.5,0.95))
brs50 <- quantile(bayes.rs50[[1]],probs=c(0.05,0.5,0.95))
brs60 <- quantile(bayes.rs60[[1]],probs=c(0.05,0.5,0.95))
brs70 <- quantile(bayes.rs70[[1]],probs=c(0.05,0.5,0.95))
brs80 <- quantile(bayes.rs80[[1]],probs=c(0.05,0.5,0.95))
brs90 <- quantile(bayes.rs90[[1]],probs=c(0.05,0.5,0.95))


bsim20 <- quantile(bayes.sim20[[1]],probs=c(0.1,0.5,0.9))
bsim30 <- quantile(bayes.sim30[[1]],probs=c(0.1,0.5,0.9))
bsim40 <- quantile(bayes.sim40[[1]],probs=c(0.1,0.5,0.9))
bsim50 <- quantile(bayes.sim50[[1]],probs=c(0.1,0.5,0.9))
bsim60 <- quantile(bayes.sim60[[1]],probs=c(0.1,0.5,0.9))
bsim70 <- quantile(bayes.sim70[[1]],probs=c(0.1,0.5,0.9))
bsim80 <- quantile(bayes.sim80[[1]],probs=c(0.1,0.5,0.9))
bsim90 <- quantile(bayes.sim90[[1]],probs=c(0.1,0.5,0.9))


btsim20 <- quantile(boot.sim20[[1]],probs=c(0.1,0.5,0.9))
btsim30 <- quantile(boot.sim30[[1]],probs=c(0.1,0.5,0.9))
btsim40 <- quantile(boot.sim40[[1]],probs=c(0.1,0.5,0.9))
btsim50 <- quantile(boot.sim50[[1]],probs=c(0.1,0.5,0.9))
btsim60 <- quantile(boot.sim60[[1]],probs=c(0.1,0.5,0.9))
btsim70 <- quantile(boot.sim70[[1]],probs=c(0.1,0.5,0.9))
btsim80 <- quantile(boot.sim80[[1]],probs=c(0.1,0.5,0.9))
btsim90 <- quantile(boot.sim90[[1]],probs=c(0.1,0.5,0.9))

ors20 <- quantile(orig.rs20[[1]],probs=c(0.05,0.5,0.95))
ors30 <- quantile(orig.rs30[[1]],probs=c(0.05,0.5,0.95))
ors40 <- quantile(orig.rs40[[1]],probs=c(0.05,0.5,0.95))
ors50 <- quantile(orig.rs50[[1]],probs=c(0.05,0.5,0.95))
ors60 <- quantile(orig.rs60[[1]],probs=c(0.05,0.5,0.95))
ors70 <- quantile(orig.rs70[[1]],probs=c(0.05,0.5,0.95))
ors80 <- quantile(orig.rs80[[1]],probs=c(0.05,0.5,0.95))
ors90 <- quantile(orig.rs90[[1]],probs=c(0.05,0.5,0.95))

b2rs20 <- quantile(bayes.2rs20[[1]],probs=c(0.1,0.5,0.9))
b2rs30 <- quantile(bayes.2rs30[[1]],probs=c(0.1,0.5,0.9))
b2rs40 <- quantile(bayes.2rs40[[1]],probs=c(0.1,0.5,0.9))
b2rs50 <- quantile(bayes.2rs50[[1]],probs=c(0.1,0.5,0.9))
b2rs60 <- quantile(bayes.2rs60[[1]],probs=c(0.1,0.5,0.9))
b2rs70 <- quantile(bayes.2rs70[[1]],probs=c(0.1,0.5,0.9))
b2rs80 <- quantile(bayes.2rs80[[1]],probs=c(0.1,0.5,0.9))
b2rs90 <- quantile(bayes.2rs90[[1]],probs=c(0.1,0.5,0.9))



b220 <- quantile(bayes.cer20[[1]],probs=c(0.05,0.5,0.95))
b230 <- quantile(bayes.cer30[[1]],probs=c(0.05,0.5,0.95))
b240 <- quantile(bayes.cer40[[1]],probs=c(0.05,0.5,0.95))
b250 <- quantile(bayes.cer50[[1]],probs=c(0.05,0.5,0.95))
b260 <- quantile(bayes.cer60[[1]],probs=c(0.05,0.5,0.95))
b270 <- quantile(bayes.cer70[[1]],probs=c(0.05,0.5,0.95))
b280 <- quantile(bayes.cer80[[1]],probs=c(0.05,0.5,0.95))
b290 <- quantile(bayes.cer90[[1]],probs=c(0.05,0.5,0.95))


btrs20 <- quantile(boot.rs20[[1]],probs=c(0.05,0.5,0.95))
btrs30 <- quantile(boot.rs30[[1]],probs=c(0.05,0.5,0.95))
btrs40 <- quantile(boot.rs40[[1]],probs=c(0.05,0.5,0.95))
btrs50 <- quantile(boot.rs50[[1]],probs=c(0.05,0.5,0.95))
btrs60 <- quantile(boot.rs60[[1]],probs=c(0.05,0.5,0.95))
btrs70 <- quantile(boot.rs70[[1]],probs=c(0.05,0.5,0.95))
btrs80 <- quantile(boot.rs80[[1]],probs=c(0.05,0.5,0.95))
btrs90 <- quantile(boot.rs90[[1]],probs=c(0.05,0.5,0.95))

sim.d2.90 <- list()
for (i in 1:nsim) {sim.d2.90[[i]] <- sim.mat(diric.2rs90[[i]])}
sim.d2.80 <- list()
for (i in 1:nsim) {sim.d2.80[[i]] <- sim.mat(diric.2rs80[[i]])}
sim.d2.70 <- list()
for (i in 1:nsim) {sim.d2.70[[i]] <- sim.mat(diric.2rs70[[i]])}
sim.d2.60 <- list()
for (i in 1:nsim) {sim.d2.60[[i]] <- sim.mat(diric.2rs60[[i]])}
sim.d2.50 <- list()
for (i in 1:nsim) {sim.d2.50[[i]] <- sim.mat(diric.2rs50[[i]])}
sim.d2.40 <- list()
for (i in 1:nsim) {sim.d2.40[[i]] <- sim.mat(diric.2rs40[[i]])}
sim.d2.30 <- list()
for (i in 1:nsim) {sim.d2.30[[i]] <- sim.mat(diric.2rs30[[i]])}
sim.d2.20 <- list()
for (i in 1:nsim) {sim.d2.20[[i]] <- sim.mat(diric.2rs20[[i]])}

for (i in 1:nsim) {cor90.2[i] <- cor(rowSums(sim.d2.90[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor80.2[i] <- cor(rowSums(sim.d2.80[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor70.2[i] <- cor(rowSums(sim.d2.70[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor60.2[i] <- cor(rowSums(sim.d2.60[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor50.2[i] <- cor(rowSums(sim.d2.50[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor40.2[i] <- cor(rowSums(sim.d2.40[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor30.2[i] <- cor(rowSums(sim.d2.30[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor20.2[i] <- cor(rowSums(sim.d2.20[[i]]),rowSums(sim))}

sim.b.90 <- list()
for (i in 1:nsim) {sim.b.90[[i]] <- sim.mat(out.list90[[i]])}
sim.b.80 <- list()
for (i in 1:nsim) {sim.b.80[[i]] <- sim.mat(out.list80[[i]])}
sim.b.70 <- list()
for (i in 1:nsim) {sim.b.70[[i]] <- sim.mat(out.list70[[i]])}
sim.b.60 <- list()
for (i in 1:nsim) {sim.b.60[[i]] <- sim.mat(out.list60[[i]])}
sim.b.50 <- list()
for (i in 1:nsim) {sim.b.50[[i]] <- sim.mat(out.list50[[i]])}
sim.b.40 <- list()
for (i in 1:nsim) {sim.b.40[[i]] <- sim.mat(out.list40[[i]])}
sim.b.30 <- list()
for (i in 1:nsim) {sim.b.30[[i]] <- sim.mat(out.list30[[i]])}
sim.b.20 <- list()
for (i in 1:nsim) {sim.b.20[[i]] <- sim.mat(out.list20[[i]])}

for (i in 1:nsim) {cor90b[i] <- cor(rowSums(sim.b.90[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor80b[i] <- cor(rowSums(sim.b.80[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor70b[i] <- cor(rowSums(sim.b.70[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor60b[i] <- cor(rowSums(sim.b.60[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor50b[i] <- cor(rowSums(sim.b.50[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor40b[i] <- cor(rowSums(sim.b.40[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor30b[i] <- cor(rowSums(sim.b.30[[i]]),rowSums(sim))}
for (i in 1:nsim) {cor20b[i] <- cor(rowSums(sim.b.20[[i]]),rowSums(sim))}


bins <- c(20,30,40,50,60,70,80,90)

vvrs2 <- rbind(brs20,brs30,brs40,brs50,brs60,brs70,brs80,brs90)
vv2rs2 <- rbind(sqrt(b2rs20),sqrt(b2rs30),sqrt(b2rs40),sqrt(b2rs50),sqrt(b2rs60),sqrt(b2rs70),sqrt(b2rs80),sqrt(b2rs90))
btrs2 <- rbind(btrs20,btrs30,btrs40,btrs50,btrs60,btrs70,btrs80,btrs90)
ors2 <- rbind(ors20,ors30,ors40,ors50,ors60,ors70,ors80,ors90)
vv2 <- rbind(b220,b230,b240,b250,b260,b270,b280,b290)



btsim <- rbind(btsim20,btsim30,btsim40,btsim50,btsim60,btsim70,btsim80,btsim90)
bsim <- rbind(bsim20,bsim30,bsim40,bsim50,bsim60,bsim70,bsim80,bsim90)

plot(bins,bsim[,2],ylim=c(0.07,0.23),type='l',xlab='Sampling Fraction',ylab='Root Mean Squared-Error',col='blue',lwd=2)
points(bins,btsim[,2],type='l',lwd=2,col='red')
polygon(c(bins,rev(bins)),c(btsim[,1],rev(btsim[,3])),col=rgb(1,0,0,0.25),border=NA)
polygon(c(bins,rev(bins)),c(bsim[,1],rev(bsim[,3])),col=rgb(0,0,1,0.25),border=NA)
legend('topright',c('Bootstrap Resampling','Empirical Bayesian Estimation','95% Confidence Intervals'),cex=0.75,bty='n',lwd=c(2,2,10),col=c('red','blue',rgb(1,0,0,0.25)))


plot(bins,vv2[,2],ylim=c(0.001,0.008),type='l',xlab='Sampling Fraction',ylab='Mean Squared-Error',col='blue',lwd=2)
points(bins,btrs2[,2],type='l',lwd=2,col='red')
polygon(c(bins,rev(bins)),c(btrs2[,1],rev(btrs2[,3])),col=rgb(1,0,0,0.25),border=NA)
polygon(c(bins,rev(bins)),c(vv2[,1],rev(vv2[,3])),col=rgb(0,0,1,0.25),border=NA)
legend('topright',c('Bootstrap Resampling','Empirical Bayesian Estimation','95% Confidence Intervals'),cex=0.75,bty='n',lwd=c(2,2,10),col=c('red','blue',rgb(1,0,0,0.25)))

plot(logcer,log(rowMeans(boot.st)),type='l',col='red',xlab='Log Sample Size',ylab='Log Mean Squared Error')
points(logcer,log10(colMeans(zz2)),type='l',col='blue')

cor2.95q <- c(quantile(cor90.2,probs=c(0.95)),quantile(cor80.2,probs=c(0.95)),quantile(cor70.2,probs=c(0.95)),quantile(cor60.2,probs=c(0.95)),quantile(cor50.2,probs=c(0.95)),quantile(cor40.2,probs=c(0.95)),quantile(cor30.2,probs=c(0.95)),quantile(cor20.2,probs=c(0.95)))
cor2.5q <- c(quantile(cor90.2,probs=c(0.05)),quantile(cor80.2,probs=c(0.05)),quantile(cor70.2,probs=c(0.05)),quantile(cor60.2,probs=c(0.05)),quantile(cor50.2,probs=c(0.05)),quantile(cor40.2,probs=c(0.05)),quantile(cor30.2,probs=c(0.05)),quantile(cor20.2,probs=c(0.05)))
