source("EmpiricalBayesianEstimation.R")

out <- ebe(cer.list[[1]], coord.list[[1]])
sim2 <- sim.mat(out$data)


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

test <- simulation2(cer.list[[1]], 0.5)

sim.test <- sim.mat(test[[1]])

nums <- 100
colscale <- colorRampPalette(c('purple4','blue','skyblue'))(nums)
col2 <- as.numeric(cut(cer.sum,c(classIntervals(cer.sum,nums,style='quantile')$brks,Inf),right=F))

plot(logcer,logab,pch=16,col=colscale[col2],xlab='Log of Sample Size',ylab='Log of Posterior Difference')
fit <- lm(logcer~logab)
lm1 <- predict(fit)
lines(sort(lm1),logab[order(lm1)],col='red',lty=2)

# label outliers
outlier <- names(outlierTest(fit)$rstudent)
temp <- NULL
for (i in 1:length(outlier)) {
  temp[i] <- which(row.names(cer)==outlier[i])}
outlier <- temp
text(logcer[outlier],logab[outlier],labels=row.names(cer)[outlier],cex=0.5,pos=1)