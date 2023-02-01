library(VGAM)
library(broom)
library(tidyr)
library(dplyr)
library(DirichletMultinomial)
library(classInt)
library(car)
library(lmtest)

sim.mat <- function(x) {
  names <- row.names(x)
  x <- na.omit(x)
  x <- prop.table(as.matrix(x),1)*100
  rd <- dim(x)[1]
  results <- matrix(0,rd,rd)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- as.numeric(x[s1, ])
      x2Temp <- as.numeric(x[s2, ])
      results[s1,s2] <- 200 - (sum(abs(x1Temp - x2Temp)))}}
  row.names(results) <- names
  colnames(results) <- names
  results <- results/200
  results <- round(results,3)
  return(results)}


ebe <- function(cer,attr) {

cer.prop <- prop.table(as.matrix(cer),1)
cer.sum <- rowSums(cer)
sim2 <- sim.mat(cer)


final <- matrix(0,nrow(cer),ncol(cer))
colnames(final) <- colnames(cer)
row.names(final) <- row.names(cer)
ddist <- as.matrix(dist(attr[,3:4]))
ddist <- max(ddist)-ddist
row.names(ddist) <- row.names(cer)

ddist <- (ddist-min(ddist))/(max(ddist)-min(ddist))
mult1 <- matrix(NA,nrow(ddist),ncol(ddist))
for (i in 1:length(sim2)) {
  mult1[i] <- weighted.mean(x=c(sim2[i],ddist[i]),w=c(0.5,0.5))}
diag(mult1) <- 0

for (m in 1:nrow(cer)) {
  cer.a <- as.matrix(cer)
  row.names(cer.a) <- row.names(cer)
  colnames(cer.a) <- colnames(cer)
  
  cer.a <- cer*mult1[m,]
  cer.a <- cer.a[-m,]
  cer.a <- ceiling(cer.a)
  #####
  cer.a <- cer[which(mult1[m,]>quantile(mult1[m,],probs=c(0.5))),]  #### 0.5 default
  cer.a <- ceiling(cer.a)
  cer.ret2 <- which(rowSums(cer.a)>5)
  cer.a <- cer.a[cer.ret2,]
  cer.b <- as.matrix(cer.a)
  row.names(cer.b) <- row.names(cer.a)
  colnames(cer.b) <- colnames(cer.a)
  cer.a <- cer.b
  
  dm_fit <- DirichletMultinomial::dmn(cer.a,1)
  
  tidy.DMN <- function(x, ...) {
    ret <- as.data.frame(x@fit)
    tibble::as_tibble(fix_data_frame(ret, c("conf.low", "estimate", "conf.high")))}
  
  dm_params <- tidy(dm_fit)

  par.e <- dm_params$estimate
  par_total <- sum(par.e)
  
  item_mean <- dm_params$estimate / sum(dm_params$estimate)

  cer.p <- prop.table(cer.a,1)
  out.type <- cer.a
  for (i in 1:ncol(cer.a)) {
    out.type[,i] <- (cer.a[,i]+par.e[i])/(rowSums(cer.a)+par_total)}
  
  site <- t(as.matrix(as.matrix(cer)[m,]))
  out.site <- site
  out.site <- prop.table(as.matrix(out.site),1)
  for (i in 1:ncol(site)) {
    out.site[1,i] <- (site[1,i]+par.e[i])/(rowSums(site)+par_total)}
  row.names(out.site) <- v1 <- row.names(cer)[m]
  colnames(out.site) <- v2 <- colnames(cer)
  
  for(i in 1:ncol(final)) {
    final[m,i] <- out.site[1,i]
    row.names(final)[m] <- row.names(out.site)}
}


ab <- abs(final-cer.prop)

logcer <- log(rowSums(cer))
logab <- log(rowSums(ab))
fit <- lm(logcer~logab)

# label outliers
outlier <- names(outlierTest(fit)$rstudent)


final <- final*rowSums(cer)
final[outlier,] <- cer[outlier,]

output <- list(data = final, cer_count = rowSums(cer), difference = rowSums(ab), fit = fit, outlier = outlier)

return(output)
}

