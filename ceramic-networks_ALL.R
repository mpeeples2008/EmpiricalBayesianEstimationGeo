library(tidyverse)
require(statnet)
library(igraph)
library(vegan)
library(reshape)


load('apportioned.RData')

############ SETUP INTERVALS

wareapp3 <- wareapp

per.list2 <- c('AD900_950','AD950_1000','AD1000_1050','AD1050_1100',
               'AD1100_1150','AD1150_1200','AD1200_1250','AD1250_1300')
lookup.list <- list(c(25:26),c(27:28),c(29:30),c(31:32),c(33:34),c(35:36),
                    c(37:38),c(39:40))

for (i in 1:length(per.list2)) {
  warevar <- wareapp3[,1:2]
  cer <- wareapp3[,lookup.list[[i]]+2]
  if (!is.null(ncol(cer))) {cer <- rowSums(cer)}
  cer <- cbind(warevar,cer)
  out <- cast(cer,Site~SWSN_Ware,fun.aggregate = sum)
  out[is.na(out)] <- 0
  rownames(out) <- out[,1]
  
  assign(paste(per.list2[i]),out)}


####################################################################################################
###TRIM AND RENAME CERAMIC AND ATTRIBUTE DATA FRAMES BY PERIOD######################################
####################################################################################################

cer.output <- list(AD900_950,AD950_1000,AD1000_1050,AD1050_1100,
                     AD1100_1150,AD1150_1200,AD1200_1250,AD1250_1300)

attr <- read.csv('SWSN_Sites.csv',header=T)

attr.new <- attr[,c(1,2,10,11)]
attr.comb2 <- attr.new[which(attr.new[,3]>5000),]

write.table(attr.comb2,file='attr_all.csv',sep=',',row.names=F)
attr.comb <- read.table('attr_all.csv',header=T,sep=',')

for (i in 1:length(per.list2)) {
  cer <- (cer.output[[i]])
  cer <- cer[,-1]
  trim <- which(rowSums(as.matrix(cer))>=10)
  cer <- cer[trim,]
  cer2 <- as.matrix(cer[,which(colSums(cer)>0)])
  colnames(cer2) <- colnames(cer)[which(colSums(cer)>0)]
  row.names(cer2) <- row.names(cer)
    attr.temp <- attr.comb[which(as.vector(attr.comb[,2]) %in% row.names(cer2)),]
    cer2 <- as.matrix(cer2[which(row.names(cer2) %in% attr.temp[,2]),])
    attr.temp2 <- as.data.frame(matrix(NA,nrow(cer2),4))
    for (j in 1:nrow(attr.temp2)) {
      attr.temp2[j,] <- attr.temp[which(row.names(cer2)[j]==attr.temp[,2])[1],]
    }
    row.names(attr.temp2) <- row.names(cer2)
    attr.temp2[,1] <- row.names(attr.temp2)
  assign(paste((per.list2[i]),"cer",sep=""),cer2)
  assign(paste((per.list2[i]),"attr",sep=""),attr.temp2)}


coord.list <- list(AD900_950attr,AD950_1000attr,AD1000_1050attr,
                   AD1050_1100attr,AD1100_1150attr,AD1150_1200attr,AD1200_1250attr,AD1250_1300attr)
cer.list <- list(AD900_950cer,AD950_1000cer,AD1000_1050cer,
                 AD1050_1100cer,AD1100_1150cer,AD1150_1200cer,AD1200_1250cer,AD1250_1300cer)

####################################################################################################
###EMPIRICAL BAYESIAN ESTIMATION####################################################################
####################################################################################################

source('EmpiricalBayesianEstimation.R')


AD900_950cer2 <- ebe(cer.list[[1]],coord.list[[1]])
AD950_1000cer2 <- ebe(cer.list[[2]],coord.list[[2]])
AD1000_1050cer2 <- ebe(cer.list[[3]],coord.list[[3]])
AD1050_1100cer2 <- ebe(cer.list[[4]],coord.list[[4]])
AD1100_1150cer2 <- ebe(cer.list[[5]],coord.list[[5]])
AD1150_1200cer2 <- ebe(cer.list[[6]],coord.list[[6]])
AD1200_1250cer2 <- ebe(cer.list[[7]],coord.list[[7]])
AD1250_1300cer2 <- ebe(cer.list[[8]],coord.list[[8]])


####################################################################################################
### CREATE SIMILARITY MATRICES AND NETWORKS BY PERIOD  #############################################
####################################################################################################

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



AD900sim <- sim.mat(AD900_950cer2)
AD950sim <- sim.mat(AD950_1000cer2)
AD1000sim <- sim.mat(AD1000_1050cer2)
AD1050sim <- sim.mat(AD1050_1100cer2)
AD1100sim <- sim.mat(AD1100_1150cer2)
AD1150sim <- sim.mat(AD1150_1200cer2)
AD1200sim <- sim.mat(AD1200_1250cer2)
AD1250sim <- sim.mat(AD1250_1300cer2)

####################################################################################################
###NETWORK DEFINITION###############################################################################
####################################################################################################

set.thresh <- 0.749 ### CREATE SIMILARITY THRESHOLD FOR DEFINING TIE AS PRESENT

AD900net <- network(event2dichot(AD900sim,method='absolute',thresh=set.thresh),directed=F)
AD900net %v% 'vertex.names' <- row.names(AD900sim)
AD950net <- network(event2dichot(AD950sim,method='absolute',thresh=set.thresh),directed=F)
AD950net %v% 'vertex.names' <- row.names(AD950sim)
AD1000net <- network(event2dichot(AD1000sim,method='absolute',thresh=set.thresh),directed=F)
AD1000net %v% 'vertex.names' <- row.names(AD1000sim)
AD1050net <- network(event2dichot(AD1050sim,method='absolute',thresh=set.thresh),directed=F)
AD1050net %v% 'vertex.names' <- row.names(AD1050sim)
AD1100net <- network(event2dichot(AD1100sim,method='absolute',thresh=set.thresh),directed=F)
AD1100net %v% 'vertex.names' <- row.names(AD1100sim)
AD1150net <- network(event2dichot(AD1150sim,method='absolute',thresh=set.thresh),directed=F)
AD1150net %v% 'vertex.names' <- row.names(AD1150sim)
AD1200net <- network(event2dichot(AD1200sim,method='absolute',thresh=set.thresh),directed=F)
AD1200net %v% 'vertex.names' <- row.names(AD1200sim)
AD1250net <- network(event2dichot(AD1250sim,method='absolute',thresh=set.thresh),directed=F)
AD1250net %v% 'vertex.names' <- row.names(AD1250sim)


net.list <- list(AD900net,AD950net,AD1000net,AD1050net,AD1100net,AD1150net,
                 AD1200net,AD1250net)


save.image('apportioned.RData')



