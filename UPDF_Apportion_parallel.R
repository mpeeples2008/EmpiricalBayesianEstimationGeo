require(sqldf)
require(compiler)
require(reshape2)
require(msm)
require(snowfall)
require(parallel)

# Multicore setup
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)

###############################################################
###############################################################
###############################################################

updf <- function(site,cer.type,ct,start,end,chron,interval=5,cutoff=0.1,min.period=25) {

## Remove any type with a count of 0 from all vectors
  trim <- which(ct>0)
  cer.type <- as.vector(cer.type[trim])
  site <- as.vector(site[1])
  start <- start[trim]
  end <- end[trim]
  ct <- ct[trim]
  chron <- chron[trim]

 
## define function for rounding by base defined in "interval" argument and round start and end dates for each type
  mround <- function(x,base){base*round(x/base)}
  start <- mround(start,interval)
  end <- mround(end,interval)

## find minimal ceramic intervals based on rounded type date overlaps
  years <- unique(sort(c(start,end)))
    period <- do.call(rbind,foreach(i=1:(length(years)-1)) %dopar% rbind(c(years[i],years[i+1]))) 
    # Ensure that no period is shorter than the specificed min.period argument and combine periods to satisfy this requirement
    if (length(period)>2) {
    per.short <- foreach(i=1:(nrow(period)-1)) %dopar% (period[i+1,1]-period[i,1])
    comb.per <- which(per.short<min.period)
      if (length(comb.per)>0){ period <- period[-comb.per,]
       for (i in 1:(nrow(period)-1)) {period[i,2] <- period[i+1,1]}}
    }  
    
## Define lists of ceramic period lengths and site period lengths
  cer.per <- foreach(i=1:length(start)) %dopar% (seq((start[i]:end[i]))+(start[i]-1))
  per.len <- foreach(i=1:nrow(period)) %dopar% (seq(period[i,1]:period[i,2])+(period[i,1])-1)
 
## Create Uniform Distance dataframe via parallel function
  ucalc <- function(a,b) {
    out <- (length(intersect(a,b))-1)/(length(a)-1)
      if (out<0) {out <- 0}
        return(out)}
  
  udist <- foreach(a=cer.per,.combine='rbind') %:% foreach(b=per.len,.combine='c') %dopar% (ucalc(a,b))
  if (length(udist)==1) {udist <- as.matrix(udist)}
  colnames(udist) <- paste('AD',period[,1],'-',period[,2],sep='')
  
## Create dataframe of prior values and sum to calculate prior probabilities
  prior.tab <- udist*ct
  if (ncol(prior.tab)>1) {
  if (length(prior.tab[which(chron==1),])==ncol(prior.tab)) {prior <- prior.tab[which(chron==1),]/sum(prior.tab[which(chron==1),])} else
  {prior <- colSums(prior.tab[which(chron==1),])/sum(prior.tab[which(chron==1),])}} else {
    prior <- udist[1,]
  }
  prior[is.na(prior)] <- 0
  
## Create dataframe of count probabilities
  pij <- sweep(prior.tab,2,colSums(prior.tab),'/')
  pij[is.na(pij)] <- 0
  
## Create dataframe of probabilities based on uniform distribution
  uij <- sweep(udist,2,colSums(udist),'/')
  uij[is.na(uij)] <- 0
  
## Create dataframe of standard deviations of uniform distribution probabilities
  sd.t <- sqrt(sweep((uij*(1-uij)),2,colSums(ceiling(udist)),'/'))
  
## Create dataframe of conditionals and calculate conditional probabilities
  cij <- dnorm(pij,uij,sd.t)
  cij[which(uij-sd.t==1)] <- dnorm(1,1,1) # For intervals with only one ceramic type, set to standard normal
  is.na(cij) <- !is.finite(cij) 

  # calculate conditional probability and rescale from 0-1
  conditional <- apply(cij,2,mean,na.rm=T)
  conditional[is.na(conditional)] <- 0
  if (sum(conditional)>0) {conditional <- conditional/sum(conditional)}

## Calculate posterior proportions and remove any generated NA
  posterior <- foreach(i=1:length(prior),.combine='c') %dopar% ((prior[i]*conditional[i]))
  posterior <- posterior/sum(posterior)
  posterior[is.na(posterior)] <- 0
  
## Deal with edge cases where sum of posterior probabilities = 0 or conditional = prior
  if(sum(posterior)==0) {posterior <- prior}
  if (identical(conditional,prior)) {posterior <- prior}

## Create dataframe with ceramics apportioned by uniform distribution (udist.out) and posterior estimates (post.out)
  cer.lab <- as.data.frame(cbind(rep(site,length(ct)),cer.type,ct,start,end,chron))
  colnames(cer.lab) <- c('Site','Type','Count','StartDate','EndDate','Chronology')
  
  udist.out <- cbind(cer.lab,(udist*ct))
  
  post.tab <- (((sweep(prior.tab,2,posterior,'*'))/rowSums(sweep(prior.tab,2,posterior,'*')))*ct)
  post.tab[is.na(post.tab)] <- 0
  colnames(post.tab) <- paste('AD',period[,1],'-',period[,2],sep='')
  post.out <- cbind(cer.lab,post.tab)
  
# calculated beginning (lwr) and ending (upr) dates based on the first and last period with posterior probability about the selected threshold
  lwr <- period[min(which(posterior>cutoff*max(posterior))),1]
  upr <- period[max(which(posterior>cutoff*max(posterior))),2]
  
# Reduced ceramic list for occupation period as selected by cutoff
  extra <- as.matrix(post.tab[,-c((min(which(posterior>cutoff*max(posterior)))):max(which(posterior>cutoff*max(posterior))))])
  colnames(extra) <- colnames(post.tab)[-c((min(which(posterior>cutoff*max(posterior)))):max(which(posterior>cutoff*max(posterior))))]
  occ.cer <- as.matrix(post.tab[,c((min(which(posterior>cutoff*max(posterior)))):max(which(posterior>cutoff*max(posterior))))])
  colnames(occ.cer) <- colnames(post.tab)[c((min(which(posterior>cutoff*max(posterior)))):max(which(posterior>cutoff*max(posterior))))]

  if (length(extra)>0) {
  occ.cer.out <- occ.cer/rowSums(occ.cer)*(rowSums(extra))+occ.cer} else {
    occ.cer.out <- occ.cer
  }
  
  occ.cer.out[is.na(occ.cer.out)] <- 0
  
  occ.cer.out <- cbind(cer.lab,occ.cer.out)
  
  per.occ <- period[c((min(which(posterior>cutoff*max(posterior)))):max(which(posterior>cutoff*max(posterior)))),]
  
## Create output list and return
  out.list <- list(site=site,prior=prior,posterior=posterior,conditional=conditional,period=period,samp.size=sum(ct),udist=udist.out,postdist=post.out,occ=cbind(lwr,upr),per.occ=per.occ,cer.occ=occ.cer.out)
  return(out.list)
}


##############################################
##############################################
##############################################
# Read in ceramic and attribute data in long format and check for missing values

dat.long <- read.csv('data_all.csv',header=T)
sites <- unique(dat.long$Site)

#Create plot of Bayesian procedure by site/Intra-site context
out.cer <- list()
beg.per <- 300
end.per <- 2000
interval <- 25
full.per <- seq(beg.per,end.per,by=interval)
full.lab <- do.call(rbind,foreach(i=1:(length(full.per)-1)) %dopar% rbind(c(full.per[i],full.per[i+1])))


for (i in 1:length(sites)) {
  qv <- dat.long[which(dat.long$Site==sites[i]),]
  vv <- qv[which(qv$cyberSW==1),]
  if(nrow(vv)>1) {
      up <- updf(site=qv$Site,cer.type=qv$Type,ct=qv$Count,start=qv$Begin,end=qv$End,chron=qv$cyberSW,interval=25,cutoff=0.25,min.period=10)
  }
  per.tab <- as.matrix(up$per.occ)
  if (length(per.tab)==2) {per.tab <- t(per.tab)}
  int.brk <- (per.tab[,2]-per.tab[,1])/interval
  tmp <- as.matrix(up$cer.occ[,7:ncol(up$cer.occ)])
  app <- matrix(0,nrow(tmp),sum(int.brk))
  app.per <- seq(min(up$per.occ),max(up$per.occ),by=interval)
  
  j <- 0
  for (m in 1:length(int.brk)) {
    j <- int.brk[m]+j
    app[,(j-int.brk[m]+1:int.brk[m])] <- tmp[,m]/int.brk[m]}
  full.app <- matrix(0,nrow(app),length(full.per)-1)
  colnames(full.app) <- paste('AD',full.lab[,1],sep='')
  
  if (min(up$per.occ)<beg.per) {app <- app[,min(which(app.per %in% seq(beg.per,end.per,by=interval))):ncol(app)]
  app.per <- app.per[-(1:min(which(app.per %in% seq(beg.per,end.per,by=interval)))-1)]}
  
  if (max(up$per.occ)>end.per) {app <- app[,1:(max(which(app.per %in% seq(beg.per,end.per,by=interval)))-1)]
  app.per <- app.per[-((max(which(app.per %in% seq(beg.per,end.per,by=interval)))+1):length(app.per))]}
  
  col.start <- which(full.lab[,1]==min(app.per))
  
  full.app[,col.start:(col.start+dim(app)[2]-1)] <- app
  out.cer[[i]] <- cbind(up$cer.occ[,1:6],full.app)
  
  out.cer
}


final <- do.call(rbind,out.cer) 


####################################################################################################
###CONVERT TYPE DATA TO WARES ######################################################################
####################################################################################################

Ceramic_type_master <- read.table(file='Ceramic_type_master.csv',sep=',',header=T)
a.lab <- paste('AD',full.lab[,1],sep='')

sqltext <- 'SELECT final.Site, Ceramic_type_master.SWSN_Ware, '
for (i in 1:length(a.lab)) {
  sqltext <- paste(sqltext,"sum(final.",a.lab[i],")*1 AS ",a.lab[i],", ",sep="")}
sqltext <- substr(sqltext,1,nchar(sqltext)-2)
sqltext <- paste(sqltext,"FROM Ceramic_type_master INNER JOIN final ON Ceramic_type_master.SWSN_Type = final.Type
                 WHERE (Ceramic_type_master.Decoration='bichrome' Or Ceramic_type_master.Decoration='polychrome' Or Ceramic_type_master.Decoration='undifferentiated dec') AND (Ceramic_type_master.SWSN_Ware Not Like 'Undiff%')
                 GROUP BY final.Site, Ceramic_type_master.SWSN_Ware")

wareapp <- sqldf(sqltext)

rm(app, cl, full.app, full.lab, out.cer, per.tab, qv, tmp, up, vv, a.lab, app.per, beg.per, col.start, end.per, full.per,
   i, int.brk, interval, j, m, sqltext)

save.image('apportioned.RData')