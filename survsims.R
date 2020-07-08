## 4/17/2019
## sims to look at method with aft models
library(survival)
library(ggplot2)
library(dplyr)


NSIMS = 1000
n = seq(30, 130,by=10)
len = length(n)

alpha = 0.5
beta = 0.25
sigmax = 1
s.shape = 4 

#function to calc min n needed for 70% power
powercomp_surv <- function(n, gam, NSIMS){
  len <- length(n)
  #data <- dat
  
  pow.pred<-rep(NA,len)
  pow.nopred<-rep(NA,len)
  
  ti.pred <- rep(NA,len)
  ti.nopred <- rep(NA,len)
  
  status <- rep(NA, len)
  
  set.seed(12345)
  pb <- txtProgressBar(0,len*NSIMS,style=3)
  for(j in 1:len){
    pvals.nopred.p<-rep(NA,NSIMS)
    pvals.pred.p <- rep(NA,NSIMS)
    
    pvals.nopred.ti<-rep(NA,NSIMS)
    pvals.pred.ti <- rep(NA,NSIMS)
    
    status.ti.perc <- rep(NA, NSIMS)
    
    for(i in 1:NSIMS){
      #report total trial size
      n.tx <- n[j]/2
      
      #treatment arm 1, ctl
      A.tx.1<-rep(0, n.tx) #trial ctl
      X.tx.1 <- rnorm(n.tx, 0, sigmax)
      surv_time.ti <- rweibull(n.tx, shape = s.shape, scale = exp(alpha + beta*X.tx.1))
      surv_time.pow <- rweibull(n.tx, shape = s.shape, scale = exp(alpha + beta*X.tx.1 + gam*A.tx.1))
      cens_time = 36
      
      status.ti.1 <-  ifelse (surv_time.ti <= cens_time, 1, 0) # event indicator
      status.pow.1 <- ifelse (surv_time.pow <= cens_time, 1, 0) # event indicator
      
      Y.tx.1.ti <- ifelse(surv_time.ti <= cens_time, surv_time.ti, cens_time)
      Y.tx.1.p <- ifelse(surv_time.pow <= cens_time, surv_time.pow, cens_time)
      
      #treatment arm 2, tx
      A.tx.2<-rep(1,n.tx) #trialtx
      X.tx.2 <- rnorm(n.tx, 0, sigmax)
      
      surv_time.ti <- rweibull(n.tx, shape = s.shape, scale = exp(alpha + beta*X.tx.2))
      surv_time.pow <- rweibull(n.tx, shape = s.shape, scale = exp(alpha + beta*X.tx.2 + gam*A.tx.2))
      cens_time = 36
      
      status.ti.2 <-  ifelse (surv_time.ti <= cens_time, 1, 0) # event indicator
      status.pow.2 <- ifelse (surv_time.pow <= cens_time, 1, 0) # event indicator
      
      Y.tx.2.ti <-ifelse(surv_time.ti <= cens_time, surv_time.ti, cens_time)
      Y.tx.2.p <- ifelse(surv_time.pow <= cens_time, surv_time.pow, cens_time)
      
      #collecting pvals
      Y.tx.ti <- c(Y.tx.1.ti, Y.tx.2.ti)
      status.ti <- c(status.ti.1, status.ti.2)
      Y.tx.p <- c(Y.tx.1.p, Y.tx.2.p)
      status.p <- c(status.pow.1, status.pow.2)
      X.tx <- c(X.tx.1, X.tx.2)
      A.tx <- c(A.tx.1, A.tx.2)
      
      pvals.nopred.p[i] <- summary(survreg(Surv(time = Y.tx.p, status.p) ~ A.tx, dist = "loglogistic"))$table[2,4]
      pvals.pred.p[i] <- summary(survreg(Surv(time = Y.tx.p, event = status.p) ~ A.tx + X.tx, dist = "loglogistic"))$table[2,4]
      
      pvals.nopred.ti[i] <- summary(survreg(Surv(time = Y.tx.ti, event = status.ti) ~ A.tx, dist = "loglogistic"))$table[2,4]
      pvals.pred.ti[i] <- summary(survreg(Surv(time = Y.tx.ti, event = status.ti) ~ A.tx + X.tx, dist = "loglogistic"))$table[2,4]
      
      status.ti.perc[i] <- sum(status.ti)/length(status.ti)
      
      setTxtProgressBar(pb, (j-1)*NSIMS+i)
    }
    status[j] <- 1-mean(status.ti.perc)
    pow.pred[j] <- sum(pvals.pred.p < 0.05)/NSIMS 
    pow.nopred[j] <- sum(pvals.nopred.p < 0.05)/NSIMS 
    
    ti.pred[j] <- sum(pvals.pred.ti < 0.05)/NSIMS 
    ti.nopred[j] <- sum(pvals.nopred.ti < 0.05)/NSIMS 
  }
  close(pb)
  
  min.n <- NULL
  min.n <- c(min.n, n[which(pow.nopred ==  pow.nopred[pow.nopred >= 0.8][1])]) #no HCs
  min.n <- c(min.n, n[which(pow.pred ==  pow.pred[pow.pred >= 0.8][1])]) #with HCs
  
  return(min.n) 
}
