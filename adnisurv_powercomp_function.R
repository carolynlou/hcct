powercomp_surv_adni <- function(NSIMS, n, gam, pow = 0.8, data){
  len <- length(n)
  data <- dat
  
  pow.pred<-rep(NA,len)
  pow.nopred<-rep(NA,len)
  
  ti.pred <- rep(NA,len)
  ti.nopred <- rep(NA,len)
  
  cens <- rep(NA, len)
  
  set.seed(12345)
  pb <- txtProgressBar(0,len*NSIMS,style=3)
  for(j in 1:len){
    pvals.nopred.p<-rep(NA,NSIMS)
    pvals.pred.p <- rep(NA,NSIMS)
    
    pvals.nopred.ti<-rep(NA,NSIMS)
    pvals.pred.ti <- rep(NA,NSIMS)
    
    cens.ti.perc <- rep(NA, NSIMS)
    
    for(i in 1:NSIMS){
      #split ADNI data into 2
      tx1 <- data[sample(nrow(data),nrow(data)/2,replace=F),]
      tx2 <- suppressMessages(anti_join(data, tx1))
      
      n.tx <- n[j]/2
      
      #treatment arm 1. tx
      ftx.1 <- tx1[sample(nrow(tx1),n.tx,replace=F),] #"bootstrap" a future tx sample of MCI pts, tx arm 1
      A.tx.1 <- rep(1,n.tx) #trial tx'd
      X.tx.1 <- ftx.1$SPARE.AD
      
      Y.tx.1.p <- pmin(gam*ftx.1$TimeToConv, 36)
      cens.1.p <- 1*(Y.tx.1.p < 36) 
      
      Y.tx.1.ti <- ftx.1$TimeToConv
      cens.1.ti <- 1*(Y.tx.1.ti < 36) 
      
      #treatment arm 2, ctl
      ftx.2 <- tx2[sample(nrow(tx2),n.tx,replace=F),] #a future tx sample of MCI pts, tx arm 2
      A.tx.2 <- rep(0, n.tx) #trial controls
      X.tx.2 <- ftx.2$SPARE.AD
      
      Y.tx.2.p <- 1*ftx.2$TimeToConv
      cens.2.p <- 1*(Y.tx.2.p < 36) 
      
      Y.tx.2.ti <- ftx.2$TimeToConv
      cens.2.ti <- 1*(Y.tx.2.ti < 36) 
      
      #collecting pvals
      Y.tx.ti <- c(Y.tx.1.ti, Y.tx.2.ti)
      cens.ti <- c(cens.1.ti, cens.2.ti)
      Y.tx.p <- c(Y.tx.1.p, Y.tx.2.p)
      cens.p <- c(cens.1.p, cens.2.p)
      X.tx <- c(X.tx.1, X.tx.2)
      A.tx <- c(A.tx.1, A.tx.2)
      
      pvals.nopred.p[i] <- summary(survreg(Surv(time = Y.tx.p, cens.p) ~ A.tx, dist = "loglogistic"))$table[2,4]
      pvals.pred.p[i] <- summary(survreg(Surv(time = Y.tx.p, event = cens.p) ~ A.tx + X.tx, dist = "loglogistic"))$table[2,4]
      
      pvals.nopred.ti[i] <- summary(survreg(Surv(time = Y.tx.ti, event = cens.ti) ~ A.tx, dist = "loglogistic"))$table[2,4]
      pvals.pred.ti[i] <- summary(survreg(Surv(time = Y.tx.ti, event = cens.ti) ~ A.tx + X.tx, dist = "loglogistic"))$table[2,4]
      
      cens.ti.perc[i] <- sum(cens.ti)/length(cens.ti)
      
      setTxtProgressBar(pb, (j-1)*NSIMS+i)
    }
    cens[j] <- 1-mean(cens.ti.perc)
    pow.pred[j] <- sum(pvals.pred.p < 0.05)/NSIMS 
    pow.nopred[j] <- sum(pvals.nopred.p < 0.05)/NSIMS 
    
    ti.pred[j] <- sum(pvals.pred.ti < 0.05)/NSIMS 
    ti.nopred[j] <- sum(pvals.nopred.ti < 0.05)/NSIMS 
  }
  close(pb)
  
  min.n <- NULL
  min.n <- c(min.n, n[which(pow.nopred ==  pow.nopred[pow.nopred >= pow][1])]) #no HCs
  min.n <- c(min.n, n[which(pow.pred ==  pow.pred[pow.pred >= pow][1])]) #with HCs
  
  return(min.n) 
}
