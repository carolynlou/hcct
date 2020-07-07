twoarmcomp <- function(NSIMS, n = seq(50, 150, by = 50), gam, pow = 0.8, data){
  len <- length(n)
  
  pow.pred<-rep(NA,len)
  pow.nopred<-rep(NA,len)
  pow.lin <- rep(NA,len)
  
  ti.pred <- rep(NA,len)
  ti.nopred <- rep(NA,len)
  ti.lin <- rep(NA,len)
  
  set.seed(4529)
  pb <- txtProgressBar(0,len*NSIMS,style=3)
  for(j in 1:len){
    pvals.pred.p<-rep(NA,NSIMS)
    pvals.nopred.p<-rep(NA,NSIMS)
    pvals.lin.p <- rep(NA,NSIMS)
    
    pvals.pred.ti <- rep(NA,NSIMS)
    pvals.nopred.ti<-rep(NA,NSIMS)
    pvals.lin.ti <- rep(NA,NSIMS)
    
    for(i in 1:NSIMS){
      #split ADNI data into 2
      tx1 <- data[sample(nrow(data),nrow(data)/2,replace=F),]
      tx2 <- suppressMessages(anti_join(data,tx1))
      
      n.tx <- n[j]/2
      
      #treatment arm 1
      ftx.1 <- tx1[sample(nrow(tx1),n.tx,replace=T),] #"bootstrap" a future tx sample of MCI pts, tx arm 1
      A.tx.1<-rep(1,n.tx) #trial tx'd
      X.tx.1 <- ftx.1$SPARE.AD
      Y.tx.1.p<-ftx.1$mem.ch+A.tx.1*gam #add tx effect to assess power
      Y.tx.1.ti<-ftx.1$mem.ch
      
      #treatment arm 2
      ftx.2 <- tx2[sample(nrow(tx2),n.tx,replace=T),] #a future tx sample of MCI pts, tx arm 2
      A.tx.2<-rep(0,n.tx) #trial controls
      X.tx.2 <- ftx.2$SPARE.AD
      Y.tx.2.p<-ftx.2$mem.ch+A.tx.2*gam
      Y.tx.2.ti<-ftx.2$mem.ch
      
      #collecting power pvals
      Y.tx.ti <- c(Y.tx.1.ti, Y.tx.2.ti)
      Y.tx.p <- c(Y.tx.1.p, Y.tx.2.p)
      X.tx <- c(X.tx.1, X.tx.2)
      A.tx <- c(A.tx.1, A.tx.2)
      
      pvals.nopred.p[i] <- summary(lm(Y.tx.p~A.tx))$coef[2,4]
      pvals.nopred.ti[i] <- summary(lm(Y.tx.ti~A.tx))$coef[2,4]
      
      pvals.pred.p[i] <- summary(lm(Y.tx.p~A.tx+X.tx))$coef[2,4]
      pvals.pred.ti[i] <- summary(lm(Y.tx.ti~A.tx+X.tx))$coef[2,4]
      
      setTxtProgressBar(pb, (j-1)*NSIMS+i)
    }
    pow.pred[j] <- sum(pvals.pred.p<0.05)/NSIMS 
    pow.nopred[j] <- sum(pvals.nopred.p<0.05)/NSIMS 
    
    ti.pred[j] <- sum(pvals.pred.ti<0.05)/NSIMS 
    ti.nopred[j] <- sum(pvals.nopred.ti<0.05)/NSIMS 
    
  }
  close(pb)
  
  #min N needed to get e.g. 70% power with effect size 0.24
  min.n <- NULL
  min.n <- c(min.n, n[which(pow.nopred ==  pow.nopred[pow.nopred >= pow][1])]) #no HCs. power = 0.6094
  min.n <- c(min.n, n[which(pow.pred ==  pow.pred[pow.pred >= pow][1])]) #with HCs. 0.633
  
  return(min.n)
}
