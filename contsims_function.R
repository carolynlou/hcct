# 1/17/2020

library(ggplot2)

#### continuous ####
twoarmcont <- function(NSIMS, n = seq(10,200,by=10), gam, pow = 0.8, rand = F){
  alpha <- 2
  beta<-4
  sigmax<-0.25
  sigmae<-1
  
  len <- length(n)
  pow.pred<-rep(NA,len)
  pow.nopred<-rep(NA,len)
  ti.pred <- rep(NA, len)
  ti.nopred <- rep(NA,len)
  
  set.seed(12345)
  pb <- txtProgressBar(0,len*NSIMS,style=3)
  for(j in 1:len){
    pvals.pred.p<-rep(NA,NSIMS)
    pvals.nopred.p<-rep(NA,NSIMS)
    
    pvals.pred.ti <- rep(NA,NSIMS)
    pvals.nopred.ti<-rep(NA,NSIMS)
    
    for(i in 1:NSIMS){
      
      n.tx <- n[j]/2
      
      #treatment arm 1
      A.tx.1 <- rep(1,n.tx) #trial txd
      X.tx.1 <- rnorm(n.tx,0,sigmax)
      e.tx.1 <- rnorm(n.tx,0,sigmae)
      Y.tx.1.ti <- alpha + 0*A.tx.1 + beta*X.tx.1 + e.tx.1 #no tx effect to test ti
      if(rand == T){
        Y.tx.1.p <- alpha + (gam + rnorm(n.tx, 0, gam*0.1))*A.tx.1 + beta*X.tx.1 + e.tx.1
      } else{
        Y.tx.1.p <- alpha + gam*A.tx.1 + beta*X.tx.1 + e.tx.1
      }
      
      #treatment arm 2
      A.tx.2 <- rep(0,n.tx) #trial controls
      X.tx.2 <- rnorm(n.tx,0,sigmax)
      e.tx.2 <- rnorm(n.tx,0,sigmae)
      Y.tx.2.ti <- alpha + 0*A.tx.2 + beta*X.tx.2 + e.tx.2 #no tx effect to test ti
      if(rand == T){
        Y.tx.2.p <- alpha + (gam + rnorm(n.tx, 0, gam*0.1))*A.tx.2 + beta*X.tx.2 + e.tx.2
      }else{
        Y.tx.2.p <- alpha + gam*A.tx.2 + beta*X.tx.2 + e.tx.2
      }
    
      #testing for tx effect across groups
      Y.tx.p <- c(Y.tx.1.p,Y.tx.2.p)
      Y.tx.ti <- c(Y.tx.1.ti,Y.tx.2.ti)
      X.tx <- c(X.tx.1,X.tx.2)
      A.tx <- c(A.tx.1,A.tx.2)
      
      #---collecting without pred pvals---#
      pvals.nopred.p[i] <- summary(lm(Y.tx.p ~ A.tx))$coef[2,4]
      pvals.nopred.ti[i] <- summary(lm(Y.tx.ti ~ A.tx))$coef[2,4]
      
      #---collecting with pred pvals---#
      pvals.pred.p[i] <- summary(lm(Y.tx.p ~ A.tx + X.tx))$coef[2,4]
      pvals.pred.ti[i] <- summary(lm(Y.tx.ti ~ A.tx + X.tx))$coef[2,4]
      
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
  min.n <- c(min.n, n[which(pow.nopred ==  pow.nopred[pow.nopred >= pow][1])])
  min.n <- c(min.n, n[which(pow.pred ==  pow.pred[pow.pred >= pow][1])])
  
  print(paste("pow.nopred", pow.nopred))
  print(paste("pow.pred", pow.pred))
  print(paste("ti.nopred", ti.nopred))
  print(paste("ti.pred", ti.pred))
  
  #ggplot plots
  pow.p<-c(pow.nopred, pow.pred)
  ti <- c(ti.nopred, ti.pred)
  method <- c(rep("Without Historical Controls",len), rep("With Historical Controls",len))
  
  d.plot<-data.frame( n.plot = rep(n,2), pow.p, method)
  p.p <- ggplot(aes(x=n.plot, y=pow.p, color=method),data=d.plot)+
    geom_line(size=1)+
    coord_cartesian(ylim=c(0,1))+
    ggtitle("Power")+
    xlab("Trial N")+
    theme_classic()+
    theme(text = element_text(size = 14))+
    theme(legend.position = "none")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme( axis.title.y = element_blank())
  
  d.plot <- data.frame( n.plot = rep(n,2), ti, method)
  p.ti <- ggplot(aes(x=n.plot, y=ti, color=method),data=d.plot)+
    geom_line(size=1)+
    coord_cartesian(ylim=c(0,0.5))+
    ggtitle("Type I")+
    xlab("Trial N")+
    theme_classic()+
    theme(legend.position = "none") +
    theme(text = element_text(size = 14))+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme( axis.title.y = element_blank())
  p.ti
  
  return(list(min.n, pow.p, ti, p.p, p.ti))
}