powercomp_surv_adni <- function(NSIMS, n, gam, pow = 0.8, data){
  len <- length(n)
  data <- dat
  
  pow.pred<-rep(NA,len)
  pow.nopred<-rep(NA,len)
  
  ti.pred <- rep(NA,len)
  ti.nopred <- rep(NA,len)
  
  cens <- rep(NA, len)
  
  set.seed(960113)
  start<- proc.time()
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
      ftx.1 = ftx.1 %>%
        mutate(A.tx = rep(1,n.tx),#trial tx'd
               Y.tx.p = pmin(TimeToConv*(gam + rnorm(n.tx, 0, gam*0.1)), 36),
               cens.p = 1*(Y.tx.p < 36),
               Y.tx.ti = TimeToConv,
               cens.ti = 1*(Y.tx.ti < 36))
      
      #treatment arm 2, ctl
      ftx.2 <- tx2[sample(nrow(tx2),n.tx,replace=F),] #a future tx sample of MCI pts, tx arm 2
      ftx.2 = ftx.2 %>%
        mutate(A.tx = rep(0,n.tx),#trial tx'd
               Y.tx.p = TimeToConv,
               cens.p = 1*(Y.tx.p < 36),
               Y.tx.ti = TimeToConv,
               cens.ti = 1*(Y.tx.ti < 36))
      
      #collecting pvals
      dat.final = rbind(ftx.1, ftx.2)
      
      pvals.nopred.p[i] <- summary(coxph(Surv(time = Y.tx.p, cens.p) ~ A.tx + AGE + memb, data = dat.final))$coef[1,5]
      pvals.pred.p[i] <- summary(coxph(Surv(time = Y.tx.p, event = cens.p) ~ A.tx + SPARE.AD + AGE + memb, data = dat.final))$coef[1,5]
      
      pvals.nopred.ti[i] <- summary(coxph(Surv(time = Y.tx.ti, event = cens.ti) ~ A.tx + AGE + memb, data = dat.final))$coef[1,5]
      pvals.pred.ti[i] <- summary(coxph(Surv(time = Y.tx.ti, event = cens.ti) ~ A.tx + SPARE.AD + AGE + memb, data = dat.final))$coef[1,5]
      
      cens.ti.perc[i] <- sum(dat.final$cens.ti)/length(dat.final$cens.ti)
      
      setTxtProgressBar(pb, (j-1)*NSIMS+i)
    }
    cens[j] <- 1-mean(cens.ti.perc)
    pow.pred[j] <- sum(pvals.pred.p < 0.05)/NSIMS 
    pow.nopred[j] <- sum(pvals.nopred.p < 0.05)/NSIMS 
    
    ti.pred[j] <- sum(pvals.pred.ti < 0.05)/NSIMS 
    ti.nopred[j] <- sum(pvals.nopred.ti < 0.05)/NSIMS 
  }
  close(pb)
  proc.time() - start
  
  pow.pred
  pow.nopred
  
  ti.pred
  ti.nopred
  
  min.n <- NULL
  min.n <- c(min.n, n[which(pow.nopred ==  pow.nopred[pow.nopred >= pow][1])]) #no HCs
  min.n <- c(min.n, n[which(pow.pred ==  pow.pred[pow.pred >= pow][1])]) #with HCs
  
  print(paste("pow.nopred", pow.nopred))
  print(paste("pow.pred", pow.pred))
  print(paste("ti.nopred", ti.nopred))
  print(paste("ti.pred", ti.pred))
  print(paste("perc cens", cens))
  
  #ggplot plots
  pow.p<-c(pow.nopred, pow.pred)
  ti <- c(ti.nopred, ti.pred)
  method <- c(rep("Without Historical Controls",len), rep("With Historical Controls",len))
  
  d.plot<-data.frame( n.plot = rep(n,2), pow.p, method)
  p.p <- ggplot(aes(x=n.plot, y=pow.p, color=method),data=d.plot)+
    geom_line(size=1)+
    coord_cartesian(ylim=c(0,1))+
    ggtitle("Power")+
    #ylab("Power")+
    xlab("Trial N")+
    theme_classic()+
    theme(text = element_text(size = 14))+
    theme(legend.position = "none")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme( axis.title.y = element_blank())
  # theme(legend.justification = c(0.05, 1), legend.position = c(0.05, 1)) + 
  # scale_color_discrete(name = "Method", 
  #                      labels = c("With Historical Controls", "Without Historical Controls")) 
  
  d.plot<-data.frame( n.plot = rep(n,2), ti, method)
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
  #theme(legend.justification = c(1, 1), legend.position = c(1, 0.75)) + 
  #scale_color_discrete(name = "Method", 
  #                     labels = c("With Historical Controls", "Without Historical Controls")) 
  
  p.ti
  return(list(min.n, p.p, p.ti)) 
}
