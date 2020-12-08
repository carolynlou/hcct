powercomp_adnicont <- function(NSIMS, n = seq(50, 150, by = 50), gam, pow = 0.8, data, boot = F){
  len <- length(n)
  
  pow.pred<-rep(NA,len)
  pow.nopred<-rep(NA,len)
  pow.lin <- rep(NA,len)
  
  ti.pred <- rep(NA,len)
  ti.nopred <- rep(NA,len)
  ti.lin <- rep(NA,len)
  
  set.seed(12345)
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
      ftx.1 <- tx1[sample(nrow(tx1), n.tx, replace=boot),] #"bootstrap" a future tx sample of MCI pts, tx arm 1
      ftx.1 %<>% mutate(A = rep(1, n.tx),
                        Y.p = mem.ch + A*(gam + rnorm(n.tx, 0, gam*0.1)),
                        Y.t1 = mem.ch)
      
      #treatment arm 2
      ftx.2 <- tx2[sample(nrow(tx2), n.tx, replace=boot),] #a future tx sample of MCI pts, tx arm 2
      ftx.2 %<>% mutate(A = rep(0, n.tx),
                        Y.p = mem.ch + A*(gam + rnorm(n.tx, 0, gam*0.1)),
                        Y.t1 = mem.ch)
      
      dat.final = rbind(ftx.1, ftx.2)
      
      #collecting pvals
      pvals.nopred.p[i] <- summary(lm(Y.p ~ A + age + memb, data = dat.final))$coef[2,4]
      pvals.nopred.ti[i] <- summary(lm(Y.t1 ~ A + age + memb, data = dat.final))$coef[2,4]
      
      pvals.pred.p[i] <- summary(lm(Y.p ~ A + SPARE.AD + age + memb, data = dat.final))$coef[2,4]
      pvals.pred.ti[i] <- summary(lm(Y.t1 ~ A + SPARE.AD + age + memb, data = dat.final))$coef[2,4]
      
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
  
  return(list(min.n, p.p, p.ti))
}
