#run code from Graphs.R up to about L156 to generate totdata

library(betareg)

totdata$Ba2 <- totdata$ModelledBa
totdata$G2 <- totdata$GrowthModelled
totdata$M2 <- totdata$MortModelled
totdata$R2 <- totdata$IngModelled
totdata$Age2 <- totdata$Age
totdata$Age3 <- as.factor(totdata$Age)

#create columns for total BA and relative PFT BA
tot.Ba2 <- aggregate(totdata$Ba2,by=list(totdata$Cell,totdata$Age),FUN=sum,na.rm=T)
totdata$TotBa2 <- totdata$StdLogBa2.CM <- tot.Ba2$x[match(paste(totdata$Cell,totdata$Age),paste(tot.Ba2[,1],tot.Ba2[,2]))]
totdata$RelBa2 <- totdata$Ba2/totdata$TotBa2

#create columns for relative PFT growth

#find the maximum PFT growth by cell
max.G2 <- aggregate(totdata$G2,by=list(totdata$Cell),FUN=max,na.rm=T) 
totdata$maxG2 <- max.G2$x[match(totdata$Cell,max.G2[,1])]
#find the next highest PFT growth by cell
max2.G2 <- aggregate(totdata$G2[totdata$G2 < totdata$maxG2],by=list(totdata$Cell[totdata$G2 < totdata$maxG2]),FUN=max,na.rm=T)
totdata$max2G2 <- max2.G2$x[match(totdata$Cell,max2.G2[,1])]

#calculate PFT growth relative to the maximum, or if it is the maximum, then second highest growth
totdata$relG2 <- totdata$G2/ifelse(totdata$G2==totdata$maxG2,totdata$max2G2,totdata$maxG2)

#repeat for mortality now
max.M2 <- aggregate(totdata$M2,by=list(totdata$Cell),FUN=max,na.rm=T)
totdata$maxM2 <- max.M2$x[match(totdata$Cell,max.M2[,1])]
max2.M2 <- aggregate(totdata$M2[totdata$M2 < totdata$maxM2],by=list(totdata$Cell[totdata$M2 < totdata$maxM2]),FUN=max,na.rm=T)
totdata$max2M2 <- max2.M2$x[match(totdata$Cell,max2.M2[,1])]

totdata$relM2 <- totdata$M2/ifelse(totdata$M2==totdata$maxM2,totdata$max2M2,totdata$maxM2)


#repeat for recruitment now
max.R2 <- aggregate(totdata$R2,by=list(totdata$Cell),FUN=max,na.rm=T)
totdata$maxR2 <- max.R2$x[match(totdata$Cell,max.R2[,1])]
max2.R2 <- aggregate(totdata$R2[totdata$R2 < totdata$maxR2],by=list(totdata$Cell[totdata$R2 < totdata$maxR2]),FUN=max,na.rm=T)
totdata$max2R2 <- max2.R2$x[match(totdata$Cell,max2.R2[,1])]

totdata$relR2 <- totdata$R2/ifelse(totdata$R2==totdata$maxR2,totdata$max2R2,totdata$maxR2)

logit <- function(x) log(x/(1-x))
alogit <- function(x) 1/(1+exp(-x))


#create matrix to save regression parameters
save.coef <- matrix(nrow=6,ncol=8)
rownames(save.coef) <- unique(totdata$Pft)
for (iPft in unique(totdata$Pft)) {
  #subset data for 1 PFT where it is at least X% of the total BA
  td <- subset(totdata,Pft==iPft & RelBa2>0.01)
  
  #log BA
  lBa <- (td$RelBa2)
    
  #log and scale demographic rates
  sG <- scale(log(td$relG2))
  sM <- scale(log(td$relM2))
  sR <- scale(log(td$relR2))
  
  sG2 <- predict(betareg(lBa~sG,family=binomial))
  sM2 <- predict(betareg(lBa~sM,family=binomial))
  sR2 <- predict(betareg(lBa~sR,family=binomial))
  
  sG3 <- scale(sG2)
  sM3 <- scale(sM2)
  sR3 <- scale(sR2)
  
  
  #plot individual rates
  layout(matrix(1:3,nrow=1))
  plot(sG,lBa)
  points(sG,sG2,col="red")
  points(sG,predict(loess(lBa~sG)),col="green")
  plot(sM,lBa)
  points(sM,sM2,col="red")
  points(sM,predict(loess(lBa~sM)),col="green")
  plot(sR,lBa)
  points(sR,sR2,col="red")
  points(sR,predict(loess(lBa~sR)),col="green")
  title(main=iPft)
  print(iPft)
  lm.mod <- betareg(lBa ~ (sG + sM + sR)*td$Age2)
  save.coef[iPft,] <- lm.mod$coef$mean
  print(summary(lm.mod))
}

#colnames(save.coef) <- names(lm.mod$coef)

layout(matrix(1:6,nrow=2,byrow=T))
x <- c(2,4,6,8)
for (iPft in unique(totdata$Pft)) {
  plot(x,save.coef[iPft,2]+save.coef[iPft,6]*x,type="b",pch="G",col="red",ylim=c(0.0,1),main=iPft,ylab="Coef",xlab="Age/10")
  points(x,save.coef[iPft,3]+save.coef[iPft,7]*x,type="b",pch="M",col="green")
  points(x,save.coef[iPft,4]+save.coef[iPft,8]*x,type="b",pch="R",col="blue")
}

##########

layout(1)
plot(sM,lBa)
points(sM,alogit(predict(lm(logit(lBa)~sM))),col="red")
points(sM,predict(loess(lBa~sM)),col="green")



plot(sM,logit(lBa))
points(sM,predict(lm(logit(lBa)~sM)),col="red")
points(sM,predict(loess(logit(lBa)~sM)),col="green")
