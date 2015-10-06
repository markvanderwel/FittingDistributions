
library(betareg)

#read and process ncdf output to create totdata
source("create totdata for relative ba regressions.R")

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
#totdata$relG2 <- totdata$G2/totdata$maxG2

#repeat for mortality now
max.M2 <- aggregate(totdata$M2,by=list(totdata$Cell),FUN=max,na.rm=T)
totdata$maxM2 <- max.M2$x[match(totdata$Cell,max.M2[,1])]
max2.M2 <- aggregate(totdata$M2[totdata$M2 < totdata$maxM2],by=list(totdata$Cell[totdata$M2 < totdata$maxM2]),FUN=max,na.rm=T)
totdata$max2M2 <- max2.M2$x[match(totdata$Cell,max2.M2[,1])]

totdata$relM2 <- totdata$M2/ifelse(totdata$M2==totdata$maxM2,totdata$max2M2,totdata$maxM2)
#totdata$relM2 <- totdata$M2/totdata$maxM2


#repeat for recruitment now
max.R2 <- aggregate(totdata$R2,by=list(totdata$Cell),FUN=max,na.rm=T)
totdata$maxR2 <- max.R2$x[match(totdata$Cell,max.R2[,1])]
max2.R2 <- aggregate(totdata$R2[totdata$R2 < totdata$maxR2],by=list(totdata$Cell[totdata$R2 < totdata$maxR2]),FUN=max,na.rm=T)
totdata$max2R2 <- max2.R2$x[match(totdata$Cell,max2.R2[,1])]

totdata$relR2 <- totdata$R2/ifelse(totdata$R2==totdata$maxR2,totdata$max2R2,totdata$maxR2)
#totdata$relR2 <- totdata$R2/totdata$maxR2

logit <- function(x) log(x/(1-x))
alogit <- function(x) 1/(1+exp(-x))


#create matrix to save regression parameters
save.coef <- matrix(nrow=6,ncol=12)
rownames(save.coef) <- unique(totdata$Pft)
for (iPft in unique(totdata$Pft)) {
  #subset data for 1 PFT where it is at least X% of the total BA
  td <- subset(totdata,Pft==iPft & RelBa2>0.01)
  
  #log BA
  lBa <- (td$RelBa2)
  
  #log and scale demographic rates
  sG <- scale(log(td$relG2),scale=T)
  sM <- scale(log(td$relM2),scale=T)
  sR <- scale(log(td$relR2),scale=T)
  
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
  lm.mod <- betareg(lBa ~ (sG + sM + sR)*(td$Age2+(I(td$Age2^2)))) #includes quadratic interactions with Age
  save.coef[iPft,] <- lm.mod$coef$mean
  print(summary(lm.mod))
}

colnames(save.coef) <- names(lm.mod$coef$mean)

#plot standardized regression coefficients by pft and stand age
layout(matrix(1:6,nrow=2,byrow=T))
x <- c(1:10)
for (iPft in unique(totdata$Pft)) {
  plot(x,save.coef[iPft,2]+save.coef[iPft,7]*x+save.coef[iPft,8]*x^2,type="b",pch="G",col="red",ylim=c(0.0,1.2),main=iPft,ylab="Coef",xlab="Age/10")
  points(x,save.coef[iPft,3]+save.coef[iPft,9]*x+save.coef[iPft,10]*x^2,type="b",pch="M",col="green")
  points(x,save.coef[iPft,4]+save.coef[iPft,11]*x+save.coef[iPft,12]*x^2,type="b",pch="R",col="blue")
  
  points(x,save.coef[iPft,1]+save.coef[iPft,5]*x+save.coef[iPft,6]*x^2,type="b",pch="I",col="black")
  
}

#################
###Graph coefficients relative BA vs. demography
#x11(4.5,6.5)
pdf("Figure 7 coefs relative BA vs. demography.pdf",width=4.5,height=6.5)
par(mfrow=c(3,2),mar=c(0.5,0.5,0.5,0.5),oma=c(3.75,4.25,1.75,0.5))

#G=red, M=blue, R=grey

x <- c(1:10)

iPft<-"BC"
plot(x,save.coef[iPft,2]+save.coef[iPft,7]*x+save.coef[iPft,8]*x^2,type="b",pch=21,col="red",bg="red",
     xlim=c(0,10),ylim=c(-0.4,1.3),las=1,cex=1.4,cex.axis=1.15,xaxt="n")
axis(side=1,at=c(0,2,4,6,8,10),labels=F)
points(x,save.coef[iPft,3]+save.coef[iPft,9]*x+save.coef[iPft,10]*x^2,type="b",pch=22,col="blue",
       bg="blue",cex=1.4)
points(x,save.coef[iPft,4]+save.coef[iPft,11]*x+save.coef[iPft,12]*x^2,type="b",pch=24,col="grey60",
       bg="grey60",cex=1.4)
text(0.75,1.25,iPft,cex=1.2)
text(9,0.1,"G",cex=1.2,col="red")
text(2,0,"L",cex=1.2,col="blue")
text(3,0.85,"R",cex=1.2,col="grey60")

iPft<-"BH"
plot(x,save.coef[iPft,2]+save.coef[iPft,7]*x+save.coef[iPft,8]*x^2,type="b",pch=21,col="red",bg="red",
     xlim=c(0,10),ylim=c(-0.4,1.3),las=1,cex=1.4,cex.axis=1.15,xaxt="n",yaxt="n")
axis(side=1,at=c(0,2,4,6,8,10),labels=F)
axis(side=2,at=seq(0.0,1.5,0.5),labels=F)
points(x,save.coef[iPft,3]+save.coef[iPft,9]*x+save.coef[iPft,10]*x^2,type="b",pch=22,col="blue",
       bg="blue",cex=1.4)
points(x,save.coef[iPft,4]+save.coef[iPft,11]*x+save.coef[iPft,12]*x^2,type="b",pch=24,col="grey60",
       bg="grey60",cex=1.4)
text(0.75,1.25,iPft,cex=1.2)
text(9,0.8,"G",cex=1.2,col="red")
text(7,0.15,"L",cex=1.2,col="blue")
text(2,1,"R",cex=1.2,col="grey60")

iPft<-"NC"
plot(x,save.coef[iPft,2]+save.coef[iPft,7]*x+save.coef[iPft,8]*x^2,type="b",pch=21,col="red",bg="red",
     xlim=c(0,10),ylim=c(-0.4,1.3),las=1,cex=1.4,cex.axis=1.15,xaxt="n")
axis(side=1,at=c(0,2,4,6,8,10),labels=F)
points(x,save.coef[iPft,3]+save.coef[iPft,9]*x+save.coef[iPft,10]*x^2,type="b",pch=22,col="blue",
       bg="blue",cex=1.4)
points(x,save.coef[iPft,4]+save.coef[iPft,11]*x+save.coef[iPft,12]*x^2,type="b",pch=24,col="grey60",
       bg="grey60",cex=1.4)
text(0.75,1.25,iPft,cex=1.2)
text(9,0.75,"G",cex=1.2,col="red")
text(6,0.7,"L",cex=1.2,col="blue")
text(8,0.1,"R",cex=1.2,col="grey60")

iPft<-"NH"
plot(x,save.coef[iPft,2]+save.coef[iPft,7]*x+save.coef[iPft,8]*x^2,type="b",pch=21,col="red",bg="red",
     xlim=c(0,10),ylim=c(-0.4,1.3),las=1,cex=1.4,cex.axis=1.15,xaxt="n",yaxt="n")
axis(side=1,at=c(0,2,4,6,8,10),labels=F)
axis(side=2,at=seq(0.0,1.5,0.5),labels=F)
points(x,save.coef[iPft,3]+save.coef[iPft,9]*x+save.coef[iPft,10]*x^2,type="b",pch=22,col="blue",
       bg="blue",cex=1.4)
points(x,save.coef[iPft,4]+save.coef[iPft,11]*x+save.coef[iPft,12]*x^2,type="b",pch=24,col="grey60",
       bg="grey60",cex=1.4)
text(0.75,1.25,iPft,cex=1.2)
text(7,1.2,"G",cex=1.2,col="red")
text(8,-0.05,"L",cex=1.2,col="blue")
text(1,0.75,"R",cex=1.2,col="grey60")

iPft<-"SC"
plot(x,save.coef[iPft,2]+save.coef[iPft,7]*x+save.coef[iPft,8]*x^2,type="b",pch=21,col="red",bg="red",
     xlim=c(0,10),ylim=c(-0.4,1.3),las=1,cex=1.4,cex.axis=1.15,xaxt="n")
axis(side=1,at=c(0,2,4,6,8,10),labels=c(0,20,40,60,80,100),cex.axis=1.15)
points(x,save.coef[iPft,3]+save.coef[iPft,9]*x+save.coef[iPft,10]*x^2,type="b",pch=22,col="blue",
       bg="blue",cex=1.4)
points(x,save.coef[iPft,4]+save.coef[iPft,11]*x+save.coef[iPft,12]*x^2,type="b",pch=24,col="grey60",
       bg="grey60",cex=1.4)
text(0.75,1.25,iPft,cex=1.2)
text(8,1.2,"G",cex=1.2,col="red")
text(5,-0.1,"L",cex=1.2,col="blue")
text(1,0.9,"R",cex=1.2,col="grey60")

iPft<-"SH"
plot(x,save.coef[iPft,2]+save.coef[iPft,7]*x+save.coef[iPft,8]*x^2,type="b",pch=21,col="red",bg="red",
     xlim=c(0,10),ylim=c(-0.4,1.3),las=1,cex=1.4,cex.axis=1.15,xaxt="n",yaxt="n")
axis(side=1,at=c(0,2,4,6,8,10),labels=c(0,20,40,60,80,100),cex.axis=1.15)
axis(side=2,at=seq(0.0,1.5,0.5),labels=F)
points(x,save.coef[iPft,3]+save.coef[iPft,9]*x+save.coef[iPft,10]*x^2,type="b",pch=22,col="blue",
       bg="blue",cex=1.4)
points(x,save.coef[iPft,4]+save.coef[iPft,11]*x+save.coef[iPft,12]*x^2,type="b",pch=24,col="grey60",
       bg="grey60",cex=1.4)
text(0.75,1.25,iPft,cex=1.2)
text(3,0.45,"G",cex=1.2,col="red")
text(7,0.22,"L",cex=1.2,col="blue")
text(2,0.9,"R",cex=1.2,col="grey60")

par(xpd=NA)
legend(-11.5,5.62,c("Growth","Longevity","Recruitment"),
       lty=1,pch=c(21,22,24),cex=1.3,pt.cex=1.4,col=c("red","blue","grey60"),
       pt.bg=c("red","blue","grey60"),bty="n",ncol=3)
mtext("Standardized effect", side=2, line=2.5, outer=TRUE, adj=0.5, cex=1)
mtext("Age (y)", side=1, line=2, outer=TRUE, adj=0.52, cex=1)

dev.off()

##########
layout(1)
#plot successional changes for average demographic parameters
matplot(t(alogit(save.coef[,1]+outer(save.coef[,5],x)+outer(save.coef[,6],x^2))),type="b",ylab="Mean relative BA")

#################
###Graph mean relative BA per PFT vs age

#x11(3.8,3.5)
pdf("Figure relative BA vs. age.pdf",width=3.8,height=3.5)
par(mar=c(0.5,0.5,2,0.5),oma=c(3,3,0,0))

x <- c(1:10)

plot(x,t(alogit(save.coef["BC",1]+outer(save.coef["BC",5],x)+outer(save.coef["BC",6],x^2))),
        type="b",pch=16, col="darkblue", ylim=c(0.03,0.4),xaxt="n",las=1,cex.axis=0.9)
axis(side=1,at=c(2,4,6,8,10),labels=c(20,40,60,80,100),cex.axis=1)
lines(x,t(alogit(save.coef["BH",1]+outer(save.coef["BH",5],x)+outer(save.coef["BH",6],x^2))),
        type="b",pch=16, col="skyblue2")
lines(x,t(alogit(save.coef["NC",1]+outer(save.coef["NC",5],x)+outer(save.coef["NC",6],x^2))),
      type="b",pch=16, col="darkgreen")
lines(x,t(alogit(save.coef["NH",1]+outer(save.coef["NH",5],x)+outer(save.coef["NH",6],x^2))),
      type="b",pch=16, col="green3")
lines(x,t(alogit(save.coef["SC",1]+outer(save.coef["SC",5],x)+outer(save.coef["SC",6],x^2))),
      type="b",pch=16, col="red2")
lines(x,t(alogit(save.coef["SH",1]+outer(save.coef["SH",5],x)+outer(save.coef["SH",6],x^2))),
      type="b",pch=16, col="pink2")

par(xpd=NA)
legend(0.9,0.495,c("BC","BH","NC","NH","SC","SH"),
       col=c("darkblue","skyblue2","darkgreen","green3","red2","pink2"),ncol=3,bty="n",
       pch=16,cex=0.9,pt.cex=1,lty=5)

mtext("Mean relative basal area", side=2, line=2, outer=TRUE, adj=0.45, cex=1)
mtext("Age (y)", side=1, line=1.5, outer=TRUE, adj=0.58, cex=1)

dev.off()
