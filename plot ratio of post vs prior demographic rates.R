
#run Graphs.R up to L70 first

library(ncdf4)

#get FIA PFT and total BA by cell (used to filter cells later)

#add up ba across ages
fiaba2<-apply(fiaba,c(1,3:5),FUN=sum)
nplots2<-apply(nplots,2:4,FUN=sum)

#divide through by nplots
fiaba3 <- sweep(fiaba2,2:4,nplots2,FUN="/")
#divide by total BA
relba <- sweep(fiaba3[-1,,,],2:4,fiaba3[1,,,],FUN="/")

#put into PFT x cell matrix
cellrelba <- matrix(nrow=6,ncol=length(cellplotclass))
for (i in 1:length(cellplotclass)) {
  cellrelba[,i] <- relba[,cellplotclass[i]+1,celllonlat[i,1],celllonlat[i,2]]
}

#calc ratios for demographic rates, dropping first row (total)
cellingratio <- cellingpost[-1,]/cellingprior[-1,]
cellmortratio <- cellmortpost[-1,]/cellmortprior[-1,]
cellgrowthratio <- cellgrowthpost[-1,]/cellgrowthprior[-1,]

#remove cells where PFT is <1% of observed BA
cellingratio[cellrelba<0.01] <- NA
cellmortratio[cellrelba<0.01] <- NA
cellgrowthratio[cellrelba<0.01] <- NA

#calc mean and sd of ratios by rate and PFT
mean.ing.ratio <- rowMeans(cellingratio,na.rm=T)
mean.mort.ratio <- rowMeans(cellmortratio,na.rm=T)
mean.growth.ratio <- rowMeans(cellgrowthratio,na.rm=T)
mean.ratios <- cbind(mean.growth.ratio,mean.mort.ratio,mean.ing.ratio)
colnames(mean.ratios) <- c("Growth","Mort","Recr")
rownames(mean.ratios) <- species[-1]

sd.ing.ratio <- apply(cellingratio,1,FUN=sd,na.rm=T)
sd.mort.ratio <- apply(cellmortratio,1,FUN=sd,na.rm=T)
sd.growth.ratio <- apply(cellgrowthratio,1,FUN=sd,na.rm=T)
sd.ratios <- cbind(sd.growth.ratio,sd.mort.ratio,sd.ing.ratio)

#quick plot
plot(as.numeric(mean.ratios),ylim=c(0,2.5),ylab="Ratio",xlab="",xaxt="n")
arrows(1:18,as.numeric(mean.ratios+sd.ratios),
       y1=as.numeric(mean.ratios-sd.ratios),angle=90,length=0.1,code=3)
abline(1,0,lty=2)
abline(v=c(6.5,12.5),lty=1)
text(c(3.5,9.5,15.5),0.2,labels=c("Growth","Longevity","Recruitment"),adj=0.5)

cols<-c("darkblue","skyblue2","darkgreen","green3","red2","pink2")

#x11(6.5,2.5)
pdf("Figure 5 ratio posterior and prior demographic rates.pdf",width=6.5,height=2.5)
par(mfrow=c(1,3),mar=c(0.5,0.5,0.5,0.5),oma=c(2,4,1.5,0))

plot(as.numeric(mean.ratios[1:6]),ylim=c(0.5,2.5),ylab="",xlab="",xaxt="n",
     log="y",cex.axis=1.2,las=1,pch=NA,xlim=c(0.75,6.25))
axis(side=1,at=c(1:6),tick=T,labels=row.names(mean.ratios),cex.axis=1.17)
arrows(0.5,1,6.5,1,lty=2,length=0)
arrows(1:6,as.numeric(mean.ratios[1:6]+sd.ratios[1:6]),
       y1=as.numeric(mean.ratios[1:6]-sd.ratios[1:6]),angle=90,length=0,col=cols)
points(as.numeric(mean.ratios[1:6]),pch=16,col=cols,cex=1.8)

plot(as.numeric(mean.ratios[7:12]),ylim=c(0.5,2.5),ylab="",xlab="",xaxt="n",
     yaxt="n",log="y",pch=NA,xlim=c(0.75,6.25))
axis(side=1,at=c(1:6),tick=T,labels=row.names(mean.ratios),cex.axis=1.17)
axis(side=2,at=c(0.5,1,1.5,2.5),tick=T,labels=F)
arrows(0.5,1,6.5,1,lty=2,length=0)
arrows(1:6,as.numeric(mean.ratios[7:12]+sd.ratios[7:12]),
       y1=as.numeric(mean.ratios[7:12]-sd.ratios[7:12]),angle=90,length=0,col=cols)

points(as.numeric(mean.ratios[7:12]),pch=16,col=cols,cex=1.8)

plot(as.numeric(mean.ratios[13:18]),ylim=c(0.5,2.5),ylab="",xlab="",xaxt="n",
     yaxt="n",log="y",pch=NA,xlim=c(0.75,6.25))
axis(side=1,at=c(1:6),tick=T,labels=row.names(mean.ratios),cex.axis=1.17)
axis(side=2,at=c(0.5,1,1.5,2.5),tick=T,labels=F)
arrows(0.5,1,6.5,1,lty=2,length=0)
arrows(1:6,as.numeric(mean.ratios[13:18]+sd.ratios[13:18]),
       y1=as.numeric(mean.ratios[13:18]-sd.ratios[13:18]),angle=90,length=0,col=cols)
points(as.numeric(mean.ratios[13:18]),pch=16,col=cols,cex=1.8)

#legend("topright",pch=16,col=cols,pt.cex=1.4,c(row.names(mean.ratios)),ncol=3,bty="n")

par(xpd=NA)
mtext("Growth",side=3,line=0.5,cex=1,adj=-2.5)
mtext("Longevity",side=3,line=0.5,cex=1,adj=-1.25)
mtext("Recruitment",side=3,line=0.5,cex=1,adj=0.5)
mtext("Posterior mean / prior mean",side=2,line=33.25,adj=0.5)

dev.off()
