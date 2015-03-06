
cols=brewer.pal(6,"Dark2")
cols.light=rgb(t(col2rgb(cols)),alpha=128,maxColorValue=255)


x11(15,10)
layout(matrix(1:6,nrow=2,byrow=F))

for (iDemo in demo) {
  pardata.sub <- subset(pardata,Demo==iDemo)
  
  #loess.t <- with(pardata.sub,predict(loess(LogTrendSurfaceMean~Temp)))
  with(pardata.sub,plot(Temp,(LogTrendSurfaceMean),col=cols[Pft],main=paste("Post",iDemo,sep=" ")))  
  with(pardata.sub,plot(Temp,log(Prior),col=cols[Pft],main=paste("Prior",iDemo,sep=" ")))

}


x11(15,10)
layout(matrix(1:6,nrow=2,byrow=F))

for (iDemo in demo) {
  pardata.sub <- subset(pardata,Demo==iDemo)
  
  #loess.t <- with(pardata.sub,predict(loess(LogTrendSurfaceMean~Temp)))  
  with(pardata.sub,plot(Precip,(LogTrendSurfaceMean),col=cols[Pft],main=paste("Post",iDemo,sep=" ")))
  with(pardata.sub,plot(Precip,log(Prior),col=cols[Pft],main=paste("Prior",iDemo,sep=" ")))

}


x11(15,15)
par(mar=c(2,2,0,0))
layout(matrix(1:18,nrow=3,byrow=T))

y.lim = list(
  G = c(exp(-4),exp(1.4)),
  M = c(exp(2),exp(6.2)),
  R = c(exp(-2),exp(6.2))
  )

x.lim = list(
  BC = c(2,9),
  BH = c(2,10),
  NC = c(2,14),
  NH = c(2,20),
  SC = c(7,23),
  SH = c(6,22)
  )

col.range.fade <- function(pft,vals){
  orig.col = cols[pft]
  fade.col = rgb(t(col2rgb(orig.col)),alpha=20,maxColorValue=255)
  #fade.col="white"
  #print(ifelse(vals>x.lim[[pft]][1] & vals<x.lim[[pft]][2],orig.col,fade.col))
  return (ifelse(vals>x.lim[[pft]][1] & vals<x.lim[[pft]][2],orig.col,fade.col))
}

#temperature graphs
for (iDemo in demo) {
  for (iSp in species[-1]) {
  
    pardata.sub <- subset(pardata,Demo==iDemo & Pft==iSp)
      
    with(pardata.sub,plot(Temp,(Prior),col=col.range.fade(which(iSp==species[-1]),Temp),main=paste("\n",iSp,iDemo,sep=" "),ylim=y.lim[[iDemo]],xlab="",ylab="",log="y"))
    #with(pardata.sub,points(Temp,exp(LogTrendSurfaceMean),col=rgb(0,0,0,0.1),main=paste("Post",iDemo,sep=" ")))  
    with(pardata.sub,points(Temp,exp(LogValue),col=rgb(0,0,0,0.3),main=paste("Post",iDemo,sep=" ")))  
   
    
    indlatlon <- unique(pardata.sub[,3:4])
    pardata.sub$Index <- mapply(FUN=function(x1,x2) 
                                {which(x1==indlatlon[,1] & x2==indlatlon[,2])},
                                pardata.sub$Lat,
                                pardata.sub$Lon)
    
    pd.horz <- data.frame(
      Index = unique(pardata.sub$Index),
      LogValueMin = tapply(pardata.sub$LogValue,pardata.sub$Index,FUN=min),
      LogValueMed = tapply(pardata.sub$LogValue,pardata.sub$Index,FUN=median),
      LogValueMax = tapply(pardata.sub$LogValue,pardata.sub$Index,FUN=max),
      Temp = tapply(pardata.sub$Temp,pardata.sub$Index,FUN=mean),
      Precip = tapply(pardata.sub$Precip,pardata.sub$Index,FUN=mean)
    )
    
    class.cols = c("red","darkgreen","blue")
    
    for (iClass in 1:3) {
      p <- exp(predict(loess(pd.horz[,iClass+1] ~ pd.horz$Temp)))
      lines(sort(pd.horz$Temp),p[order(pd.horz$Temp)],col=class.cols[iClass],lwd=2)
      
    } #for pd.horz classes
    
    
  }
}


#precip graphs
x11(15,15)
par(mar=c(2,2,0,0))
layout(matrix(1:18,nrow=3,byrow=T))

for (iDemo in demo) {
  for (iSp in species[-1]) {
    
    pardata.sub <- subset(pardata,Demo==iDemo & Pft==iSp)
    
    with(pardata.sub,plot(Precip,log(Prior),col=col.range.fade(which(iSp==species[-1]),Temp),main=paste("\n",iSp,iDemo,sep=" "),ylim=y.lim[[iDemo]],xlab="",ylab=""))
    with(pardata.sub,points(Precip,(LogTrendSurfaceMean),col=rgb(0,0,0,0.1),main=paste("Post",iDemo,sep=" ")))  
    
  }
}




