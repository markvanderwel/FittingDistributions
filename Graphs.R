##GRAPHS

library(ncdf4)
library(colorRamps)
library(maps)
library(RColorBrewer)
library(filzbach)

nc <- nc_open("FittingOutput.nc", readunlim=F)
#print(nc)

nc.orig <- nc_open("SimOutWithClimData.nc", readunlim=F)

ba.orig <- ncvar_get(nc.orig,"Ba")

ba<-ncvar_get(nc,"Ba")
modelledba<-ba #modelledba<-ncvar_get(nc,"ModelledBa")
fiaba<-ncvar_get(nc,"FiaBa")

nplots <- ncvar_get(nc,"FiaNPlots")

lat<-ncvar_get(nc,"Lat")
lon<-ncvar_get(nc,"Lon")
plotclass<-ncvar_get(nc,"PlotClass")
species<-c("Tot","BC","BH","NC","NH","SC","SH")  #ncvar_get(nc,"Species")
Time<- 0:15 #ncvar_get(nc,"Time")

temperature<-ncvar_get(nc,"Temperature")
precipitation<-ncvar_get(nc,"Precipitation")

celltemperature<-ncvar_get(nc,"CellTemperature")
cellprecipitation<-ncvar_get(nc,"CellPrecip")

cellplotclass<-ncvar_get(nc,"CellPlotClass")

cellingprior <- ncvar_get(nc,"CellIngPrior")
cellmortprior <- ncvar_get(nc,"CellMortPrior")
cellgrowthprior <- ncvar_get(nc,"CellGrowthPrior")

cellingpost <- ncvar_get(nc,"CellIngPost")
cellmortpost <- ncvar_get(nc,"CellMortPost")
cellgrowthpost <- ncvar_get(nc,"CellGrowthPost")

cellmodelleding <- cellingpost #cellmodelleding <- ncvar_get(nc,"CellModelledIng")
cellmodelledmort <- cellmortpost #cellmodelledmort <- ncvar_get(nc,"CellModelledMort")
cellmodelledgrowth <- cellgrowthpost #cellmodelledgrowth <- ncvar_get(nc,"CellModelledGrowth")

cellingupper95 <- ncvar_get(nc,"CellIngUpper95")
cellmortupper95 <- ncvar_get(nc,"CellMortUpper95")
cellgrowthupper95 <- ncvar_get(nc,"CellGrowthUpper95")

cellinglower95 <- ncvar_get(nc,"CellIngLower95")
cellmortlower95 <- ncvar_get(nc,"CellMortLower95")
cellgrowthlower95 <- ncvar_get(nc,"CellGrowthLower95")

nc_close(nc)

###############################################

celllonlat<- matrix(nrow=length(celltemperature),ncol=2)
for (i in 1:length(celltemperature)) {
  celllonlat[i,1:2] <- which(temperature==celltemperature[i] & precipitation == cellprecipitation[i],arr.ind=T)
}

lonlatclasscell<-array(dim=c(length(lon),length(lat),3),dimnames=c("lon","lat","class"))
for(i in 1:nrow(celllonlat)) {
  lonlatclasscell[celllonlat[i,1],celllonlat[i,2],cellplotclass[i]+1]<-i
}

df.cell <- integer(0)
df.age<- integer(0)
df.lon <- numeric(0)
df.lat <- numeric(0)
df.pft <- character(0)
df.plotclass <- numeric(0)
df.fiaba <- numeric(0)
df.origba <- numeric(0)
df.bestba <- numeric(0)
df.modelledba <- numeric(0)
df.growthprior <- numeric(0)
df.mortprior <- numeric(0)
df.ingprior <- numeric(0)
df.growthbest <- numeric(0)
df.mortbest <- numeric(0)
df.ingbest <- numeric(0)
df.growthmodelled <- numeric(0)
df.mortmodelled <- numeric(0)
df.ingmodelled <- numeric(0)
df.growthupper95 <- numeric(0)
df.mortupper95 <- numeric(0)
df.ingupper95 <- numeric(0)
df.growthlower95 <- numeric(0)
df.mortlower95 <- numeric(0)
df.inglower95 <- numeric(0)
df.celltemperature <- numeric(0)
df.cellprecipitation <- numeric(0)

i <- 0

for (iLon in 1:31) {
  for (iLat in 1:27) {
    for (c in 1:3) {  #plotclass
      
      if (!is.na(lonlatclasscell[iLon,iLat,c])) {    
        for (t in 1:16) {    
    
          if (nplots[t,c,iLon,iLat] > 9) {         
            for (s in 2:7) {      
              
              i <- i + 1
              
              df.cell[i] <- lonlatclasscell[iLon,iLat,c]
              df.age[i] <- Time[t]
              df.lon[i] <- lon[iLon]
              df.lat[i] <- lat[iLat]
              df.plotclass[i] <- plotclass[c]
              df.pft[i] <- species[s]
              df.fiaba[i] <- fiaba[s,t,c,iLon,iLat]
              df.origba[i] <- ba.orig[t*2,s,1,iLon,iLat]
              df.bestba[i] <- ba[t,s,c,iLon,iLat]
              df.modelledba[i] <- modelledba[t,s,c,iLon,iLat]
              df.growthprior[i] <- cellgrowthprior[s,df.cell[i]]
              df.mortprior[i] <- cellmortprior[s,df.cell[i]]
              df.ingprior[i] <- cellingprior[s,df.cell[i]]
              df.growthbest[i] <- cellgrowthpost[s,df.cell[i]]
              df.mortbest[i] <- cellmortpost[s,df.cell[i]]
              df.ingbest[i] <- cellingpost[s,df.cell[i]]
              df.growthmodelled[i] <- cellmodelledgrowth[s,df.cell[i]]
              df.mortmodelled[i] <- cellmodelledmort[s,df.cell[i]]
              df.ingmodelled[i] <- cellmodelleding[s,df.cell[i]]
              df.growthupper95[i] <- cellgrowthupper95[s,df.cell[i]]
              df.mortupper95[i] <- cellmortupper95[s,df.cell[i]]
              df.ingupper95[i] <- cellingupper95[s,df.cell[i]]
              df.growthlower95[i] <- cellgrowthlower95[s,df.cell[i]]
              df.mortlower95[i] <- cellmortlower95[s,df.cell[i]]
              df.inglower95[i] <- cellinglower95[s,df.cell[i]]
              df.celltemperature[i] <- celltemperature[df.cell[i]]
              df.cellprecipitation[i] <- cellprecipitation[df.cell[i]]
              
              
            } #s
          } #nPlots>9
        } #c
      } #t
    } #not water
  } #lat
} #lon
totdata <- data.frame(Cell=df.cell,Age=df.age,Lon=df.lon,Lat=df.lat,PlotClass=df.plotclass,
                   Pft=df.pft,FiaBa=df.fiaba,OrigBa=df.origba,BestBa=df.bestba,
                   ModelledBa=df.modelledba,GrowthPrior=df.growthprior,MortPrior=df.mortprior,
                   IngPrior=df.ingprior,GrowthBest=df.growthbest,MortBest=df.mortbest,
                   IngBest=df.ingbest,GrowthModelled=df.growthmodelled,
                   MortModelled=df.mortmodelled,IngModelled=df.ingmodelled,
                   GrowthUpper=df.growthupper95,MortUpper=df.mortupper95,IngUpper=df.ingupper95,
                   GrowthLower=df.growthlower95,MortLower=df.mortlower95,IngLower=df.inglower95,
                   Temperature=df.celltemperature,Precipitation=df.cellprecipitation)

#Add color code to dataset.
totdata$col<-ifelse(totdata$Pft=="BC","darkblue",
                    ifelse(totdata$Pft=="BH","skyblue2",
                           ifelse(totdata$Pft=="NC","darkgreen",
                                  ifelse(totdata$Pft=="NH","lightgreen",
                                         ifelse(totdata$Pft=="SC","darkred","red")))))

#Add unique code for cell
totdata$Lat_Lon<-paste(totdata$Lat,totdata$Lon,sep="_")

##############################################################
#FIGURE: Goodness of fit: observed vs predicted PFT basal area
#Now per functional type.

x11(6.5,8)
par(mfrow=c(3,2),mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,0.1,0.1))

plot(totdata$ModelledBa[totdata$Pft=="BC"],totdata$FiaBa[totdata$Pft=="BC"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(2,28,"BC",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="BH"],totdata$FiaBa[totdata$Pft=="BH"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",yaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
axis(side=2,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(2,28,"BH",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="NC"],totdata$FiaBa[totdata$Pft=="NC"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(2,28,"NC",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="NH"],totdata$FiaBa[totdata$Pft=="NH"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",yaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
axis(side=2,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(2,28,"NH",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="SC"],totdata$FiaBa[totdata$Pft=="SC"],
     xlim=c(0,30),ylim=c(0,30),las=1,col="grey60",pch=16,cex=1.4,cex.axis=1.2)
abline(0,1,lwd=2,lty=2)
text(2,28,"SC",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="SH"],totdata$FiaBa[totdata$Pft=="SH"],
     xlim=c(0,30),ylim=c(0,30),las=1,yaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(2,28,"SH",cex=1.4)

mtext(expression("Observed basal area (m"^2*" ha"^-1*")"),
      side=2, line=2, outer=TRUE, adj=0.5, cex=1)
mtext(expression("Predicted basal area (m"^2*" ha"^-1*")"),
      side=1, line=2.5, outer=TRUE, adj=0.5, cex=1)

##############################################################
#FIGURE: Observed vs. predicted size distribution
##Compare modelled diameter distributions to FIA data
##Use same ages as for BA maps?

library(ncdf4)
nc <- nc_open("FittingOutput.nc", readunlim=F)

fiadbh<-ncvar_get(nc,"FiaDbhDens")
modeldbh<-ncvar_get(nc,"DbhDens")
modeldbh[modeldbh==-999] <- NA
fian<-ncvar_get(nc,"FiaNPlots")
dbh<-ncvar_get(nc,"Dbh")

fiadbh2 <- matrix(nrow=10,ncol=21)
modeldbh2 <- fiadbh2

fian2 <- apply(fian,1,sum)      #sum over time

for (i in 1:21) {
  
  #weighted mean (by age) per age class
  fiadbh2[,i] <- apply(fiadbh[,i,,,],1,FUN=weighted.mean,w=fian[i,,,],na.rm=T)    
  modeldbh2[,i] <- apply(modeldbh[,i,,,],1,FUN=weighted.mean,w=fian[i,,,],na.rm=T)
}

bins <- numeric(21)
bins[1:5] <- 1
bins[6:7] <- 2
bins[8:21] <- 3

fiadbh3 <- matrix(nrow=10,ncol=4)
modeldbh3 <- fiadbh2

for (i in 1:3) {
  
  fiadbh3[,i] <- apply(fiadbh2[,bins==i],1,FUN=weighted.mean,w=fian2[bins==i],na.rm=T)
  modeldbh3[,i] <- apply(modeldbh2[,bins==i],1,FUN=weighted.mean,w=fian2[bins==i],na.rm=T)
}
fiadbh3[,4] <- apply(fiadbh2,1,FUN=weighted.mean,w=fian2,na.rm=T)
modeldbh3[,4] <- apply(modeldbh2,1,FUN=weighted.mean,w=fian2,na.rm=T)


#x11(15,15)
#layout(matrix(1:4,nrow=2,ncol=2,byrow=T))

x11(8,3.5)
par(mfrow=c(1,4),mar=c(2,2,2,0.5),oma=c(3,3,0.5,0.5))

binlabs <- c("0-50 y","60-70 y",">80 y","All ages")

for (i in 1:4) {
  
  dbhcomp <- cbind(fiadbh3[,i],modeldbh3[,i])
  
  matplot(dbh-5,
          dbhcomp,
          type="b",
          pch=19,
          lty=1,
          col=c("black","grey"),
          log="y",
          cex=1.2,
          ylim=c(4e-2,5e3),
          xlim=c(5,95),
          xlab="DBH midpoint (cm)",
          ylab="Density (/ha)",
          main=binlabs[i]
  )
  
}

legend("topright",
       c("Observed (FIA)","Predicted (CAIN)"),
       lty=1,
       pch=19,
       col=c("black","grey60"),
       bty="n"
)

mtext("Density (/ha)",side=2, line=1.5, outer=TRUE, adj=0.5, cex=0.8)
mtext("DBH midpoint (cm)",side=1, line=1.5, outer=TRUE, adj=0.5, cex=0.8)

############################################################################
############################################################################
#FIGURE: Predicted PFT BA against predicted potential growth, mortality, and
#recruitment among and within grid cells

tot.slope<-read.table("Overall slopes.txt",h=T)
cell.slope<-read.table("Cell slopes.txt",h=T)

#select age
data<-totdata[totdata$Age==5,]
data$Lat_Lon<-factor(data$Lat_Lon)
gdata<-data

x11(width=6.5,height=8)
par(mfrow=c(3,2),mar=c(1.5,4.5,1,1),oma=c(3,0,2,0))

plot(gdata$ModelledBa,gdata$GrowthBest,pch=NA,las=1,cex.axis=1.2)
for (cell in unique(gdata$Lat_Lon)){
      celldata<-gdata[gdata$Lat_Lon==cell,]
    for (pft in unique(celldata$Pft)){
      pftdata<-celldata[celldata$Pft==pft,]
      pftdata<-pftdata[order(pftdata$ModelledBa),]
      points(pftdata$ModelledBa,pftdata$GrowthBest,col=pftdata$col,pch=21,cex=1.1)
      lines(pftdata$ModelledBa,pftdata$GrowthBest,col=pftdata$col,lwd=2)
   }
}

xvals.c<-c(1.2,2.2,3.2,4.2,5.2,6.2)
xvals.t<-c(0.8,1.8,2.8,3.8,4.8,5.8)
yvals.c<-c(cell.slope[cell.slope$DemoRate=="GrowthBest",]$mean.slope)
yvals.t<-c(tot.slope[tot.slope$DemoRate=="GrowthBest",]$slope50)
y_upper.c<-c(cell.slope[cell.slope$DemoRate=="GrowthBest",]$max.slope)
y_upper.t<-c(tot.slope[tot.slope$DemoRate=="GrowthBest",]$slope97.5)
y_lower.c<-c(cell.slope[cell.slope$DemoRate=="GrowthBest",]$min.slope)
y_lower.t<-c(tot.slope[tot.slope$DemoRate=="GrowthBest",]$slope2.5)
plot(xvals.t,yvals.t,pch=21,col="black",bg="black",xaxt="n",ylim=c(-0.1,0.3),
     xlim=c(0.5,6.5),las=1,cex=1.5,cex.axis=1.2)
abline(0,0,lty=2)
arrows(xvals.c,y_lower.c,xvals.c,y_upper.c,code=3,length=0)
arrows(xvals.t,y_lower.t,xvals.t,y_upper.t,code=3,length=0)
points(xvals.c,yvals.c,pch=21,cex=1.5,col="black",bg="white")
axis(side=1,at=c(1:6),tick=T,labels=c("BC","BH","NC","NH","SC","SH"),cex.axis=1.2)
legend("topleft",c("Among grid cells","Within grid cells"),pch=21,
       pt.bg=c("black","white"),bty="n",pt.cex=1.5,cex=1.3)

#Predicted potential mortality

#select age
data<-totdata2[totdata2$Age==5,]
mdata<-data

plot(mdata$ModelledBa,mdata$MortBest,pch=NA,las=1,cex.axis=1.2)
for (cell in unique(mdata$Lat_Lon)){
  celldata<-mdata[mdata$Lat_Lon==cell,]
  for (pft in unique(celldata$Pft)){
    pftdata<-celldata[celldata$Pft==pft,]
    pftdata<-pftdata[order(pftdata$ModelledBa),]
    points(pftdata$ModelledBa,pftdata$MortBest,col=pftdata$col,pch=21,cex=1.1)
    lines(pftdata$ModelledBa,pftdata$MortBest,col=pftdata$col,lwd=2)
  }
}

yvals.c<-c(cell.slope[cell.slope$DemoRate=="MortBest",]$mean.slope)
yvals.t<-c(tot.slope[tot.slope$DemoRate=="MortBest",]$slope50)
y_upper.c<-c(cell.slope[cell.slope$DemoRate=="MortBest",]$max.slope)
y_upper.t<-c(tot.slope[tot.slope$DemoRate=="MortBest",]$slope97.5)
y_lower.c<-c(cell.slope[cell.slope$DemoRate=="MortBest",]$min.slope)
y_lower.t<-c(tot.slope[tot.slope$DemoRate=="MortBest",]$slope2.5)
plot(xvals.t,yvals.t,pch=21,col="black",bg="black",xaxt="n",ylim=c(-30,60),
     xlim=c(0.5,6.5),las=1,cex=1.5,cex.axis=1.2)
abline(0,0,lty=2)
arrows(xvals.c,y_lower.c,xvals.c,y_upper.c,code=3,length=0)
arrows(xvals.t,y_lower.t,xvals.t,y_upper.t,code=3,length=0)
points(xvals.c,yvals.c,pch=21,cex=1.5,col="black",bg="white")
axis(side=1,at=c(1:6),tick=T,labels=c("BC","BH","NC","NH","SC","SH"),cex.axis=1.2)

#Predicted potential recruitment

#select age
data<-totdata2[totdata2$Age==5,]
rdata<-data

plot(rdata$ModelledBa,rdata$IngBest,pch=NA,las=1,cex.axis=1.2)
for (cell in unique(rdata$Lat_Lon)){
  celldata<-rdata[rdata$Lat_Lon==cell,]
  for (pft in unique(celldata$Pft)){
    pftdata<-celldata[celldata$Pft==pft,]
    pftdata<-pftdata[order(pftdata$ModelledBa),]
    points(pftdata$ModelledBa,pftdata$IngBest,col=pftdata$col,pch=21,cex=1.1)
    lines(pftdata$ModelledBa,pftdata$IngBest,col=pftdata$col,lwd=2)
  }
}

yvals.c<-c(cell.slope[cell.slope$DemoRate=="IngBest",]$mean.slope)
yvals.t<-c(tot.slope[tot.slope$DemoRate=="IngBest",]$slope50)
y_upper.c<-c(cell.slope[cell.slope$DemoRate=="IngBest",]$max.slope)
y_upper.t<-c(tot.slope[tot.slope$DemoRate=="IngBest",]$slope97.5)
y_lower.c<-c(cell.slope[cell.slope$DemoRate=="IngBest",]$min.slope)
y_lower.t<-c(tot.slope[tot.slope$DemoRate=="IngBest",]$slope2.5)
plot(xvals.t,yvals.t,pch=21,col="black",bg="black",xaxt="n",ylim=c(-30,80),
     xlim=c(0.5,6.5),las=1,cex=1.5,cex.axis=1.2)
abline(0,0,lty=2)
arrows(xvals.c,y_lower.c,xvals.c,y_upper.c,code=3,length=0)
arrows(xvals.t,y_lower.t,xvals.t,y_upper.t,code=3,length=0)
points(xvals.c,yvals.c,pch=21,cex=1.5,col="black",bg="white")
axis(side=1,at=c(1:6),tick=T,labels=c("BC","BH","NC","NH","SC","SH"),cex.axis=1.2)

par(xpd=NA)
legend(-8.8,380,c("BC","BH","NC","NH","SC","SH"),
       col=c("darkblue","skyblue2","darkgreen","lightgreen","darkred","red"),
       pch=c(21,21,21,21,21,21),ncol=6,bty="n",pt.bg="white",pt.cex=1.2,cex=1.4,lty=1,lwd=2)

##############################################################
#FIGURE: Predicted potential growth, mortality, and
#recruitment (incl upper and lower bounds against climate)
#RUN "kriging testing 2" script to generate predictions ("pardata")

pardata$Lat_Lon<-paste(pardata$Lat,pardata$Lon,sep="_")

x11(6.5,6.5)
par(mar=c(0.5,0.5,0.5,0.5),oma=c(4,5.5,1.75,0.5))
layout(matrix(1:6,nrow=3,byrow=T))

y.lim = list(
  G = c(exp(-4.5),exp(1.4)),
  M = c(exp(2.5),exp(6.2)),
  R = c(exp(-2),exp(6.2))
)

x.lim.Temp = list(
  BC = c(2,9),
  BH = c(2,10),
  NC = c(2,14),
  NH = c(2,20),
  SC = c(7,23),
  SH = c(6,22)
)

#This is a guess, what values should be included?
#1 and 99 percentiles result in different values for temperature?
x.lim.Precip = list(
  BC = c(45,118),
  BH = c(44,115),
  NC = c(46,152),
  NH = c(46,143),
  SC = c(68,144),
  SH = c(62,150)
)

#Graph panels: first Temp
for (iDemo in demo) {
      
    pardata.sub <- subset(pardata,Demo==iDemo)
    
    with(pardata.sub,plot(Temp,exp(LogValue),pch=NA,las=1,xaxt="n",
                          ylim=y.lim[[iDemo]],xlab="",ylab="",log="y",cex.axis=1.2))
        
    if(iDemo=="R") axis(side=1,at=c(5,10,15,20),tick=T,labels=T,cex.axis=1.2) else axis(side=1,at=c(5,10,15,20),tick=T,labels=F)
          
    for (iSp in species[-1]) {
       
      pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
      
      col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","skyblue2",ifelse(iSp=="NC","darkgreen",
                                    ifelse(iSp=="NH","lightgreen",ifelse(iSp=="SC","darkred","red")))))
      col.light=rgb(t(col2rgb(col)),alpha=128,maxColorValue=255)
                  
        pd.horz <- data.frame(
        Lat_Lon = unique(pardata.sub.Sp$Lat_Lon),
        LogValueMin = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Lat_Lon,FUN=min),
        LogValueMed = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Lat_Lon,FUN=median),
        LogValueMax = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Lat_Lon,FUN=max),
        Temp = tapply(pardata.sub.Sp$Temp,pardata.sub.Sp$Lat_Lon,FUN=mean),
        Pft = iSp,
        Demo = iDemo  
        )
        
        #add predictions
        pd.horz$PredMin = exp(predict(loess(pd.horz$LogValueMin~pd.horz$Temp)))
        pd.horz$PredMed = exp(predict(loess(pd.horz$LogValueMed~pd.horz$Temp)))
        pd.horz$PredMax = exp(predict(loess(pd.horz$LogValueMax~pd.horz$Temp)))
        pd.horz2<-pd.horz[pd.horz$Temp>x.lim.Temp[[iSp]][1] & pd.horz$Temp<x.lim.Temp[[iSp]][2],]
 
        #add lines
        pd.horz2<-pd.horz2[order(pd.horz2$Temp),]
        lines(pd.horz2$Temp,pd.horz2$PredMin,col=col.light,lwd=2)
        lines(pd.horz2$Temp,pd.horz2$PredMed,col=col,lwd=2)
        lines(pd.horz2$Temp,pd.horz2$PredMax,col=col.light,lwd=2)
        
  }
  
  ################
  #Same for Precip
  with(pardata.sub,plot(Precip,exp(LogValue),pch=NA,ylim=y.lim[[iDemo]],
                        xaxt="n",yaxt="n",xlab="",ylab="",log="y"))
  axis(side=1,at=c(40,60,80,100,120,140),tick=T,labels=F)
  
  if(iDemo=="G") axis(side=2,at=c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5),tick=T,labels=F)
  if(iDemo=="M") axis(side=2,at=c(10,20,50,100,200,500),tick=T,labels=F) else axis(side=2,at=c(0.1,0.5,1,5,10,50,100,500),tick=T,labels=F)
  
  if(iDemo=="R") axis(side=1,at=c(40,60,80,100,120,140),tick=T,labels=T,cex.axis=1.2) else axis(side=1,at=c(5,10,15,20),tick=T,labels=F)
    
  for (iSp in species[-1]) {
    
    pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
    
    col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","skyblue2",ifelse(iSp=="NC","darkgreen",
                                                ifelse(iSp=="NH","lightgreen",ifelse(iSp=="SC","darkred","red")))))
    col.light=rgb(t(col2rgb(col)),alpha=128,maxColorValue=255)
    
    pd.horz <- data.frame(
      Lat_Lon = unique(pardata.sub.Sp$Lat_Lon),
      LogValueMin = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Lat_Lon,FUN=min),
      LogValueMed = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Lat_Lon,FUN=median),
      LogValueMax = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Lat_Lon,FUN=max),
      Precip = tapply(pardata.sub.Sp$Precip,pardata.sub.Sp$Lat_Lon,FUN=mean),
      Pft = iSp,
      Demo = iDemo    
    )
    
    #add predictions
    pd.horz$PredMin = exp(predict(loess(pd.horz$LogValueMin~pd.horz$Precip)))
    pd.horz$PredMed = exp(predict(loess(pd.horz$LogValueMed~pd.horz$Precip)))
    pd.horz$PredMax = exp(predict(loess(pd.horz$LogValueMax~pd.horz$Precip)))
    pd.horz2<-pd.horz[pd.horz$Precip>x.lim.Temp[[iSp]][1] & pd.horz$Precip<x.lim.Precip[[iSp]][2],]
    
    #add lines
    pd.horz2<-pd.horz2[order(pd.horz2$Precip),]
    lines(pd.horz2$Precip,pd.horz2$PredMin,col=col.light,lwd=2)
    lines(pd.horz2$Precip,pd.horz2$PredMed,col=col,lwd=2)
    lines(pd.horz2$Precip,pd.horz2$PredMax,col=col.light,lwd=2)
    
  }
}

par(xpd=NA)
legend(-75,exp(27),c("BC","BH","NC","NH","SC","SH"),
       col=c("darkblue","skyblue2","darkgreen","lightgreen","darkred","red"),ncol=6,bty="n",
       lty=1,lwd=2,cex=1.4)

mtext("Growth", side=2, line=3.5, outer=TRUE, adj=0.85, cex=1)
mtext("Mortality", side=2, line=3.5, outer=TRUE, adj=0.5, cex=1)
mtext("Recruitment", side=2, line=3.5, outer=TRUE, adj=0.1, cex=1)
mtext("Mean annual precipitation (mm)", side=1, line=2.5, 
      outer=TRUE, adj=0.9, cex=1)
mtext(expression(paste("Mean annual temperature (",degree,"C)")), 
      side=1, line=2.75, outer=TRUE, adj=0.1, cex=1)

################################################
#FIGURE: Observed vs. predicted BA per stand age
#Total BA may be best?
#Use same age classes, or just ages, as for size distributions.

cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)
states=c("Minn","Wisc","Iowa","Ill","Indian","Missou","Arkan","Mississ","Kentu","Tenn","Alaba",
         "Flor","Georg","South Car","North Car","Virgini","West Vir","Ohio","Penns","Maryl",
         "Delaw","New Jer","New Yo","Connect","Rhod","Massach","Vermo","New Ham","Main","Michi")


panel.fiaba <- function(time,iPft) {
  
  #x <- fiaba[iPft,time,2,,] #not correct, looks at one plot class only...
  x <- apply(fiaba[iPft,time,,,],c(2,3),FUN=mean) #MCV
  #x <- fiaba[iPft,time,mean(0:2),,] #average of plot classes
  #x[nplots[time,2,,]<5] <- NA  #not correct, one plot class only...
  x[apply(nplots[time,,,],c(2,3),FUN=sum)<5] <- NA #MCV
  x[x>20] <- 20
  image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)
  
}

panel.modelledba <- function(time,iPft) {
  
  #x <- modelledba[time,iPft,mean(0:2),,]
  #x <- modelledba[time,iPft,2,,]
  x <- apply(modelledba[time,iPft,,,],c(2,3),FUN=mean) #MCV
  x[x>20] <- 20
  image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)

}

x11(6.5,4.5)
par(mfrow=c(3,4),mar=c(0,0,0,0),oma=c(0,4,3,0))
panel.fiaba(2,5)
panel.modelledba(2,5)
panel.fiaba(6,5)
panel.modelledba(6,5)
panel.fiaba(2,6)
panel.modelledba(2,6)
panel.fiaba(6,6)
panel.modelledba(6,6)
panel.fiaba(2,7)
panel.modelledba(2,7)
panel.fiaba(6,7)
panel.modelledba(6,7)

mtext("20 years",side=3, line=1.5, outer=TRUE, adj=0.2, cex=1.2)
mtext("60 years",side=3, line=1.5, outer=TRUE, adj=0.77, cex=1.2)
mtext("observed",side=3, line=-0.05, outer=TRUE, adj=0.06, cex=1)
mtext("predicted",side=3, line=-0.05, outer=TRUE, adj=0.35, cex=1)
mtext("observed",side=3, line=-0.05, outer=TRUE, adj=0.63, cex=1)
mtext("predicted",side=3, line=-0.05, outer=TRUE, adj=0.9, cex=1)
mtext("NH",side=1, line=-28, outer=TRUE, adj=-0.07, cex=1.2)
mtext("SC",side=1, line=-17, outer=TRUE, adj=-0.07, cex=1.2)
mtext("SH",side=1, line=-7, outer=TRUE, adj=-0.07, cex=1.2)

############################################################################
############################################################################
#OLD FIGURE: Another option for comparing among vs. within grid cell effects.
#Grid cell mean BA against grid cell mean demographic rate per PFT
#Second series of panels, plot class BA - mean BA, vs. 
#plot class demography - mean demography
#Illegible.

#Calculate means per grid cell and functional type
Cellmean<-aggregate(totdata[,c("ModelledBa","GrowthBest","MortBest","IngBest")],
                    list(totdata$Lon_Lat,totdata$Pft),mean,na.rm=T)
names(Cellmean)<-c("Lon_Lat","Pft","meanModelledBa","meanGrowthBest",
                   "meanMortBest","meanIngBest")

#add means to totdata and add within-grid cell changes
totdata2<-merge(totdata,Cellmean,all.x=T)
totdata2$changeModelledBa<-totdata2$ModelledBa-totdata2$meanModelledBa
totdata2$changeGrowthBest<-totdata2$GrowthBest-totdata2$meanGrowthBest
totdata2$changeMortBest<-totdata2$MortBest-totdata2$meanMortBest
totdata2$changeIngBest<-totdata2$IngBest-totdata2$meanIngBest

x11(6.5,8)
par(mfrow=c(3,2),mar=c(4,4,1,1))
plot(totdata2$meanGrowthBest,totdata2$meanModelledBa,col=totdata2$col,pch=16,cex=0.8)
for (pft in unique(totdata2$Pft)){
  data<-totdata2[totdata2$Pft==pft,]
  abline(lm(data$meanModelledBa~data$meanGrowthBest),lwd=2,col=data$col)
}
plot(totdata2$changeGrowthBest,totdata2$changeModelledBa,col=totdata2$col,pch=16,cex=0.8)
for (pft in unique(totdata2$Pft)){
  data<-totdata2[totdata2$Pft==pft,]
  abline(lm(data$changeModelledBa~data$changeGrowthBest),lwd=2,col=data$col)
}
plot(totdata2$meanMortBest,totdata2$meanModelledBa,col=totdata2$col,pch=16,cex=0.8)
for (pft in unique(totdata2$Pft)){
  data<-totdata2[totdata2$Pft==pft,]
  abline(lm(data$meanModelledBa~data$meanMortBest),lwd=2,col=data$col)
}
plot(totdata2$changeMortBest,totdata2$changeModelledBa,col=totdata2$col,pch=16,cex=0.8)
for (pft in unique(totdata2$Pft)){
  data<-totdata2[totdata2$Pft==pft,]
  abline(lm(data$changeModelledBa~data$changeMortBest),lwd=2,col=data$col)
}
plot(totdata2$meanIngBest,totdata2$meanModelledBa,col=totdata2$col,pch=16,cex=0.8)
for (pft in unique(totdata2$Pft)){
  data<-totdata2[totdata2$Pft==pft,]
  abline(lm(data$meanModelledBa~data$meanIngBest),lwd=2,col=data$col)
}
plot(totdata2$changeIngBest,totdata2$changeModelledBa,col=totdata2$col,pch=16,cex=0.8)
for (pft in unique(totdata2$Pft)){
  data<-totdata2[totdata2$Pft==pft,]
  abline(lm(data$changeModelledBa~data$changeIngBest),lwd=2,col=data$col)
}
