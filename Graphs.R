##GRAPHS

library(ncdf4)
library(colorRamps)
library(maps)
library(RColorBrewer)
library(filzbach)

nc <- nc_open("FittingOutput.nc", readunlim=F)

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
          #for (t in c(2,4,6,8)) {
    
          if (nplots[t,c,iLon,iLat] > 9) {         
          #if (T) {
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

#Add unique code for cell and for region (north vs. south)
totdata$Lat_Lon<-paste(totdata$Lat,totdata$Lon,sep="_")
totdata$region<-ifelse(totdata$Lat<36,"south",ifelse(totdata$Lat<44,"mid","north"))

##############################################################
#FIGURE: Goodness of fit: observed vs predicted PFT basal area
#Now per functional type.
#col=rgb(0,0,0,0.2),pch=16

#x11(6.5,8)
pdf("Figure 3 goodness of fit.pdf",width=6.5,height=8)
par(mfrow=c(3,2),mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,0.1,0.1))

plot(totdata$ModelledBa[totdata$Pft=="BC"],totdata$FiaBa[totdata$Pft=="BC"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",cex=1.8,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(5.75,29.5,"Boreal Conifer",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="BH"],totdata$FiaBa[totdata$Pft=="BH"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",yaxt="n",cex=1.8,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
axis(side=2,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(6.75,29.5,"Boreal Hardwood",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="NC"],totdata$FiaBa[totdata$Pft=="NC"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",cex=1.8,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(11.25,29.5,"Northern Temperate Conifer",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="NH"],totdata$FiaBa[totdata$Pft=="NH"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",yaxt="n",cex=1.8,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
axis(side=2,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(12.25,29.5,"Northern Temperate Hardwood",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="SC"],totdata$FiaBa[totdata$Pft=="SC"],
     xlim=c(0,30),ylim=c(0,30),las=1,cex=1.8,cex.axis=1.2)
abline(0,1,lwd=2,lty=2)
text(11.5,29.5,"Southern Temperate Conifer",cex=1.4)
plot(totdata$ModelledBa[totdata$Pft=="SH"],totdata$FiaBa[totdata$Pft=="SH"],
     xlim=c(0,30),ylim=c(0,30),las=1,yaxt="n",cex=1.8,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
text(12.5,29.5,"Southern Temperate Hardwood",cex=1.4)

mtext(expression("Observed basal area (m"^2*" ha"^-1*")"),
      side=2, line=2, outer=TRUE, adj=0.5, cex=1)
mtext(expression("Predicted basal area (m"^2*" ha"^-1*")"),
      side=1, line=2.5, outer=TRUE, adj=0.5, cex=1)

dev.off()

##############################################################
##FIGURE: Observed vs. predicted size distribution
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

#x11(7.75,3.5)
pdf("Figure 4 size distributions.pdf",width=7.75,height=3.5)
par(mfrow=c(1,4),mar=c(0.5,0.5,0.5,0.5),oma=c(4,5,1.5,0.1))

dbhcomp <- cbind(fiadbh3[,1],modeldbh3[,1])
matplot(dbh-5,dbhcomp,type="b",pch=16,lty=1,col=c("black","grey60"),
          log="y",cex=1.4,ylim=c(3e-2,5e3),xlim=c(0,100),yaxt="n",cex.axis=1.1)
axis(side=2,at=c(0.1,1,10,100,1000),labels=c(0.1,1,10,100,1000),las=1,cex.axis=1.1)

dbhcomp <- cbind(fiadbh3[,2],modeldbh3[,2])
matplot(dbh-5,dbhcomp,type="b",pch=16,lty=1,col=c("black","grey60"),
        log="y",cex=1.4,ylim=c(3e-2,5e3),xlim=c(0,100),yaxt="n",cex.axis=1.1)
axis(side=2,at=c(0.1,1,10,100,1000),labels=F)

dbhcomp <- cbind(fiadbh3[,3],modeldbh3[,3])
matplot(dbh-5,dbhcomp,type="b",pch=16,lty=1,col=c("black","grey60"),
        log="y",cex=1.4,ylim=c(3e-2,5e3),xlim=c(0,100),yaxt="n",cex.axis=1.1)
axis(side=2,at=c(0.1,1,10,100,1000),labels=F)

dbhcomp <- cbind(fiadbh3[,4],modeldbh3[,4])
matplot(dbh-5,dbhcomp,type="b",pch=16,lty=1,col=c("black","grey60"),
        log="y",cex=1.4,ylim=c(3e-2,5e3),xlim=c(0,100),yaxt="n",cex.axis=1.1)
axis(side=2,at=c(0.1,1,10,100,1000),labels=F)

legend("topright",c("Observed","Predicted"),
       lty=1,pch=19,cex=1.3,pt.cex=1.4,col=c("black","grey60"),bty="n")

mtext(expression("Density (stems ha"^-1*")"),side=2, line=3, outer=TRUE, adj=0.5, cex=1)
mtext("DBH midpoint (cm)",side=1, line=2.5, outer=TRUE, adj=0.5, cex=1)

mtext("0-54 y",side=3, line=0.1, outer=TRUE, adj=0.11, cex=0.9)
mtext("55-74 y",side=3, line=0.1, outer=TRUE, adj=0.37, cex=0.9)
mtext(expression("">= 75 *y),side=3, line=-0.15, outer=TRUE, adj=0.63, cex=0.9)
mtext("All",side=3, line=0.1, outer=TRUE, adj=0.88, cex=0.9)

dev.off()

##############################################################
#FIGURE: Predicted potential growth, mortality, and
#recruitment (incl upper and lower bounds against climate)
#RUN "kriging testing 2" script to generate predictions ("pardata")

#x11(6.5,6.5)
pdf("Figure 6 demography-climate.pdf",width=6.5,height=6.5)
par(mar=c(0.5,0.5,0.5,0.5),oma=c(4,5.5,1.75,0.5))
layout(matrix(1:6,nrow=3,byrow=T))

y.lim = list(
  G = c(0.2,3),
  M = c(15,500),
  R = c(2,500)
)

# Limits calculated as the 1st and 99th percentiles of grid cells 
# where PFT had BA >1.0 in FIA plots

x.lim.Temp = list(
  BC = c(2.3,7.25),
  BH = c(2.3,8.4),
  NC = c(2.3,14.2),
  NH = c(2.3,19.1),
  SC = c(10.8,19.6),
  SH = c(5.6,19.3)
)

x.lim.Precip = list(
  BC = c(49,100),
  BH = c(49,97),
  NC = c(49,137),
  NH = c(50,134),
  SC = c(85,136),
  SH = c(62,137)
)

#Graph panels: first Temp
for (iDemo in demo) {
      
    pardata.sub <- subset(pardata,Demo==iDemo)
    
    with(pardata.sub,plot(Temp,exp(LogValue),pch=NA,las=1,xaxt="n",yaxt="n",
                          ylim=y.lim[[iDemo]],xlab="",ylab="",log="y",xlim=c(2,20)))
    
    if(iDemo=="G") axis(side=2,at=c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5),tick=T,labels=T,cex.axis=1.2,las=1)
    if(iDemo=="M") axis(side=2,at=c(10,20,50,100,200,500),tick=T,labels=T,cex.axis=1.2,las=1)
    if(iDemo=="R") axis(side=2,at=c(2,5,20,50,200,500),tick=T,labels=T,cex.axis=1.2,las=1)

    if(iDemo=="R") axis(side=1,at=c(5,10,15,20),tick=T,labels=T,cex.axis=1.2) else axis(side=1,at=c(5,10,15,20),tick=T,labels=F)
          
    for (iSp in species[-1]) {
       
      pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
      pardata.sub.Sp$Lat_Lon<-paste(pardata.sub.Sp$Lat,pardata.sub.Sp$Lon,sep="_")
      
      col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","skyblue2",ifelse(iSp=="NC","darkgreen",
                                    ifelse(iSp=="NH","green3",ifelse(iSp=="SC","red2","pink2")))))
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
        #lines(pd.horz2$Temp,pd.horz2$PredMin,col=col.light,lwd=2)
        lines(pd.horz2$Temp,pd.horz2$PredMin,col=col,lwd=2,lty=2)
        lines(pd.horz2$Temp,pd.horz2$PredMed,col=col,lwd=2)
        #lines(pd.horz2$Temp,pd.horz2$PredMax,col=col.light,lwd=2)
        lines(pd.horz2$Temp,pd.horz2$PredMax,col=col,lwd=2,lty=2)
        
  }
  
  ################
  #Same for Precip
  with(pardata.sub,plot(Precip,exp(LogValue),pch=NA,ylim=y.lim[[iDemo]],
                        xaxt="n",yaxt="n",xlab="",ylab="",log="y", xlim=c(50,140)))
  axis(side=1,at=c(40,60,80,100,120,140),tick=T,labels=F)
  
  if(iDemo=="G") axis(side=2,at=c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5),tick=T,labels=F)
  if(iDemo=="M") axis(side=2,at=c(10,20,50,100,200,500),tick=T,labels=F)
  if(iDemo=="R") axis(side=2,at=c(2,5,20,50,200,500),tick=T,labels=F,cex.axis=1.2)
  
  if(iDemo=="R") axis(side=1,at=c(40,60,80,100,120,140),tick=T,labels=T,cex.axis=1.2)
    
  for (iSp in species[-1]) {
    
    pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
    pardata.sub.Sp$Lat_Lon<-paste(pardata.sub.Sp$Lat,pardata.sub.Sp$Lon,sep="_")
    
    col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","skyblue2",ifelse(iSp=="NC","darkgreen",
                                                ifelse(iSp=="NH","green3",ifelse(iSp=="SC","red2","pink2")))))
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
    pd.horz2<-pd.horz[pd.horz$Precip>x.lim.Precip[[iSp]][1] & pd.horz$Precip<x.lim.Precip[[iSp]][2],]
    
    #add lines
    pd.horz2<-pd.horz2[order(pd.horz2$Precip),]
    #lines(pd.horz2$Precip,pd.horz2$PredMin,col=col.light,lwd=2)
    lines(pd.horz2$Precip,pd.horz2$PredMin,col=col,lwd=2,lty=2)
    lines(pd.horz2$Precip,pd.horz2$PredMed,col=col,lwd=2)
    #lines(pd.horz2$Precip,pd.horz2$PredMax,col=col.light,lwd=2)
    lines(pd.horz2$Precip,pd.horz2$PredMax,col=col,lwd=2,lty=2)
    
  }
}

par(xpd=NA)
legend(-59,exp(20.3),c("BC","BH","NC","NH","SC","SH"),
       col=c("darkblue","skyblue2","darkgreen","green3","red2","pink2"),ncol=6,bty="n",
       lty=1,lwd=2,cex=1.4)

mtext(expression("Growth (cm y"^-1*")"), side=2, line=3, outer=TRUE, adj=0.92, cex=1)
mtext("Longevity (y)", side=2, line=3.25, outer=TRUE, adj=0.52, cex=1)
mtext(expression("Recruitment (ha"^-1*" yr"^-1*")"), side=2, line=3, outer=TRUE, adj=0.04, cex=1)
mtext("Mean annual precipitation (mm)", side=1, line=2.5, 
      outer=TRUE, adj=0.92, cex=1)
mtext(expression(paste("Mean annual temperature (",degree,"C)")), 
      side=1, line=2.75, outer=TRUE, adj=0.1, cex=1)

dev.off()

################################################
#FIGURE: Observed vs. predicted BA per stand age, and BA trajectories

library(maps)
library(mapdata)

##PANEL A
cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)
states=c("Minn","Wisc","Iowa","Ill","Indian","Missou","Arkan","Mississ","Kentu","Tenn","Alaba",
         "Flor","Georg","South Car","North Car","Virgini","West Vir","Ohio","Penns","Maryl",
         "Delaw","New Jer","New Yo","Connect","Rhod","Massach","Vermo","New Ham","Main","Michi")

panel.fiaba <- function(iPft) {
  
  x <- apply(fiaba[iPft,,,,]*nplots,c(3,4),FUN=sum)/apply(nplots,c(3,4),FUN=sum)
  x[apply(nplots,c(3,4),FUN=sum)<5] <- NA
  x[x>20] <- 20
  image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)
  
}

panel.modelledba <- function(iPft) {

  x <- apply(modelledba[,iPft,,,]*nplots,c(3,4),FUN=sum)/apply(nplots,c(3,4),FUN=sum)
  x[x>20] <- 20
  image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)

}

panel.modelledba <- function(iPft) {
  
  x <- apply(modelledba[,iPft,,,]*nplots,c(3,4),FUN=sum)/apply(nplots,c(3,4),FUN=sum)
  x[x>20] <- 20
  image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)
  
}

##Small map icons
#x11(5,3)
pdf("Small map to format.pdf",width=5,height=3)
par(mar=c(0,0,0,0))
plot(1,1,pch=NA,xlim=c(-120,-62),ylim=c(25,55))
map(database="state",regions=states,interior=F,add=T)
rect(-98,20,-65,36)
rect(-98,44,-65,36)
rect(-98,44,-65,64)
dev.off()

#Function for plotting legend
color.bar <- function(cols, min, max, nticks=5, ticks=seq(min, max, len=nticks)) {
  scale = (length(cols)-1)/(max-min)
  
  plot(c(min,max), c(0,1), type="n", bty="n", xaxt="n", xlab=expression("Basal area (m"^2*" ha"^-1*")"), 
       yaxt="n", ylab="",ylim=c(0,1))
  axis(1, ticks, las=1)
  for (i in 1:(length(cols)-1)) {
    x = (i-1)/scale + min
    rect(x,0,x+1/scale,0.75, col=cols[i], border=NA)
  }
}

#x11(6.5,7.5)
pdf("Figure 2 obs vs pred BA maps.pdf",width=6.5,height=7.5)
m<-rbind(c(1,2,2,3,3,4,4,5),c(6,7,7,8,8,9,9,10),c(11,12,12,13,13,14,14,15),
         c(16,16,16,17,17,18,18,18),c(19,19,20,20,20,20,21,21))
layout(m,widths=c(rep(c(1.4,0.6,0.8,0.2,0.2,0.8,0.6,1.4),5)),
       heights=c(1.5,1.5,1.5,0.5,2.5))
par(mar=c(0,0,0,0),oma=c(4,4,1.5,0))

panel.fiaba(2)
panel.modelledba(2)
plot(1,type="n",axes=F,xlab="",ylab="")
panel.fiaba(3)
panel.modelledba(3)
panel.fiaba(4)
panel.modelledba(4)
plot(1,type="n",axes=F,xlab="",ylab="")
panel.fiaba(5)
panel.modelledba(5)
panel.fiaba(6)
panel.modelledba(6)
plot(1,type="n",axes=F,xlab="",ylab="")
panel.fiaba(7)
panel.modelledba(7)

##Add legend
par(mar=c(2,0.5,0,0.5))
color.bar(colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000), 0, 20)
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")

##PANEL B
##Based on totdata: calculate weighted average per plot age.

#create data by region
regions <- ifelse(lat>44,"North",ifelse(lat<36,"South","Mid"))
ba[ba==-999] <- NA

#create arrays to hold results
pred.traj <- array(dim=c(3,dim(ba)[1],length(species[-1])))
dimnames(pred.traj) <- list(c("North","Mid","South"),NULL,species[-1])
obs.traj <- pred.traj

#fill arrays
for (iRegion in c("North","Mid","South")) {
  for (iTime in 1:dim(pred.traj)[2]) {
    for (iPft in species[-1]) {
      pred.traj[iRegion,iTime,iPft] <- weighted.mean(
        ba[iTime,species==iPft,,,regions==iRegion],
        nplots[iTime,,,regions==iRegion],
        na.rm=T)
      obs.traj[iRegion,iTime,iPft] <- weighted.mean(
        fiaba[species==iPft,iTime,,,regions==iRegion],
        nplots[iTime,,,regions==iRegion],
        na.rm=T)
    }
  }
}
#replace NaNs with NA, and set age=0 to 0
pred.traj[is.nan(pred.traj)] <- NA
obs.traj[is.nan(obs.traj)] <- NA
pred.traj[,1,] <- 0
obs.traj[,1,] <- 0

#Calculate means per region, age and functional type
avgdata<-aggregate(totdata[,c("ModelledBa","FiaBa")],
                    list(totdata$region,totdata$Age,totdata$Pft),mean,na.rm=T)
names(avgdata)<-c("region","Age","Pft","meanModelledBa","meanFiaBa")

avgdata2<-avgdata[avgdata$Age<11,]
n.data<-avgdata2[avgdata2$region=="north",]
m.data<-avgdata2[avgdata2$region=="mid",]
s.data<-avgdata2[avgdata2$region=="south",]

##Add in panels with BA-age trajectories
par(mar=c(0.5,0.5,3.5,0.5))

spp<-c("BC","BH","NC","NH")

plot(n.data$Age,n.data$meanModelledBa,pch=NA,xlim=c(0,10),ylim=c(0,20),xlab="Age (y)",
     ylab=expression("Basal area (m"^2*" ha"^-1*")"),las=1,xaxt="n",cex.axis=1.1)
axis(side=1,at=seq(0,10,2),tick=T,labels=c(0,20,40,60,80,100),cex.axis=1.1)
for (iSp in spp){
  col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","skyblue2",ifelse(iSp=="NC","darkgreen",
             ifelse(iSp=="NH","green3","pink2"))))
  
  data<-n.data[n.data$Pft==iSp,]
  data2<-rbind(data,c("north",0,iSp,0,0))
  data2$Age<-as.numeric(data2$Age)
  data2<-data2[order(data2$Pft,data2$Age),]

  #lines(data2$Age,data2$meanFiaBa,col=col,lwd=2)
  #lines(data2$Age,data2$meanModelledBa,col=col,lty=5,lwd=2)
  lines(0:10,obs.traj["North",1:11,match(iSp,species[-1])],col=col,lty=5,lwd=2)
  lines(0:10,pred.traj["North",1:11,match(iSp,species[-1])],col=col,lwd=2)
  
}

spp1<-c("NC","NH","SH")

plot(m.data$Age,m.data$meanModelledBa,pch=NA,xlim=c(0,10),ylim=c(0,20),xlab="Age (y)",
     ylab=expression("Basal area (m"^2*" ha"^-1*")"),las=1,yaxt="n",xaxt="n")
axis(side=1,at=seq(0,10,2),tick=T,labels=c(0,20,40,60,80,100),cex.axis=1.1)
axis(side=2,at=seq(0,25,5),tick=T,labels=F)
for (iSp in spp1){
  col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","skyblue2",ifelse(iSp=="NC","darkgreen",
                                                                       ifelse(iSp=="NH","green3","pink2"))))
  
  data<-m.data[m.data$Pft==iSp,]
  data2<-rbind(data,c("mid",0,iSp,0,0))
  data2$Age<-as.numeric(data2$Age)
  data2<-data2[order(data2$Pft,data2$Age),]
  
  #lines(data2$Age,data2$meanFiaBa,col=col,lwd=2)
  #lines(data2$Age,data2$meanModelledBa,col=col,lty=5,lwd=2)
  lines(0:10,obs.traj["Mid",1:11,match(iSp,species[-1])],col=col,lty=5,lwd=2)
  lines(0:10,pred.traj["Mid",1:11,match(iSp,species[-1])],col=col,lwd=2)
  
}

spp2<-c("NH","SC","SH")

plot(s.data$Age,s.data$meanModelledBa,pch=NA,xlim=c(0,10),ylim=c(0,20),xlab="Age (y)",
     ylab=expression("Basal area (m"^2*" ha"^-1*")"),las=1,yaxt="n",xaxt="n")
axis(side=1,at=seq(0,10,2),tick=T,labels=c(0,20,40,60,80,100),cex.axis=1.1)
axis(side=2,at=seq(0,25,5),tick=T,labels=F)
for (iSp in spp2){
  col = ifelse(iSp=="SC","red2",ifelse(iSp=="NH","green3",ifelse(iSp=="NC","darkgreen","pink2")))
  
  data<-s.data[s.data$Pft==iSp,]
  data2<-rbind(data,c("south",0,iSp,0,0))
  data2$Age<-as.numeric(data2$Age)
  data2<-data2[order(data2$Pft,data2$Age),]
    
  #lines(data2$Age,data2$meanFiaBa,col=col,lwd=2)
  #lines(data2$Age,data2$meanModelledBa,col=col,lty=5,lwd=2)
  lines(0:10,obs.traj["South",1:11,match(iSp,species[-1])],col=col,lty=5,lwd=2)
  lines(0:10,pred.traj["South",1:11,match(iSp,species[-1])],col=col,lwd=2)
  
}

par(xpd=NA)
legend(-24.5,26,c("BC (obs)","BC (pred)","BH (obs)","BH (pred)","NC (obs)","NC (pred)",
                "NH (obs)","NH (pred)","SC (obs)","SC (pred)","SH (obs)","SH (pred)"),                
                col=c("darkblue","darkblue","skyblue2","skyblue2","darkgreen","darkgreen",
                      "green3","green3","red2","red2","pink2","pink2"),ncol=6,bty="n",
                lty=c(rep(c(2,1),6)),lwd=2,cex=1.1)

mtext("observed",side=3, line=-0.05, outer=TRUE, adj=0.06, cex=1)
mtext("predicted",side=3, line=-0.05, outer=TRUE, adj=0.32, cex=1)
mtext("observed",side=3, line=-0.05, outer=TRUE, adj=0.66, cex=1)
mtext("predicted",side=3, line=-0.05, outer=TRUE, adj=0.93, cex=1)
mtext("BC",side=1, line=-47, outer=TRUE, adj=-0.045, cex=1.2)
mtext("NC",side=1, line=-37, outer=TRUE, adj=-0.045, cex=1.2)
mtext("SC",side=1, line=-27, outer=TRUE, adj=-0.045, cex=1.2)
mtext("BH",side=1, line=-47, outer=TRUE, adj=0.515, cex=1.2)
mtext("NH",side=1, line=-37, outer=TRUE, adj=0.515, cex=1.2)
mtext("SH",side=1, line=-27, outer=TRUE, adj=0.515, cex=1.2)
mtext(expression("Basal area (m"^2*" ha"^-1*")"),side=2,line=32.5,adj=0.4,cex=0.9)
mtext("Age (y)",side=1,line=2.5,adj=-0.95,cex=0.9)
mtext("(a)",side=1, line=-52, outer=TRUE, adj=-0.08, cex=1.3)
mtext("(b)",side=1, line=-15.5, outer=TRUE, adj=-0.08, cex=1.3)
mtext(expression("Basal area (m"^2*" ha"^-1*")"),side=1,line=-19.25,adj=-2.5,cex=0.9)

dev.off()
