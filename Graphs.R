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
                    ifelse(totdata$Pft=="BH","lightblue",
                           ifelse(totdata$Pft=="NC","darkgreen",
                                  ifelse(totdata$Pft=="NH","green",
                                         ifelse(totdata$Pft=="SC","darkred","red")))))

#Add unique code for cell
totdata$Lon_Lat<-paste(totdata$Lon,totdata$Lat,sep="_")

##############################################################
#FIGURE: Goodness of fit: observed vs predicted PFT basal area
#Now one panel with symbols per PFT. Six panels instead?

#One panel: pretty illegible...
plot(totdata$ModelledBa,totdata$FiaBa,xlab="Predicted basal area m2/ha",
     ylab="Observed basal area m2/ha",col=totdata$col,pch=16,cex=0.8)
abline(0,1,lwd=2,lty=2)

#Per functional type instead: 6 panels.
#Or one panel per age class, with colorcode for pft's?

x11(6.5,8)
#pdf("FIG2b.pdf",width=6.5,height=8)
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
##Compare modelled diameter distributions to FIA data (adjust still...)

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
#par(mar=c(3,3,1,1))

x11(6.5,2.5)
#pdf("FIG3.pdf",width=6.5,height=3.5)
par(mfrow=c(1,4),mar=c(4,4,2,1))

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

############################################################################
############################################################################
#FIGURE: Predicted PFT BA against predicted potential growth, mortality, and
#recruitment among and within grid cells

tot.slope<-read.table("Overall slopes.txt",h=T)
cell.slope<-read.table("Cell slopes.txt",h=T)

#select age
data<-totdata[totdata$Age==5,]
data$Lon_Lat<-factor(data$Lon_Lat)
gdata<-data[order(data$GrowthBest,data$Cell,data$Pft),]

x11(width=6.5,height=8)
par(mfrow=c(3,2),mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,0.1,0.1))

plot(gdata$ModelledBa,gdata$GrowthBest,pch=NA,log="")
for (cell in 1:nrow(unique(celllonlat))){
      inclcells <- lonlatclasscell[unique(celllonlat)[cell,1],unique(celllonlat)[cell,2],]
      celldata<-gdata[gdata$Cell %in% inclcells,]
    for (pft in unique(celldata$Pft)){
      pftdata<-celldata[celldata$Pft==pft,]
      lines(pftdata$ModelledBa,pftdata$GrowthBest,col=pftdata$col,type="l")
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
plot(xvals.t,yvals.t,pch=16,xaxt="n",ylim=c(-0.2,0.3),xlim=c(0.5,6.5),cex=1.4)
abline(0,0,lty=2)
points(xvals.c,yvals.c,pch=1,cex=1.4)
arrows(xvals.c,y_lower.c,xvals.c,y_upper.c,code=3,length=0)
arrows(xvals.t,y_lower.t,xvals.t,y_upper.t,code=3,length=0)
axis(side=1,at=c(1:6),tick=T,labels=c("BC","BH","NC","NH","SC","SH"))

#Predicted potential mortality

#select age
data<-totdata2[totdata2$Age==5,]
mdata<-data[order(data$MortBest,data$Cell,data$Pft),]

plot(mdata$ModelledBa,mdata$MortBest,pch=NA,log="")
for (cell in 1:nrow(unique(celllonlat))){
  inclcells <- lonlatclasscell[unique(celllonlat)[cell,1],unique(celllonlat)[cell,2],]
  celldata<-mdata[mdata$Cell %in% inclcells,]
  for (pft in unique(celldata$Pft)){
    pftdata<-celldata[celldata$Pft==pft,]
    lines(pftdata$ModelledBa,pftdata$MortBest,col=pftdata$col,type="l")
  }
}

xvals<-c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2,4.8,5.2,5.8,6.2,6.8,7.2)
yvals<-c(rep(NA,14))
y_upper<-c(rep(NA,14))
y_lower<-c(rep(NA,14))
plot(xvals,yvals,xaxt="n",ylim=c(-0.02,0.2),xlim=c(0.5,6.5))
arrows(xvals,y_lower,xvals,y_upper,code=3,length=0)
axis(side=1,at=c(1:6),tick=T,labels=c("BC","BH","NC","NH","SC","SH"))

#Predicted potential recruitment

#select age
data<-totdata2[totdata2$Age==5,]
rdata<-data[order(data$IngBest,data$Cell,data$Pft),]

plot(rdata$ModelledBa,rdata$IngBest,pch=NA,log="")
for (cell in 1:nrow(unique(celllonlat))){
  inclcells <- lonlatclasscell[unique(celllonlat)[cell,1],unique(celllonlat)[cell,2],]
  celldata<-rdata[rdata$Cell %in% inclcells,]
  for (pft in unique(celldata$Pft)){
    pftdata<-celldata[celldata$Pft==pft,]
    lines(pftdata$ModelledBa,pftdata$IngBest,col=pftdata$col,type="l")
  }
}

xvals<-c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2,4.8,5.2,5.8,6.2,6.8,7.2)
yvals<-c(rep(NA,14))
y_upper<-c(rep(NA,14))
y_lower<-c(rep(NA,14))
plot(xvals,yvals,xaxt="n",ylim=c(-0.02,0.2),xlim=c(0.5,6.5))
arrows(xvals,y_lower,xvals,y_upper,code=3,length=0)
axis(side=1,at=c(1:6),tick=T,labels=c("BC","BH","NC","NH","SC","SH"))

############################################################################
############################################################################
#FIGURE: Another option for comparing among vs. within grid cell effects.
#Grid cell mean BA against grid cell mean demographic rate per PFT
#Second series of panels, plot class BA - mean BA, vs. 
#plot class demography - mean demography 

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
plot(totdata2$changeGrowthBest,totdata2$changeModelledBa,col=totdata2$col,pch=16,cex=0.8)
plot(totdata2$meanMortBest,totdata2$meanModelledBa,col=totdata2$col,pch=16,cex=0.8)
plot(totdata2$changeMortBest,totdata2$changeModelledBa,col=totdata2$col,pch=16,cex=0.8)
plot(totdata2$meanIngBest,totdata2$meanModelledBa,col=totdata2$col,pch=16,cex=0.8)
plot(totdata2$changeIngBest,totdata2$changeModelledBa,col=totdata2$col,pch=16,cex=0.8)

##############################################################
#FIGURE: Predicted potential growth, mortality, and
#recruitment (incl upper and lower bounds against climate)
#RUN "kriging testing 2" script to generate predictions ("pardata")

#cols=brewer.pal(6,"Dark2")
#cols=c("darkblue","lightblue","darkgreen","green","darkred","red")

cols = list(
  BC = "darkblue",
  BH = "lightblue",
  NC = "darkgreen",
  NH = "green",
  SC = "darkred",
  SH = "red"
)

cols.light=rgb(t(col2rgb(cols)),alpha=128,maxColorValue=255)

x11(6,6)
par(mar=c(2,2,0,0))
layout(matrix(1:6,nrow=3,byrow=T))

y.lim = list(
  G = c(exp(-4),exp(1.4)),
  M = c(exp(2),exp(6.2)),
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

col.range.fade <- function(pft,vals){
  orig.col = cols[pft]
  fade.col = rgb(t(col2rgb(orig.col)),alpha=20,maxColorValue=255)
  #fade.col="white"
  return (ifelse(vals>x.lim[[pft]][1] & vals<x.lim[[pft]][2],orig.col,fade.col))
}

#Graph panels: first Temp
for (iDemo in demo) {
      
    pardata.sub <- subset(pardata,Demo==iDemo)
    
    with(pardata.sub,plot(Temp,exp(LogValue),main=paste("\n",iDemo,sep=" "),pch=NA,ylim=y.lim[[iDemo]],xlab="",ylab="",log="y"))
    #with(pardata.sub,plot(Temp,(Prior),col=col.range.fade(which(iSp==species[-1]),Temp),main=paste("\n",iSp,iDemo,sep=" "),ylim=y.lim[[iDemo]],xlab="",ylab="",log="y"))
    #with(pardata.sub,points(Temp,exp(LogValue),col=rgb(0,0,0,0.3),main=paste("Post",iDemo,sep=" ")))  
          
    for (iSp in species[-1]) {
       
      pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
      
      col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","lightblue",ifelse(iSp=="NC","darkgreen",
                                    ifelse(iSp=="NH","green",ifelse(iSp=="SC","darkred","red")))))
      
      indlatlon <- unique(pardata.sub.Sp[,3:4])
      pardata.sub.Sp$Index <- mapply(FUN=function(x1,x2) 
      {which(x1==indlatlon[,1] & x2==indlatlon[,2])},
      pardata.sub.Sp$Lat,
      pardata.sub.Sp$Lon)
    
      pd.horz <- data.frame(
        Index = unique(pardata.sub.Sp$Index),
        LogValueMin = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Index,FUN=min),
        LogValueMed = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Index,FUN=median),
        LogValueMax = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Index,FUN=max),
        Temp = tapply(pardata.sub.Sp$Temp,pardata.sub.Sp$Index,FUN=mean),
        Precip = tapply(pardata.sub.Sp$Precip,pardata.sub.Sp$Index,FUN=mean),
        Pft = iSp
      )
    
      #class.cols = c("red","darkgreen","blue")
      #correct for pft range (min(x.lim.Temp[[iSp]]),max(x.lim.Temp[[iSp]]))
    
      for (iClass in 1:3) {
        p <- exp(predict(loess(pd.horz[,iClass+1] ~ pd.horz$Temp)))          
        #lines(sort(pd.horz$Temp),p[order(pd.horz$Temp)],col=class.cols[iClass],lwd=2)
        lines(sort(pd.horz$Temp),p[order(pd.horz$Temp)],col=col,lwd=2)
      
      } 
        
  }
  
  #Same for Precip
  with(pardata.sub,plot(Precip,exp(LogValue),main=paste("\n",iDemo,sep=" "),pch=NA,ylim=y.lim[[iDemo]],xlab="",ylab="",log="y"))
    
  for (iSp in species[-1]) {
    
    pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
    
    col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","lightblue",ifelse(iSp=="NC","darkgreen",
                 ifelse(iSp=="NH","green",ifelse(iSp=="SC","darkred","red")))))
        
    indlatlon <- unique(pardata.sub.Sp[,3:4])
    pardata.sub.Sp$Index <- mapply(FUN=function(x1,x2) 
    {which(x1==indlatlon[,1] & x2==indlatlon[,2])},
    pardata.sub.Sp$Lat,
    pardata.sub.Sp$Lon)
    
    pd.horz <- data.frame(
      Index = unique(pardata.sub.Sp$Index),
      LogValueMin = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Index,FUN=min),
      LogValueMed = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Index,FUN=median),
      LogValueMax = tapply(pardata.sub.Sp$LogValue,pardata.sub.Sp$Index,FUN=max),
      Temp = tapply(pardata.sub.Sp$Temp,pardata.sub.Sp$Index,FUN=mean),
      Precip = tapply(pardata.sub.Sp$Precip,pardata.sub.Sp$Index,FUN=mean),
      Pft = iSp,
      col = ifelse(iSp=="BC","darkblue",ifelse(iSp=="BH","lightblue",ifelse(iSp=="NC","darkgreen",
                    ifelse(iSp=="NH","green",ifelse(iSp=="SC","darkred","red")))))
    )
    
    #class.cols = c("red","darkgreen","blue")
    
    for (iClass in 1:3) {
      p <- exp(predict(loess(pd.horz[,iClass+1] ~ pd.horz$Precip)))
      #lines(sort(pd.horz$Precip),p[order(pd.horz$Precip)],col=class.cols[iClass],lwd=2)
      lines(sort(pd.horz$Precip),p[order(pd.horz$Precip)],col=col,lwd=2)
      
    } 
    
  }
}

################################################
#FIGURE: Observed vs. predicted BA per stand age
#Total BA may be best?
#Use same age classes, or just ages, as for size distributions.

cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)
states=c("Minn","Wisc","Iowa","Ill","Indian","Missou","Arkan","Mississ","Kentu","Tenn","Alaba",
         "Flor","Georg","South Car","North Car","Virgini","West Vir","Ohio","Penns","Maryl",
         "Delaw","New Jer","New Yo","Connect","Rhod","Massach","Vermo","New Ham","Main","Michi")

#Calculate mean observed and predicted BA per grid cell and stand age
sumBA<-aggregate(totdata[,c("FiaBa","ModelledBa")], 
                 list(totdata$Cell,totdata$Age,totdata$Lon,totdata$Lat,totdata$PlotClass),mean,na.rm=T)
names(sumBA)<-c("Cell","Age","Lon","Lat","PlotClass","FiaBa","ModelledBa")

meanBA<-aggregate(totdata[,c("FiaBa","ModelledBa")], 
                 list(totdata$Cell,totdata$Age,totdata$Lon,totdata$Lat),mean,na.rm=T)
names(meanBA)<-c("Cell","Age","Lon","Lat","meanFiaBa","meanModelledBa")

##DOES NOT WORK YET

par(mfrow=c(8,2),mar=c(0,0,0,0),oma=c(5,6,0.1,0.1))

for(time in c(5)){
  data <- meanBA[meanBA$Age==time,c("Lon","Lat","meanFiaBa")]
  names(data)<-c("x","y","z")
  data$z[data$z>20]<-20
  
  m<-list(x=data$x,y=data$y,z=data$z)
  
  #problem: also y needs to be in ascending order...
  image(m$x,m$y,m$z,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",
        xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)
  
}

#####################################################################################
#####################################################################################
#old for BA maps: ADJUSTED
make.frame2 <- function(time) {
  
  layout(matrix(c(1,2,3),nrow=3,byrow=T),heights=c(1,1,0.3))
  par(mar=c(0,0,0,0))
  
  iPft <- 5
  
  x <- fiaba[iPft,time,2,,] #middle plot class only
  x[nplots[time,2,,]<10] <- NA
  x[x>20] <- 20
  image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)
  
  x <- ba.orig[time*2,iPft,1,,]
  x[x>20] <- 20
  image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)
  
  
  plot(NULL,NULL,xlim=0:1,ylim=0:1,xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
  text(0.25,0.5, label=paste("Year",Time[time]*10,sep=" "), cex=2.5, adj=0)
  
}

x11(1.6667,4)
make.frame2(8)

#####

#df.fiaba[i] <- fiaba[s,t,c,iLon,iLat]
#df.modelledba[i] <- modelledba[t,s,c,iLon,iLat]

x <- fiaba[,5,,,]
x[nplots[5,,,]<10] <- NA
x2 <- apply(x,2,FUN=sum,na.rm=T)
x3 <- apply(x,3,FUN=mean,na.rm=T)
x3[x3>20] <- 20
image(lon,lat,x3,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
map(database="state",regions=states,interior=F,add=T)



