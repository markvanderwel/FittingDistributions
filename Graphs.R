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

#Number of plot classes per grid cell
countClass<-aggregate(totdata[,c("PlotClass")], 
                      list(totdata$Cell,totdata$Age),function(x) length(unique(x)))
names(countClass)<-c("Cell","Age","count")
class3<-countClass[countClass$count==3,]
#At 60 years?

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
pdf("FIG2a.pdf",width=3.3,height=3.3)

plot(totdata$ModelledBa,totdata$FiaBa,xlab="Predicted basal area m2/ha",
     ylab="Observed basal area m2/ha",col=totdata$col,pch=16,cex=0.8)
abline(0,1,lwd=2,lty=2)

dev.off()

#Per functional type instead: 6 panels.
pdf("FIG2b.pdf",width=6.5,height=8)
par(mfrow=c(3,2),mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,0.1,0.1))

plot(totdata$ModelledBa[totdata$Pft=="BC"],totdata$FiaBa[totdata$Pft=="BC"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
plot(totdata$ModelledBa[totdata$Pft=="BH"],totdata$FiaBa[totdata$Pft=="BH"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",yaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
axis(side=2,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
plot(totdata$ModelledBa[totdata$Pft=="NC"],totdata$FiaBa[totdata$Pft=="NC"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
plot(totdata$ModelledBa[totdata$Pft=="NH"],totdata$FiaBa[totdata$Pft=="NH"],
     xlim=c(0,30),ylim=c(0,30),las=1,xaxt="n",yaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
axis(side=2,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)
plot(totdata$ModelledBa[totdata$Pft=="SC"],totdata$FiaBa[totdata$Pft=="SC"],
     xlim=c(0,30),ylim=c(0,30),las=1,col="grey60",pch=16,cex=1.4,cex.axis=1.2)
abline(0,1,lwd=2,lty=2)
plot(totdata$ModelledBa[totdata$Pft=="SH"],totdata$FiaBa[totdata$Pft=="SH"],
     xlim=c(0,30),ylim=c(0,30),las=1,yaxt="n",col="grey60",pch=16,cex=1.4,cex.axis=1.2)
axis(side=1,at=seq(0,30,5),labels=F)
abline(0,1,lwd=2,lty=2)

mtext(expression("Observed basal area (m"^2*" ha"^-1*")"),
      side=2, line=2, outer=TRUE, adj=0.5, cex=1)
mtext(expression("Predicted basal area (m"^2*" ha"^-1*")"),
      side=1, line=2.5, outer=TRUE, adj=0.5, cex=1)

dev.off()

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

fian2 <- apply(fian,1,sum)

for (i in 1:21) {
  
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

pdf("FIG3.pdf",width=6.5,height=3.5)
par(mfrow=c(1,4),mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,0.1,0.1))

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
          ylim=c(5e-2,5e3),
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

dev.off()

############################################################################
#FIGURE: Predicted PFT BA against predicted potential growth, mortality, and
#recruitment
#For now at 60 years.

#Predicted potential growth
data60<-totdata[totdata$Age==6,]
data60$Lon_Lat<-factor(data60$Lon_Lat)
gdata60<-data60[order(data60$GrowthBest,data60$Cell,data60$Pft),]

#Estimate overall slope, and slopes per grid cell
#prediction function for species richness
g_pred<-function(int,slope){
  return(int + slope*gdata60$ModelledBa)
}

#Likelihood function
g_ll<-function(ints,slopes,ints_mean,ints_sd,slopes_mean,slopes_sd,sigma){
  
  pred<-g_pred(ints[gdata60$Lon_Lat],slopes[gdata60$Lon_Lat])
  
  #likelihood
  loglike<-sum(dnorm(gdata60$GrowthBest,pred,sigma,log=T))
  
  #parameter hierarchy
  log_int_hier<-sum(dnorm(ints,ints_mean,ints_sd,log=T))
  log_slope_hier<-sum(dnorm(slopes,slopes_mean,slopes_sd,log=T))
  
  return(loglike + log_int_hier + log_slope_hier)
}

#Retrieve MCMC output
fb.pars.g<-list(
  ints=c(-10,10,1,0,0,1,158),
  slopes = c(-10,10,1,0,0,1,158),
  ints_mean=c(-10,10,2,0,0,1),
  ints_sd=c(1e-6,10,2,1,0,1),
  slopes_mean=c(-10,10,2,0,0,1),
  slopes_sd=c(1e-6,10,2,1,0,1),
  sigma=c(1e-6,10,1,1,0,1)
)

fb.out.g<-filzbach(150000,150000,g_ll,nrow(gdata60),fb.pars.g)

#Converged (mmm...)?
g_llvec<-function(x) g_ll(x[1:158],x[159:316],x[317],x[318],x[319],x[320],x[321])
fb.out.g.ll2<-apply(fb.out.g,1,g_llvec)
plot(fb.out.g.ll2,type="l")

#Calculate goodness of fit
fb.pm.g<-colMeans(fb.out.g)
pred<-g_pred((fb.pm.g[1:158])[gdata60$Lon_Lat],(fb.pm.g[159:316])[gdata60$Lon_Lat])
plot(pred,gdata60$GrowthBest)
abline(0,1)

#Calculate credible intervals
fb.ci.g<-apply(fb.out.g,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci.g
cell.slopes.g<-fb.ci.g[2,159:316]
tot.slope.g<-fb.ci.g[,319]

##################
##################
#Make graph panels

plot(gdata60$ModelledBa,gdata60$GrowthBest,pch=NA,log="")
for (cell in 1:nrow(unique(celllonlat))){
      inclcells <- lonlatclasscell[unique(celllonlat)[cell,1],unique(celllonlat)[cell,2],]
      celldata<-gdata60[gdata60$Cell %in% inclcells,]
    for (pft in unique(celldata$Pft)){
      pftdata<-celldata[celldata$Pft==pft,]
      lines(pftdata$ModelledBa,pftdata$GrowthBest,col=pftdata$col,type="l")
   }
}

#Predicted potential mortality
data60<-totdata[totdata$Age==6,]
mdata60<-data60[order(data60$MortBest,data60$Cell,data60$Pft),]

plot(mdata60$ModelledBa,mdata60$MortBest,pch=NA,log="")
for (cell in 1:nrow(unique(celllonlat))){
  inclcells <- lonlatclasscell[unique(celllonlat)[cell,1],unique(celllonlat)[cell,2],]
  celldata<-mdata60[mdata60$Cell %in% inclcells,]
  for (pft in unique(celldata$Pft)){
    pftdata<-celldata[celldata$Pft==pft,]
    lines(pftdata$ModelledBa,pftdata$MortBest,col=pftdata$col,type="l")
  }
}

#Predicted potential recruitment
data60<-totdata[totdata$Age==6,]
rdata60<-data60[order(data60$IngBest,data60$Cell,data60$Pft),]

plot(rdata60$ModelledBa,rdata60$IngBest,pch=NA,log="")
for (cell in 1:nrow(unique(celllonlat))){
  inclcells <- lonlatclasscell[unique(celllonlat)[cell,1],unique(celllonlat)[cell,2],]
  celldata<-rdata60[rdata60$Cell %in% inclcells,]
  for (pft in unique(celldata$Pft)){
    pftdata<-celldata[celldata$Pft==pft,]
    lines(pftdata$ModelledBa,pftdata$IngBest,col=pftdata$col,type="l")
  }
}

##############################################################
#FIGURE: Predicted potential growth, mortality, and
#recruitment (incl upper and lower bounds against climate)
#RUN "kriging testing 2" script to generate predictions ("pardata")

cols=brewer.pal(6,"Dark2")
cols.light=rgb(t(col2rgb(cols)),alpha=128,maxColorValue=255)

x11(6,6)
par(mar=c(2,2,0,0))
layout(matrix(1:6,nrow=3,byrow=T))

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
  return (ifelse(vals>x.lim[[pft]][1] & vals<x.lim[[pft]][2],orig.col,fade.col))
}

#Graph panels: first Temp
for (iDemo in demo) {
      
    pardata.sub <- subset(pardata,Demo==iDemo)
    
    with(pardata.sub,plot(Temp,exp(LogValue),main=paste("\n",iSp,iDemo,sep=" "),pch=NA,ylim=y.lim[[iDemo]],xlab="",ylab="",log="y"))
    #with(pardata.sub,plot(Temp,(Prior),col=col.range.fade(which(iSp==species[-1]),Temp),main=paste("\n",iSp,iDemo,sep=" "),ylim=y.lim[[iDemo]],xlab="",ylab="",log="y"))
    #with(pardata.sub,points(Temp,exp(LogValue),col=rgb(0,0,0,0.3),main=paste("Post",iDemo,sep=" ")))  
          
    for (iSp in species[-1]) {
       
      pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
      
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
        Precip = tapply(pardata.sub.Sp$Precip,pardata.sub.Sp$Index,FUN=mean)
      )
    
      class.cols = c("red","darkgreen","blue")
    
      for (iClass in 1:3) {
        p <- exp(predict(loess(pd.horz[,iClass+1] ~ pd.horz$Temp)))
        lines(sort(pd.horz$Temp),p[order(pd.horz$Temp)],col=class.cols[iClass],lwd=2)
      
      } 
        
  }
  
  #Same for Precip
  with(pardata.sub,plot(Precip,exp(LogValue),main=paste("\n",iSp,iDemo,sep=" "),pch=NA,ylim=y.lim[[iDemo]],xlab="",ylab="",log="y"))
    
  for (iSp in species[-1]) {
    
    pardata.sub.Sp <- subset(pardata.sub, Pft==iSp)
    
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
      Precip = tapply(pardata.sub.Sp$Precip,pardata.sub.Sp$Index,FUN=mean)
    )
    
    class.cols = c("red","darkgreen","blue")
    
    for (iClass in 1:3) {
      p <- exp(predict(loess(pd.horz[,iClass+1] ~ pd.horz$Precip)))
      lines(sort(pd.horz$Precip),p[order(pd.horz$Precip)],col=class.cols[iClass],lwd=2)
      
    } 
    
  }
}

################################################
#FIGURE: Observed vs. predicted BA per stand age
#Total BA (instead of per PFT) may be best?
#Use same age classes, or just ages, as for size distributions.

cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)
states=c("Minn","Wisc","Iowa","Ill","Indian","Missou","Arkan","Mississ","Kentu","Tenn","Alaba",
         "Flor","Georg","South Car","North Car","Virgini","West Vir","Ohio","Penns","Maryl",
         "Delaw","New Jer","New Yo","Connect","Rhod","Massach","Vermo","New Ham","Main","Michi")

#Calculate sum of observed and predicted BA per grid cell and stand age (sums over plot classes... correct?)
sumBA<-aggregate(totdata[,c("FiaBa","ModelledBa")], 
                 list(totdata$Cell,totdata$Age,totdata$Lon,totdata$Lat),sum)
names(sumBA)<-c("Cell","Age","Lon","Lat","sumFiaBa","sumModelledBa")

#Make plots

##DOES NOT WORK YET, DOUBLECHECK DATA (SUM PER PLOT/PLOT CLASS?).
#Original version works fine, but seems to have only one plot class.
#Add total basal area to ncdf file?
pdf("FIG1.pdf",width=3.15,height=12)

par(mfrow=c(8,2),mar=c(0,0,0,0),oma=c(5,6,0.1,0.1))

for(time in c(1,5)){
  data <- sumBA[sumBA$Age==time,c("Lon","Lat","sumFiaBa")]
  names(data)<-c("x","y","z")
  data[data$z>20,]$z<-20
  
  m <- with(data, tapply(z, list(x, y), I))
  image(m,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",
        xlab="",ylab="",bty="n",asp=1.2)
  map(database="state",regions=states,interior=F,add=T)
  #map(database="state",regions=states,interior=F)
  
  #data <- sumBA[sumBA$Age==time,]
  #data[data$sumModelledBa>20,]$sumModelledBa<-20
  #image(data$lon,data$lat,data$sumModelledBa,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",
  #      xlab="",ylab="",bty="n",asp=1.2)
  #map(database="state",regions=states,interior=F,add=T)
  
}

dev.off()


#####################################################################################
#####################################################################################
#old for BA maps
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