library(ncdf4)
library(colorRamps)
library(maps)
library(animation)
library(RColorBrewer)

nc <- nc_open("D:\\code\\UF\\io\\FittingOutput.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - NoNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - FirstPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - SecondPassNeighbourhoodPriors.nc", readunlim=F)

nc.orig <- nc_open("D:\\code\\UF\\io\\SimOutWithClimData.nc", readunlim=F)

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

nc_close(nc)

###############################################

celllonlat<- matrix(nrow=length(celltemperature),ncol=2)
for (i in 1:length(celltemperature)) {
  celllonlat[i,1:2] <- which(temperature==celltemperature[i] & precipitation == cellprecipitation[i],arr.ind=T)
}
lonlatcell<-matrix(nrow=length(lon),ncol=length(lat))
for(i in 1:nrow(celllonlat)) {
  lonlatcell[celllonlat[i,1],celllonlat[i,2]]<-i
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

i <- 0

for (iLon in 1:31) {
  for (iLat in 1:27) {
    
    if (!is.na(lonlatcell[iLon,iLat])) {    
      for (t in 1:16) {    
        for (c in 1:3) {
        
        if (nplots[t,c,iLon,iLat] > 9) {         
          for (s in 2:7) {      
            
            i <- i + 1
            
            df.cell[i] <- lonlatcell[iLon,iLat]
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
                      
              
            } #s
          } #nPlots>9
        } #c
      } #t
    } #not water
  } #lat
} #lon
data <- data.frame(Cell=df.cell,Age=df.age,Lon=df.lon,Lat=df.lat,PlotClass=df.plotclass,Pft=df.pft,FiaBa=df.fiaba,OrigBa=df.origba,BestBa=df.bestba,
                   ModelledBa=df.modelledba,GrowthPrior=df.growthprior,MortPrior=df.mortprior,IngPrior=df.ingprior,
                   GrowthBest=df.growthbest,MortBest=df.mortbest,IngBest=df.ingbest,GrowthModelled=df.growthmodelled,
                   MortModelled=df.mortmodelled,IngModelled=df.ingmodelled)



R2 <- function(x,y) {
  ssr=sum((x-y)^2)
  sst=sum((y-mean(y))^2)
  
  r2 = cor(x,y)^2
  
  return (list(R2=1-ssr/sst,r2=r2))
}

#############################################

x11(10,7)
par(mgp=c(0,0,0))
par(mar=c(6,1,1,1))
layout(matrix(1:18,nrow=3))

#plot obs-vs-pred ba
#(by pft)
for (s in 2:7) {
  plot(data$FiaBa[data$Pft==species[s]],data$OrigBa[data$Pft==species[s]],xaxt="n",yaxt="n",xlab="",ylab="",log="")  
  abline(0,1,col="red",lwd=2,lty="solid")
  #abline(coef(lm(data$OrigBa[data$Pft==species[s]]~data$FiaBa[data$Pft==species[s]])),col="dodgerblue",lwd=2)
  r2<-R2(data$FiaBa[data$Pft==species[s]],data$OrigBa[data$Pft==species[s]])
  #title(xlab=paste(round(r2$r2,2),round(r2$R2,2),sep=" ; "))
    
  plot(data$FiaBa[data$Pft==species[s]],data$BestBa[data$Pft==species[s]],xaxt="n",yaxt="n",,xlab="",ylab="",log="")
  abline(0,1,col="red",lwd=2)
  #abline(coef(lm(data$ModelledBa[data$Pft==species[s]]~data$FiaBa[data$Pft==species[s]])),col="dodgerblue",lwd=2)
  r2<-R2(data$FiaBa[data$Pft==species[s]],data$BestBa[data$Pft==species[s]])
  #title(xlab=paste(round(r2$r2,2),round(r2$R2,2),sep=" ; "))
    
  plot(data$FiaBa[data$Pft==species[s]],data$ModelledBa[data$Pft==species[s]],xaxt="n",yaxt="n",,xlab="",ylab="",log="")
  abline(0,1,col="red",lwd=2)
  #abline(coef(lm(data$ModelledBa[data$Pft==species[s]]~data$FiaBa[data$Pft==species[s]])),col="dodgerblue",lwd=2)
  r2<-R2(data$FiaBa[data$Pft==species[s]],data$ModelledBa[data$Pft==species[s]])
  #title(xlab=paste(round(r2$r2,2),round(r2$R2,2),sep=" ; "))
    
}
#############################################
#plot obs-vs-pred ba
#(1 plot)

cols=brewer.pal(6,"Dark2")
cols.light=rgb(t(col2rgb(cols)),alpha=128,maxColorValue=255)


x11(15,5)
layout(matrix(1:3,nrow=1))
par(mgp=c(3,0.7,0))

plot(data$OrigBa,data$FiaBa,xlab="Predicted",ylab="Observed", col=cols.light[match(data$Pft,species)-1],xaxs="i",yaxs="i",xlim=c(-0.01,40.01),ylim=c(-0.01,30.01),xaxp=c(0,40,4),yaxp=c(0,30,3),las=1,cex.axis=2,tcl=0)
title(main="OrigBa")
for (s in 2:7) {
  b <- coef(lm(data$FiaBa[data$Pft==species[s]]~data$OrigBa[data$Pft==species[s]]))
  r <- range(data$OrigBa[data$Pft==species[s]])
  segments(r[1],b[1]+b[2]*r[1],r[2],b[1]+b[2]*r[2], col=cols[s-1],lwd=3)
  segments(r[1],b[1]+b[2]*r[1],r[2],b[1]+b[2]*r[2], col=rgb(0,0,0,0.2),lwd=3) #darken lines a bit
  
  cat(c(species[s],"Orig",unlist(R2(data$FiaBa[data$Pft==species[s]],data$OrigBa[data$Pft==species[s]])),"\n"))
}
abline(0,1,lwd=5,col="grey30",lty="dashed")
box(lwd=2)

cat(c("ALL","Orig",unlist(R2(data$FiaBa,data$OrigBa)),"\n"))

plot(data$BestBa,data$FiaBa,xlab="Predicted",ylab="Observed", col=cols.light[match(data$Pft,species)-1],xaxs="i",yaxs="i",xlim=c(-0.01,30.01),ylim=c(-0.01,30.01),xaxp=c(0,30,3),yaxp=c(0,30,3),las=1,cex.axis=2,tcl=0)
title(main="BestBa")
for (s in 2:7) {
  b <- coef(lm(data$FiaBa[data$Pft==species[s]]~data$BestBa[data$Pft==species[s]]))
  r <- range(data$BestBa[data$Pft==species[s]])
  segments(r[1],b[1]+b[2]*r[1],r[2],b[1]+b[2]*r[2], col=cols[s-1],lwd=3)
  segments(r[1],b[1]+b[2]*r[1],r[2],b[1]+b[2]*r[2], col=rgb(0,0,0,0.2),lwd=3) #darken lines a bit
  
  cat(c(species[s],"Best",unlist(R2(data$FiaBa[data$Pft==species[s]],data$BestBa[data$Pft==species[s]])),"\n"))
}
abline(0,1,lwd=5,col="grey30",lty="dashed")
box(lwd=2)

cat(c("ALL","Best",unlist(R2(data$FiaBa,data$BestBa)),"\n"))

plot(data$ModelledBa,data$FiaBa,xlab="Predicted",ylab="Observed", col=cols.light[match(data$Pft,species)-1],xaxs="i",yaxs="i",xlim=c(-0.01,40.01),ylim=c(-0.01,30.01),xaxp=c(0,40,4),yaxp=c(0,30,3),las=1,cex.axis=2,tcl=0)
title(main="ModelledBa")
for (s in 2:7) {
  b <- coef(lm(data$FiaBa[data$Pft==species[s]]~data$ModelledBa[data$Pft==species[s]]))
  r <- range(data$ModelledBa[data$Pft==species[s]])
  segments(r[1],b[1]+b[2]*r[1],r[2],b[1]+b[2]*r[2], col=cols[s-1],lwd=3)
  segments(r[1],b[1]+b[2]*r[1],r[2],b[1]+b[2]*r[2], col=rgb(0,0,0,0.2),lwd=3) #darken lines a bit
  
  cat(c(species[s],"Modelled",unlist(R2(data$FiaBa[data$Pft==species[s]],data$ModelledBa[data$Pft==species[s]])),"\n"))
}
abline(0,1,lwd=5,col="grey30",lty="dashed")
box(lwd=2)

cat(c("ALL","Modelled",unlist(R2(data$FiaBa,data$ModelledBa)),"\n"))


#############################################
#plot prior-vs-modelled demography

x11(9,3)
layout(matrix(1:3,nrow=1))
par(mgp=c(1.5,0.7,0))
par(mar=c(2,1,2,1))


plot(data$GrowthPrior,data$GrowthModelled,col=cols.light[match(data$Pft,species)-1],xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
abline(0,1,lwd=4,col=rgb(0,0,0,0.7),lty="dashed")
box(lwd=2)

plot(data$MortPrior,data$MortModelled,col=cols.light[match(data$Pft,species)-1],xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
abline(0,1,lwd=4,col=rgb(0,0,0,0.7),lty="dashed")
box(lwd=2)

plot(data$IngPrior,data$IngModelled,col=cols.light[match(data$Pft,species)-1],xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
abline(0,1,lwd=4,col=rgb(0,0,0,0.7),lty="dashed")
box(lwd=2)


#############################################

cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)
states=c("Minn","Wisc","Iowa","Ill","Indian","Missou","Arkan","Mississ","Kentu","Tenn","Alaba","Flor","Georg","South Car","North Car","Virgini","West Vir","Ohio","Penns","Maryl","Delaw","New Jer","New Yo","Connect","Rhod","Massach","Vermo","New Ham","Main","Michi")
ani.options(convert="C:/Progra~1/ImageMagick-6.9.0-Q16/convert.exe")
#ani.options(convert=NULL)
ani.options(ani.width=1000,ani.height=400)
ani.options(autobrowse=F)
ani.options(interval=0.75)

setwd("I:\\analyses\\US fitting CAIN to distributions")
ani.options(outdir = getwd())


make.frame <- function(time,model) {

  layout(matrix(c(1:6,rep(14,6),7:13,rep(14,5)),nrow=4,byrow=T),heights=c(1,0.2,1,0.3))
  par(mar=c(0,0,0,0))
      
  for (iPft in 2:7) {
    
    x <- fiaba[iPft,time,2,,] #middle plot class only
    x[nplots[time,2,,]<5] <- NA
    x[x>20] <- 20
    image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
    map(database="state",regions=states,interior=F,add=T)
    
  }
  for (iPft in 2:7) {
    
    if(model=="orig")
      x <- ba.orig[time*2,iPft,1,,]
    else if (model=="best")
      x <- ba[time,iPft,2,,] #middle plot class only
    else if (model=="modelled")
      x <- modelledba[time,iPft,2,,] #middle plot class only
    
    x[x>20] <- 20
    image(lon,lat,x,col=cols,zlim=c(0,20),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
    map(database="state",regions=states,interior=F,add=T)
  } #for Pft
  
  plot(NULL,NULL,xlim=0:1,ylim=0:1,xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
  text(0.1,0.5, label=paste("Year",Time[time]*10,sep=" "), cex=3, adj=0)
  
}

x11(10,4)
make.frame(8,"orig")

saveGIF({for (i in 1:10) make.frame(i,"orig")}, movie.name="orig_ba.gif")
saveGIF({for (i in 1:10) make.frame(i,"best")}, movie.name="best_ba.gif")
saveGIF({for (i in 1:11) make.frame(i,"modelled")}, movie.name="modelled_ba.gif")

#############################################
#just print 1 animated pft

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

ani.options(ani.width=167,ani.height=368)
saveGIF({for (i in 1:10) make.frame2(i)}, movie.name="orig_ba_single.gif")


#############################################
#map with 1 degree grid overlaying eastern US

map(database="state",regions=states,interior=F,add=F,lwd=2)
x<-apply(nplots,c(3,4),FUN=sum)
for (iLon in 1:dim(x)[1]){
  for (iLat in 1:dim(x)[2]) {
    if (x[iLon,iLat]>0) {
      rect(lon[iLon]-0.5,lat[iLat]-0.5,lon[iLon]+0.5,lat[iLat]+0.5, col="grey98",)
    }
  }
}
map(database="state",regions=states,interior=F,add=T,lwd=2)

#############################################
#map example rate

x11(6,3)
par(mar=c(0,0,0,0))
layout(matrix(1:2,nrow=1))

cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)

x <- matrix(nrow=31,ncol=27)
for (i in 1:dim(cellingpost)[2]) {
  if (cellplotclass[i] == 1) { #only save rates for middle plot class
    x[celllonlat[i,1],celllonlat[i,2]] <- cellingpost[2,i]
    } #if
  } #for cell
x[x<2] <- NA
x[x>250] <- 250

image(lon,lat,(x),col=cols,zlim=c(0,(250)),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
map(database="state",regions=states,interior=F,add=T,lwd=2)


x <- matrix(nrow=31,ncol=27)
for (i in 1:dim(cellingpost)[2]) {
  if (cellplotclass[i] == 1) { #only save rates for middle plot class
    x[celllonlat[i,1],celllonlat[i,2]] <- cellmodelleding[2,i]
  } #if
} #for cell
x[x<2] <- NA
x[x>250] <- 250
image(lon,lat,(x),col=cols,zlim=c(0,(250)),add=F,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1.2)
map(database="state",regions=states,interior=F,add=T,lwd=2)



#############################################
#plot of 1 cell's pft abundance over time
x11()
par(mgp=c(1.5,0.7,0))
cols=c(rep("white",3),rgb(152,185,84,maxColorValue=255),rgb(74,126,187,maxColorValue=255),rgb(190,75,72,maxColorValue=255))
matplot(10*Time[1:9],ba[1:9,2:7,2,12,10],type="l",lwd=8,lty="solid",col=cols,las=1,cex.axis=2,tcl=0,xaxs="i",yaxs="i",ylim=c(0,16),xlab="",ylab="")
matplot(10*Time[1:9],t(fiaba[2:7,1:9,2,12,10]),type="l",lwd=8,lty="dashed",add=T,col=cols)
box(lwd="2")
