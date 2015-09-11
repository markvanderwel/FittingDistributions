
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
Time<- 0:20 #ncvar_get(nc,"Time")

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
        
        #for (t in 1:16) {    
        for (t in 1:11) {
          
          #if (nplots[t,c,iLon,iLat] > 9) {  
          if (T) {
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
