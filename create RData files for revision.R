##SAVE RData files

library(ncdf4)

###############################################
#1. Graphs
#Save RData files for each of the 10 model runs

for(set in 1:10){
  
  cat(paste0(set,"\n"))
  
  setwd(paste("D:/Documents/UofR Forest dynamics/Git repository DMAR/Model output inverse CAIN MS/CAINprogram",set,sep="")) 
  
  nc <- nc_open("FittingOutput.nc", readunlim=F)
  
  setwd("C:/Users/DMAR/Dropbox/Current projects/UofR/Git repository DMAR/FittingDistributions")
  
  nc.orig <- nc_open("SimOutWithClimData.nc", readunlim=F)
  
  ba.orig <- ncvar_get(nc.orig,"Ba")
  
  ba<-ncvar_get(nc,"Ba")
  ba <- ba[-(1:400),,,,,]
  
  #Results set1
  #ba2 <- apply(ba[655:804,,,,,],2:6,FUN=median)
  #modelledba<-ba2 #modelledba<-ncvar_get(nc,"ModelledBa")
  
  fiaba<-ncvar_get(nc,"FiaBa")
  
  nplots <- ncvar_get(nc,"FiaNPlots")
  
  lat<-ncvar_get(nc,"Lat")
  lon<-ncvar_get(nc,"Lon")
  plotclass<-ncvar_get(nc,"PlotClass")
  species<-c("Tot","BC","BH","NC","NH","SC","SH")  #ncvar_get(nc,"Species")
  sample.no<-ncvar_get(nc,"Sample")
  Time<- 0:15 #ncvar_get(nc,"Time")
  
  temperature<-ncvar_get(nc,"Temperature")
  precipitation<-ncvar_get(nc,"Precipitation")
  
  celltemperature<-ncvar_get(nc,"CellTemperature")
  cellprecipitation<-ncvar_get(nc,"CellPrecip")
  
  cellplotclass<-ncvar_get(nc,"CellPlotClass")
  
  cellingprior <- ncvar_get(nc,"CellIngPrior")
  cellmortprior <- ncvar_get(nc,"CellMortPrior")
  cellgrowthprior <- ncvar_get(nc,"CellGrowthPrior")
  
  #Results set1
  cellingpost <- ncvar_get(nc,"CellIngPost")
  #cellingpost2 <- apply(cellingpost[655:804,,],2:3,FUN=median)
  cellingpost <- cellingpost[-(1:400),,]
  
  #Results set1
  cellmortpost <- ncvar_get(nc,"CellMortPost")
  #cellmortpost2 <- apply(cellmortpost[655:804,,],2:3,FUN=median)
  cellmortpost <- cellmortpost[-(1:400),,]
  
  #Results set1
  cellgrowthpost <- ncvar_get(nc,"CellGrowthPost")
  #cellgrowthpost2 <- apply(cellgrowthpost[655:804,,],2:3,FUN=median)
  cellgrowthpost <- cellgrowthpost[-(1:400),,]
  
  #cellmodelleding <- cellingpost2 #cellmodelleding <- ncvar_get(nc,"CellModelledIng")
  #cellmodelledmort <- cellmortpost2 #cellmodelledmort <- ncvar_get(nc,"CellModelledMort")
  #cellmodelledgrowth <- cellgrowthpost2 #cellmodelledgrowth <- ncvar_get(nc,"CellModelledGrowth")
  
  cellingupper95 <- ncvar_get(nc,"CellIngUpper95")
  cellmortupper95 <- ncvar_get(nc,"CellMortUpper95")
  cellgrowthupper95 <- ncvar_get(nc,"CellGrowthUpper95")
  
  cellinglower95 <- ncvar_get(nc,"CellIngLower95")
  cellmortlower95 <- ncvar_get(nc,"CellMortLower95")
  cellgrowthlower95 <- ncvar_get(nc,"CellGrowthLower95")
  
  #For figure size distributions (COMMENT OUT IF POSSIBLE, TAKES LONG)
  fiadbh<-ncvar_get(nc,"FiaDbhDens")
  modeldbh<-ncvar_get(nc,"DbhDens")
  fian<-ncvar_get(nc,"FiaNPlots")
  dbh<-ncvar_get(nc,"Dbh")
  
  nc_close(nc)
  
  #rm(list=ls())
  
  #Save Rdata file
  setwd("D:/Documents/UofR Forest dynamics/Git repository DMAR/RData files")
  save(list = ls(all.names = TRUE), file = paste("Env_graphs_set",set,".RData",sep=""), envir = .GlobalEnv)
  
}