####
library(filzbach)
library(ncdf4)
library(colorRamps)
library(maps)
library(animation)
library(RColorBrewer)

#load color ramp
cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)


# #for looking at iterations 
# fns <- c("FittingOutput - NoNeighbourhoodPriors.nc",
#   "FittingOutput - FirstPassNeighbourhoodPriors.nc",
#   "FittingOutput - SecondPassNeighbourhoodPriors.nc",
#   "FittingOutput - ThirdPassNeighbourhoodPriors.nc",
#   "FittingOutput - FourthPass,TighterPriors.nc",
#   "FittingOutput - FifthPassNeighbourhoodPriors.nc",
#   "FittingOutput - SixthPassNeighbourhoodPriors.nc",
#   "FittingOutput - SeventhPassNeighbourhoodPriors.nc"
#   #"FittingOutput - EighthPassNeighbourhoodPriors.nc",
#   #"FittingOutput.nc"
# )
# p.sd.out<-matrix(nrow=10,ncol=54)
# for (ip in 1:7) {
#   nc <- nc_open(paste("D:\\code\\UF\\io\\",fns[ip],sep=""), readunlim=F)

  

#read in NCDF data file variables

nc <- nc_open("D:\\code\\UF\\io\\FittingOutput.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - NoNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - FirstPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - SecondPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - ThirdPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - FourthPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - FifthPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - SixthPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - SeventhPassNeighbourhoodPriors.nc", readunlim=F)
#nc <- nc_open("D:\\code\\UF\\io\\FittingOutput - EighthPassNeighbourhoodPriors.nc", readunlim=F)

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

#create hash tables between cell and lat/lon indices

celllonlat<- matrix(nrow=length(celltemperature),ncol=2)
for (i in 1:length(celltemperature)) {
  celllonlat[i,1:2] <- which(temperature==celltemperature[i] & precipitation == cellprecipitation[i],arr.ind=T)
}
lonlatcell<-matrix(nrow=length(lon),ncol=length(lat))
for(i in 1:nrow(celllonlat)) {
  lonlatcell[celllonlat[i,1],celllonlat[i,2]]<-i
}
###############################################


#put variables into data frame, one row per demographic rate/pft/plotclass/cell

demo <- c("G","M","R")

df.cell <- numeric(0)
df.pft <- numeric(0)
df.lat <- numeric(0)
df.lon <- numeric(0)
df.plotclass <- numeric(0)
df.demo <- character(0)
df.value <- numeric(0)
df.raw.value <- numeric(0)
df.temp <- numeric(0)
df.precip <- numeric(0)
df.prior <- numeric(0)


i <- 1
for (icell in 1:length(cellplotclass)) {
  
  for (ipft in 2:7) {
    
    rawcellpftdemo <- c(cellgrowthpost[ipft,icell],cellmortpost[ipft,icell],cellingpost[ipft,icell])
    priorcellpftdemo <- c(cellgrowthprior[ipft,icell],cellmortprior[ipft,icell],cellingprior[ipft,icell])
    
    for (idemo in 1:3) {
      
      df.cell[i] <- icell
      df.pft[i] <- species[ipft]
      df.lat[i] <- lat[celllonlat[icell,2]]
      df.lon[i] <- lon[celllonlat[icell,1]]
      df.plotclass[i] <- cellplotclass[icell]
      df.demo[i] <- demo[idemo]
      df.raw.value[i] <- rawcellpftdemo[idemo]
      df.temp[i] <- celltemperature[icell]
      df.precip[i] <- cellprecipitation[icell]
      df.prior[i] <- priorcellpftdemo[idemo]

      
      i <- i + 1
      
    } #idemo
  }#ipft
} #icell

pardata <- data.frame(Cell=df.cell,Pft=df.pft,Lat=df.lat,Lon=df.lon,Temp=df.temp,
                      Precip=df.precip,PlotClass=df.plotclass,Demo=df.demo,
                      RawValue=df.raw.value,Prior=df.prior)


###############################################


#add transformed/scaled variables to data frame

pardata$LogValue <- log(pardata$RawValue)

pardata$sLat <- scale(pardata$Lat)
pardata$sLon <- scale(pardata$Lon)

pardata$sTemp <- scale(pardata$Temp)
pardata$sPrecip <- scale(pardata$Precip)

pardata$LogTrendSurfaceMean <- numeric(nrow(pardata))
pardata$LogTrendSurfaceSD <- numeric(nrow(pardata))
pardata$LogZ <- numeric(nrow(pardata))

#fit trend surface to each demo/species
for (iD in demo) {
  for (iP in species[-1]) {
    
    pd.inds <- which(pardata$Demo==iD & pardata$Pft==iP)
    
    #3rd-order polynomial of lat lon
    #X <- with(pardata[pd.inds,],cbind(1,sLat,sLon,sLat^2,sLat*sLon,sLon^2,sLat^3,sLat^2*sLon,sLat*sLon^2,sLon^3,PlotClass==1,PlotClass==2))
    #3rd-order polynomial of temp precip
    #X <- with(pardata[pd.inds,],cbind(1,sTemp,sPrecip,sTemp^2,sTemp*sPrecip,sPrecip^2,sTemp^3,sTemp^2*sPrecip,sTemp*sPrecip^2,sPrecip^3,PlotClass==1,PlotClass==2))
    #2nd-order polynomial of temp precip
    X <- with(pardata[pd.inds,],cbind(1,sTemp,sPrecip,sTemp^2,sTemp*sPrecip,sPrecip^2))
    y <- pardata$LogValue[pd.inds]
    
    b <- solve(t(X)%*%X)%*%t(X)%*%y

    pardata$LogTrendSurfaceMean[pd.inds] <- X %*% b
    pardata$LogTrendSurfaceSD[pd.inds] <- sd(y - X %*% b)
    
    #standardized residual from log trend surface
    pardata$LogZ[pd.inds] <- (y - X %*% b) / sd(y - X %*% b)
  }
}

#LogZ is response (fitted values from trend surface, normalized to unit variance)
pardata$LogZ <- (pardata$LogValue - pardata$LogTrendSurfaceMean)/pardata$LogTrendSurfaceSD


# #for looking at iterations
# p.sd.out[ip,]<-unique(pardata$LogTrendSurfaceSD)
# }
# matplot(p.sd.out,type="l",lty=1,col="grey80",log="y")


####

#p.out is matrix to hold parameters used in covariance matrices (for kriging)
p.out<-matrix(nrow=18,ncol=2)


#for each demo/pft,
#1. calculate plot class correlations
#2. construct distance matrices
#3. do Filzbach inference of covariance matrix parameters
i<-0
for (iD in demo) {
  for (iP in species[-1]) {
    
    i = i + 1     
      
    pardata.test <- subset(pardata,Demo==iD & Pft==iP)
  
    #quick check of trend surface fits
    #plot(pardata.test$LogTrendSurfaceMean,pardata.test$LogValue,main=paste(iD,iP,sep=" "))
    #print(cor(pardata.test$LogTrendSurfaceMean,pardata.test$LogValue)^2)

    #quick check of plot class correlations
    #pd.wide<-reshape(pardata.test,v.names="LogZ",idvar=c("Lat","Lon"),timevar="PlotClass",direction="wide")
    #wide.cor<-cor(pd.wide[,grep("LogZ.",colnames(pd.wide))],use="p")
    #print(mean(wide.cor[upper.tri(wide.cor)]))
    
    d <- as.matrix(dist(pardata.test[,c("Lat","Lon")]))    
    same.pc <- outer(pardata.test$PlotClass,pardata.test$PlotClass,FUN="==")
    xproducts <- pardata.test$LogZ %*% t(pardata.test$LogZ)
    
    #calculate correlation among different plot classes in same location
    cor.coeff <- sum(((d==0 & !same.pc) * xproducts)/2)/(sum(d==0 & !same.pc)/2 - 1)
    #print(cor.coeff)
    #print("")
    
    d2 <- d[upper.tri(d) & same.pc]
    
    diffsq <- pardata.test$LogZ %*% t(pardata.test$LogZ)
    diffsq2 <- diffsq[upper.tri(diffsq) & same.pc]
        
    dist.grp <- floor(2*d2)/2
    
    d3 <- tapply(d2,dist.grp,FUN=mean)
    diffsq3 <- tapply(diffsq2,dist.grp,FUN=mean)
    n3<- tapply(diffsq2,dist.grp,FUN=length)

    
    #estimate correlation decay with distance using Filzbach
    LL <- function(d.par,sigma) {
      return(sum(n3*dnorm(diffsq3,exp(-d3*d.par),sigma,log=T)))
    }
    fb.pars <- list(
      d.par = c(0.01,10.0,0.5,1,0,1),
      sigma = c(0.001,100.0,0.5,1,0,1))
    
    fb.out <- filzbach(2000,5000,LL,length(d3),fb.pars,thinning=10)
    p.out[i,] <- c(cor.coeff,mean(fb.out[,1]))
    
    
  }
}

###############################################

#function to create covariance matrix 
pred.covar <- function(dist,pc1,pc2,pc.cor,dist.slope) {
  
  pc_rho = ifelse(pc1==pc2,1,pc.cor)
  dist_rho = exp(-dist.slope*dist)
    
  covar = pc_rho * dist_rho
    
  return(covar)
  
}

###############################################

#calculate krigged predictions
i.map<-0

pardata$PriorMeanLog <- numeric(nrow(pardata))
pardata$PriorSDLog <- numeric(nrow(pardata))
  

for (iD in demo) {
  for (iP in species[-1]) {
   
    i.map <- i.map + 1
    
    pd.inds <- which(pardata$Demo==iD & pardata$Pft==iP)
    
    pardata.test <- subset(pardata,Demo==iD & Pft==iP)
        
    #MVN fit
        
    d <- as.matrix(dist(pardata.test[,c("Lat","Lon")]))
    
    n <- length(pardata.test$PlotClass)
    pc <- matrix(rep(pardata.test$PlotClass,n),nrow=n,ncol=n)
    
    cov <- pred.covar(d,pc,t(pc),p.out[i.map,1],p.out[i.map,2])    
    prec <- solve(cov)
    
    mvn.fit <- function(i) {
      
      var.i <- solve(prec[i,i]) #est var
      
      known.vec = pardata.test$LogZ[-i]
      mean.vec = rep(0,nrow(prec))
      
      pred.i <- mean.vec[i] - solve(prec[i,i]) %*% prec[i,-i] %*% (known.vec - mean.vec[-i]) #est mean
      
      return (c(pred.i,var.i))
    }
    
    res <- matrix(nrow=dim(pardata.test)[1],ncol=3)
    
    for (ii in 1:dim(pardata.test)[1]) {
      res[ii,1] <- pardata.test$LogZ[ii]
      res[ii,2:3] <- mvn.fit(ii)  
    }

    z.mvn <- (res[,1]-res[,2])/sqrt(res[,3])
    print(sd(z.mvn))
    
    
    #create set of maps to compare original, trend surface, and krigged values
    #(rows are for each plot class)
    x11(18,12)
    par(mar=c(0,1,2,1))
    layout(matrix(1:9,nrow=3,ncol=3,byrow=T))
        
    
    for (iC in 0:2) {
      mat.out <- matrix(NA,nrow=length(lon),ncol=length(lat))
      mat.out2 <- mat.out
      mat.out3 <- mat.out
      
      for(i in which(pardata.test$PlotClass == iC)) {
        mat.out[celllonlat[pardata.test$Cell[i],1],celllonlat[pardata.test$Cell[i],2]] <- pardata.test$RawValue[i]
        mat.out2[celllonlat[pardata.test$Cell[i],1],celllonlat[pardata.test$Cell[i],2]] <- exp(pardata.test$LogTrendSurfaceMean[i])# + pardata.test$LogTrendSurfaceSD[i]^2/2)
        mat.out3[celllonlat[pardata.test$Cell[i],1],celllonlat[pardata.test$Cell[i],2]] <- exp(pardata.test$LogTrendSurfaceMean[i] + res[i,2])# + res[i,3]/2)

        #mat.out[celllonlat[pardata.test$Cell[i],1],celllonlat[pardata.test$Cell[i],2]] <- pardata.test$LogZ[i] 
        #mat.out2[celllonlat[pardata.test$Cell[i],1],celllonlat[pardata.test$Cell[i],2]] <- 0
        #mat.out3[celllonlat[pardata.test$Cell[i],1],celllonlat[pardata.test$Cell[i],2]] <- res[i,2]
      }
      
      #z.range<-quantile(c(mat.out,mat.out2,mat.out3),probs=c(0.01,0.99),na.rm=T)
      z.range<-quantile(mat.out2,probs=c(0.01,0.99),na.rm=T)
      
      image(lon,lat,mat.out,zlim=z.range,col=cols)
      if (iC==0) title(main=paste(iD,iP,"Original value",sep=", "))
      image(lon,lat,mat.out2,zlim=z.range,col=cols)
      if (iC==0) title(main="Trend surface pred")
      image(lon,lat,mat.out3,zlim=z.range,col=cols)
      if (iC==0) title(main="Krigged pred")
      
      }
    
    #now put appropriate priors into original data frame
    
    #prior from trend surface
    pardata$PriorMeanLog[pd.inds] <- pardata.test$LogTrendSurfaceMean
    pardata$PriorSDLog[pd.inds] <- pardata.test$LogTrendSurfaceSD
    
    #prior from krigging
    #pardata$PriorMeanLog[pd.inds] <- pardata.test$LogTrendSurfaceMean + res[,2]
    #pardata$PriorMeanSD[pd.inds] <- pardata.test$LogTrendSurfaceSD * res[i,3]
  }
}


#save output file
write.csv(pardata,file="PriorsOut.csv")


