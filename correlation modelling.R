

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

#can put transformations here if desired
#sqrt.cellgrowthpost <- (cellgrowthpost) 
#sqrt.cellmortpost <- (cellmortpost)
#sqrt.cellingpost <- (cellingpost)

i <- 1

for (icell in 1:length(cellplotclass)) {
  
  for (ipft in 2:7) {
    
    rawcellpftdemo <- c(cellgrowthpost[ipft,icell],cellmortpost[ipft,icell],cellingpost[ipft,icell])
    
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
    
    i <- i + 1
      
    } #idemo
  }#ipft
} #icell

pardata <- data.frame(Cell=df.cell,Pft=df.pft,Lat=df.lat,Lon=df.lon,Temp=df.temp,Precip=df.precip,PlotClass=df.plotclass,Demo=df.demo,RawValue=df.raw.value)

#truncate extreme values
#(doesn't work because need to truncate raw values above)
#(also, see if we can remove outliers with sqrt transformation)
#pardata$Value[pardata$Value>3] <-3

mean.array <- array(dim=c(3,6,3)) #dimensions are demo, pft, plotclass
dimnames(mean.array) <- list(demo,species[2:7],0:2)
sd.array <- mean.array

#save array of means and sd
for (iDemo in demo) {
  for (iPft in species[2:7]) {
    for (iPlotClass in 0:2) {
      
      inds = which((pardata$Demo == iDemo) & (pardata$Pft == iPft) & (pardata$PlotClass == iPlotClass))
      
      mean.array[iDemo,iPft,as.character(iPlotClass)] = mean(pardata[inds,"RawValue"])
      sd.array[iDemo,iPft,as.character(iPlotClass)] = sd(pardata[inds,"RawValue"])
      
    }
  }
}

#add mean and sd to each row
for (i in 1:dim(pardata)[1]) {
  pardata$sd[i] <- sd.array[pardata$Demo[i],pardata$Pft[i],as.character(pardata$PlotClass[i])]
  pardata$means[i] <- mean.array[pardata$Demo[i],pardata$Pft[i],as.character(pardata$PlotClass[i])]
}

#save standardized values in data frame
pardata$Value <- (pardata$RawValue - pardata$means) / pardata$sd

#data.wide.plotclass <- reshape(pardata,v.names="Value",idvar=c("Lat","Lon","Pft","Demo"),timevar="PlotClass",direction="wide")
#data.wide.pft <- reshape(pardata,v.names="Value",idvar=c("Lat","Lon","PlotClass","Demo"),timevar="Pft",direction="wide")
#data.wide.demo <- reshape(pardata,v.names="Value",idvar=c("Lat","Lon","PlotClass","Pft"),timevar="Demo",direction="wide")

#plotclass.cor <- cor(data.wide.plotclass[,7:9],use="pair")
#pft.cor <- cor(data.wide.pft[,7:12],use="pair")
#demo.cor <- cor(data.wide.demo[,7:9],use="pair")

library(ape)

same<-function(x) {
  res = matrix(nrow=length(x),ncol=length(x))
  for (i in 1:length(x))
    res[i,] = (x == x[i])
  return(res)
}

single.row.corrs <- function(i) {
  
  plotclass.nums = match(pardata$PlotClass, 2:0) # order of correlation matrix is 2, 1, 0
  plotclass.corrs = plotclass.cor[plotclass.nums[i],plotclass.nums[inds]]
  
  pft.nums = match(pardata$Pft,species)-1
  pft.corrs = 1 #pft.cor[pft.nums[i],pft.nums]
  
  demo.nums = match(pardata$Demo, demo)
  demo.corrs = 1 #demo.cor[demo.nums[i],demo.nums]
  
  spatial.corrs = exp(predict(glm.cor,newdata=data.frame(d.classes=d[,which(inds==i)])))
    
  #spatial.nums = match(d.floor[,which(inds==i)],c(0,d.classes))
  #spatial.corrs = c(1,moransi)[spatial.nums]
  #spatial.corrs[is.na(spatial.corrs)] = 0
  
  all.corrs = plotclass.corrs * pft.corrs * demo.corrs * spatial.corrs
  names(all.corrs) = NULL
  
  #best.corrs = sort(all.corrs[-i][abs(all.corrs[-i])>0.2],decreasing=T)
  
  return(all.corrs)  
}

all.corrs <- matrix(0, nrow=dim(pardata)[1],ncol=dim(pardata)[1])


for (iDemo in demo) {
  for (iPft in species[2:7]) {
    
    cat(paste("\n", iPft, iDemo ,"\n", sep=" "))
    
    inds <- which((pardata$Demo == iDemo) & (pardata$Pft == iPft))
        
    data.wide.plotclass <- reshape(pardata[inds,],v.names="Value",idvar=c("Lat","Lon","Pft","Demo"),timevar="PlotClass",direction="wide")
    plotclass.cor <- cor(data.wide.plotclass[,grep("Value.",names(data.wide.plotclass))],use="pair")
    print(plotclass.cor)
    
    d <- as.matrix(dist(pardata[inds,c("Lat","Lon")]))
    
    same.plotclass <- same(pardata[inds,"PlotClass"])
            
    #define bins by rounding down grid cell distances
    d.floor <- floor(d)
    #d <- NULL #free up some memory
    max.d <- 15
    d.classes <- 1:max.d
    moransi <- numeric(0)
    
    for (i in 1:length(d.classes)) {
      w <- 1*((d.floor==d.classes[i]) & same.plotclass)
      m = Moran.I(pardata[inds,"Value"],w)
      w <- NULL #free up some memory
      moransi[i] <- m$observed
    }
    
    glm.cor <- glm(moransi ~ d.classes + 0,family=gaussian(link="log"),start=-1)
    #moransi <- exp(predict(glm.cor))
    
    print(moransi)
    
    for (i in inds) {
      all.corrs[i,inds] <- single.row.corrs(i)
    }
    
    d.floor <- NULL #dispose of distance matrix to free up memory
    
  }#for iPft
} #for iDemo


all.covar <- pardata$sd %*% t(pardata$sd) * all.corrs #create covariance matrix
all.covar<-all.corrs

single.row.prior <- function(i,n) {
    
  #use n best-correlated rows 
  known.rows = order(abs(all.corrs[i,]),decreasing=T)[2:(n+1)]
  
  cov.mat = all.covar[c(i,known.rows),c(i,known.rows)]
  prec.mat = solve(cov.mat)
  
  #prec.mat = nearPD(prec.mat)$mat
  
  #mean.vec = pardata$means[c(i,known.rows)]
  #known.vec = pardata$RawValue[known.rows]
  mean.vec = rep(0,length(known.rows)+1)
  known.vec = pardata$Value[known.rows]
  
  ab.inv.bb = cov.mat[1,-1] %*% solve(cov.mat[-1,-1])
  
  #prior.mean = mean.vec[1] + ab.inv.bb %*% (known.vec - mean.vec[-1])
  #prior.sd = sqrt(cov.mat[1,1] - ab.inv.bb %*% cov.mat[-1,1])
  
  prior.mean = mean.vec[1] - solve(prec.mat[1,1]) %*% prec.mat[1,-1] %*% (known.vec - mean.vec[-1])
  prior.sd = sqrt(solve(prec.mat[1,1]))
    
    
  return(c(pardata$RawValue[i],prior.mean,prior.sd))
  
}

priors <- matrix(nrow=dim(pardata)[1],ncol=2)

for (i in 1:dim(pardata)[1]) {
  priors[i,] <- single.row.prior(i,25)[2:3]
}
pardata$PriorMean <- priors[,1]
pardata$PriorSd <- priors[,2]

plot(pardata$PriorMean[pardata$Demo=="G"],pardata$RawValue[pardata$Demo=="G"])
plot(pardata$PriorMean[pardata$Demo=="M"],pardata$RawValue[pardata$Demo=="M"])

plot(pardata$PriorMean[pardata$Demo=="R"],pardata$RawValue[pardata$Demo=="R"])

write.csv(pardata,file="PriorsOut.csv")
