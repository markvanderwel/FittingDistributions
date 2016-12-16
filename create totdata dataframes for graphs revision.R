####################################
##CREATE TOTDATA DATAFRAMES AND SAVE

###################
##Load R data files

#Graphs for all 10 sets of mcmc output
#(1) Figure goodness of fit (10x), and r2's
#Calculate r2 also based on priormedian (per PFT, priormedian is last sample ~403, see sample.no)

#Change working directory
setwd("F:/RData files")

##SET1
load("Env_graphs_set6.RData")

#Samples to select (last 150)
#set1-5: 251:400
#set6-10: 451:600

#For set 6: include separate dataframe based on prior median

#Results set1: select iterations
ba2 <- apply(ba[451:600,,,,,],2:6,FUN=median)
modelledba<-ba2

cellingpost2 <- apply(cellingpost[451:600,,],2:3,FUN=median)
cellmortpost2 <- apply(cellmortpost[451:600,,],2:3,FUN=median)
cellgrowthpost2 <- apply(cellgrowthpost[451:600,,],2:3,FUN=median)

cellmodelleding <- cellingpost2
cellmodelledmort <- cellmortpost2
cellmodelledgrowth <- cellgrowthpost2

###############################################
###############################################
##Make dataframe

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
#df.sample <- integer(0)                         #add sample
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
              
              #for (sample in 1001:1004) {                       #add sample selection
                
                i <- i + 1
                
                df.cell[i] <- lonlatclasscell[iLon,iLat,c]
                df.age[i] <- Time[t]
                df.lon[i] <- lon[iLon]
                df.lat[i] <- lat[iLat]
                df.plotclass[i] <- plotclass[c]
                df.pft[i] <- species[s]
                #df.sample[i] <- sample.no[sample]                  #add sample
                df.fiaba[i] <- fiaba[s,t,c,iLon,iLat]
                df.origba[i] <- ba.orig[t*2,s,1,iLon,iLat]
                df.bestba[i] <- ba2[t,s,c,iLon,iLat]
                df.modelledba[i] <- modelledba[t,s,c,iLon,iLat]
                df.growthprior[i] <- cellgrowthprior[s,df.cell[i]]
                df.mortprior[i] <- cellmortprior[s,df.cell[i]]
                df.ingprior[i] <- cellingprior[s,df.cell[i]]
                df.growthbest[i] <- cellgrowthpost2[s,df.cell[i]]
                df.mortbest[i] <- cellmortpost2[s,df.cell[i]]
                df.ingbest[i] <- cellingpost2[s,df.cell[i]]
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
#} #sample
totdata <- data.frame(Cell=df.cell,Age=df.age,Lon=df.lon,Lat=df.lat,PlotClass=df.plotclass,
                      Pft=df.pft,FiaBa=df.fiaba,
                      #sample=df.sample, #exclude
                      OrigBa=df.origba,BestBa=df.bestba,
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

#Save totdata
#setwd("E:/Documents/UofR Forest dynamics/Git repository DMAR/Model output inverse CAIN MS")
setwd("F:/")

#Results set1
write.csv(totdata,"Totdata set6.csv",row.names=F)