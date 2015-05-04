#Predicted PFT BA against predicted potential growth, mortality, and
#recruitment among and within grid cells
#For now at 50 years.
#Use "totdata" dataframe from graphs script

#Log transform demographic rates
totdata$GrowthBestLog<-log(totdata$GrowthBest)
totdata$MortBestLog<-log(totdata$MortBest)
totdata$IngBestLog<-log(totdata$IngBest)

#Add mean ModelledBa, mean GrowthBest, mean MortBest, mean IngBest
meanBa<-aggregate(totdata[,c("ModelledBa","GrowthBestLog","MortBestLog","IngBestLog")],
                  list(totdata$Lon_Lat,totdata$Pft),mean,na.rm=T)
names(meanBa)<-c("Lon_Lat","Pft","meanModelledBa","meanGrowthBestLog","meanMortBestLog","meanIngBestLog")

totdata2<-merge(totdata,meanBa,all.x=T)

###########################
#Predicted potential growth

cell.results.g<-c()
tot.results.g<-c()

x11(6,12)
par(mfrow=c(6,2))

for (pft in unique(totdata2$Pft)){

#select age and Pft
data<-totdata2[totdata2$Age==5 & totdata2$Pft==pft,]
data$Lon_Lat<-factor(data$Lon_Lat)
gdata<-data[order(data$GrowthBestLog,data$Cell),]

#Estimate overall slope, and slopes per grid cell
g_pred<-function(int,slope_all,slope_cell){
  return(int + slope_all*gdata$meanModelledBa + 
           slope_cell*(gdata$ModelledBa-gdata$meanModelledBa))
}

#Likelihood function
g_ll<-function(int,slope_all,slope_cell,slope_cellmean,slope_cellsd,sigma){
  
  pred<-g_pred(int,slope_all,slope_cell[gdata$Lon_Lat])
  
  #likelihood
  loglike<-sum(dnorm(gdata$GrowthBestLog,pred,sigma,log=T))
  
  #parameter hierarchy
  log_slope_hier<-sum(dnorm(slope_cell,slope_cellmean,slope_cellsd,log=T))
  
  return(loglike + log_slope_hier)
}

#Retrieve MCMC output
fb.pars.g<-list(
  int=c(-10,10,1,0,0,1),
  slope_all=c(-10,10,1,0,0,1),
  slope_cell = c(-10,10,1,0,0,1,158),
  slope_cellmean=c(-10,10,2,0,0,1),
  slope_cellsd=c(1e-6,10,2,1,0,1),
  sigma=c(1e-6,10,1,1,0,1)
)

fb.out.g<-filzbach(400000,400000,g_ll,nrow(gdata),fb.pars.g)

#df.fb.out.g<-as.data.frame(fb.out.g)
#write.table(df.fb.out.g,"fb.out.g.txt",row.names=F,quote=F,sep="\t")

#Converged?
g_llvec<-function(x) g_ll(x[1],x[2],x[3:160],x[161],x[162],x[163])
fb.out.g.ll2<-apply(fb.out.g,1,g_llvec)
plot(fb.out.g.ll2,type="l",main=paste(pft))

#Calculate goodness of fit
fb.pm.g<-colMeans(fb.out.g)
pred<-g_pred(fb.pm.g[1],fb.pm.g[2],(fb.pm.g[3:160])[gdata$Lon_Lat])
plot(pred,gdata$GrowthBestLog,main=paste(pft))
abline(0,1)

#Calculate credible intervals
fb.ci.g<-apply(fb.out.g,2,FUN=quantile,probs=c(0.025,0.5,0.975))
cell.slopes.g<-c(min(fb.ci.g[2,3:159]),mean(fb.ci.g[2,3:159]),max(fb.ci.g[2,3:159]),
                 pft,"GrowthBest")
#cell.slopes.g<-c(quantile(fb.ci.g[2,3:159],probs=c(0.025,0.5,0.975)),pft,"GrowthBest")
cell.results.g<-rbind(cell.results.g,cell.slopes.g)
tot.slope.g<-c(fb.ci.g[,2],pft,"GrowthBest")
tot.results.g<-rbind(tot.results.g,tot.slope.g)

}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

#Predicted potential mortality

cell.results.m<-c()
tot.results.m<-c()

x11(6,12)
par(mfrow=c(6,2))

for (pft in unique(totdata2$Pft)){

#select age and Pft
data<-totdata2[totdata2$Age==5 & totdata2$Pft==pft,]
data$Lon_Lat<-factor(data$Lon_Lat)
mdata<-data[order(data$MortBestLog,data$Cell,data$Pft),]

#Estimate overall slope, and slopes per grid cell
m_pred<-function(int,slope_all,slope_cell){
  return(int + slope_all*mdata$meanModelledBa + 
           slope_cell*(mdata$ModelledBa-mdata$meanModelledBa))
}

#Likelihood function
m_ll<-function(int,slope_all,slope_cell,slope_cellmean,slope_cellsd,sigma){
  
  pred<-m_pred(int,slope_all,slope_cell[mdata$Lon_Lat])
  
  #likelihood
  loglike<-sum(dnorm(mdata$MortBestLog,pred,sigma,log=T))
  
  #parameter hierarchy
  log_slope_hier<-sum(dnorm(slope_cell,slope_cellmean,slope_cellsd,log=T))
  
  return(loglike + log_slope_hier)
}

#Retrieve MCMC output
fb.pars.m<-list(
  int=c(-10,10,1,0,0,1),
  slope_all=c(-10,10,1,0,0,1),
  slope_cell = c(-10,10,1,0,0,1,158),
  slope_cellmean=c(-10,10,2,0,0,1),
  slope_cellsd=c(1e-6,10,2,1,0,1),
  sigma=c(1e-6,10,1,1,0,1)
)

fb.out.m<-filzbach(400000,400000,m_ll,nrow(mdata),fb.pars.m)

#df.fb.out.m<-as.data.frame(fb.out.m)
#write.table(df.fb.out.m,"fb.out.m.txt",row.names=F,quote=F,sep="\t")

#Converged?
m_llvec<-function(x) m_ll(x[1],x[2],x[3:160],x[161],x[162],x[163])
fb.out.m.ll2<-apply(fb.out.m,1,m_llvec)
plot(fb.out.m.ll2,type="l",main=paste(pft))

#Calculate goodness of fit
fb.pm.m<-colMeans(fb.out.m)
pred<-m_pred(fb.pm.m[1],fb.pm.m[2],(fb.pm.m[3:160])[mdata$Lon_Lat])
plot(pred,mdata$MortBestLog,main=paste(pft))
abline(0,1)

#Calculate credible intervals
fb.ci.m<-apply(fb.out.m,2,FUN=quantile,probs=c(0.025,0.5,0.975))
cell.slopes.m<-c(min(fb.ci.m[2,3:159]),mean(fb.ci.m[2,3:159]),max(fb.ci.m[2,3:159]),
                 pft,"MortBest")
#cell.slopes.m<-c(quantile(fb.ci.m[2,3:159],probs=c(0.025,0.5,0.975)),pft,"MortBest")
cell.results.m<-rbind(cell.results.m,cell.slopes.m)
tot.slope.m<-c(fb.ci.m[,2],pft,"MortBest")
tot.results.m<-rbind(tot.results.m,tot.slope.m)

}

#####################################################################################
#####################################################################################

#Predicted potential recruitment

cell.results.r<-c()
tot.results.r<-c()

x11(6,12)
par(mfrow=c(6,2),mar=c())

for (pft in unique(totdata2$Pft)){
  
#select age and Pft
data<-totdata2[totdata2$Age==5 & totdata2$Pft==pft,]
data$Lon_Lat<-factor(data$Lon_Lat)
rdata<-data[order(data$IngBestLog,data$Cell,data$Pft),]

#Estimate overall slope, and slopes per grid cell
r_pred<-function(int,slope_all,slope_cell){
  return(int + slope_all*rdata$meanModelledBa + 
           slope_cell*(rdata$ModelledBa-rdata$meanModelledBa))
}

#Likelihood function
r_ll<-function(int,slope_all,slope_cell,slope_cellmean,slope_cellsd,sigma){
  
  pred<-r_pred(int,slope_all,slope_cell[rdata$Lon_Lat])
  
  #likelihood
  loglike<-sum(dnorm(rdata$IngBestLog,pred,sigma,log=T))
  
  #parameter hierarchy
  log_slope_hier<-sum(dnorm(slope_cell,slope_cellmean,slope_cellsd,log=T))
  
  return(loglike + log_slope_hier)
}

#Retrieve MCMC output
fb.pars.r<-list(
  int=c(-10,10,1,0,0,1),
  slope_all=c(-10,10,1,0,0,1),
  slope_cell = c(-10,10,1,0,0,1,158),
  slope_cellmean=c(-10,10,2,0,0,1),
  slope_cellsd=c(1e-6,10,2,1,0,1),
  sigma=c(1e-6,10,1,1,0,1)
)

fb.out.r<-filzbach(400000,400000,r_ll,nrow(rdata),fb.pars.r)

#df.fb.out.r<-as.data.frame(fb.out.r)
#write.table(df.fb.out.r,"fb.out.r.txt",row.names=F,quote=F,sep="\t")

#Converged?
r_llvec<-function(x) r_ll(x[1],x[2],x[3:160],x[161],x[162],x[163])
fb.out.r.ll2<-apply(fb.out.r,1,r_llvec)
plot(fb.out.r.ll2,type="l",main=paste(pft))

#Calculate goodness of fit
fb.pm.r<-colMeans(fb.out.r)
pred<-r_pred(fb.pm.r[1],fb.pm.r[2],(fb.pm.r[3:160])[rdata$Lon_Lat])
plot(pred,rdata$IngBestLog,main=paste(pft))
abline(0,1)

#Calculate credible intervals
fb.ci.r<-apply(fb.out.r,2,FUN=quantile,probs=c(0.025,0.5,0.975))
cell.slopes.r<-c(min(fb.ci.r[2,3:159]),mean(fb.ci.r[2,3:159]),max(fb.ci.r[2,3:159]),
                 pft,"IngBest")
#cell.slopes.r<-c(quantile(fb.ci.r[2,3:159],probs=c(0.025,0.5,0.975)),pft,"IngBest")
cell.results.r<-rbind(cell.results.r,cell.slopes.r)
tot.slope.r<-c(fb.ci.r[,2],pft,"IngBest")
tot.results.r<-rbind(tot.results.r,tot.slope.r)

}

################################################################################
################################################################################
#Combine estimated slopes in dataframe for graphs

cell.results<-as.data.frame(rbind(cell.results.g,cell.results.m,cell.results.r))
names(cell.results)<-c("min.slope","mean.slope","max.slope","Pft","DemoRate")

write.table(cell.results,"Cell slopes.txt",row.names=F,quote=F,sep="\t")

tot.results<-as.data.frame(rbind(tot.results.g,tot.results.m,tot.results.r))
names(tot.results)<-c("slope2.5","slope50","slope97.5","Pft","DemoRate")

write.table(tot.results,"Overall slopes.txt",row.names=F,quote=F,sep="\t")
