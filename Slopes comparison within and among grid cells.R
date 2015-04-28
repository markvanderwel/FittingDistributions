#Predicted PFT BA against predicted potential growth, mortality, and
#recruitment among and within grid cells
#For now at 50 years.
#Use "totdata" dataframe from graphs script

#Add mean ModelledBa
meanBa<-aggregate(totdata$ModelledBa,list(totdata$Lon_Lat,totdata$Pft),mean)
names(meanBa)<-c("Lat_Lon","Pft","meanModelledBa")

totdata2<-merge(totdata,meanBa,all.x=T)

###########################
#Predicted potential growth

#select age and Pft
data<-totdata2[totdata2$Age==5 & totdata2$Pft=="NH",]
data$Lon_Lat<-factor(data$Lon_Lat)
#gdata<-data[order(data$GrowthBest,data$Cell,data$Pft),]
gdata<-data[order(data$GrowthBest,data$Cell),]

#Estimate overall slope, and slopes per grid cell
g_pred<-function(int,slope_all,slope_cell){
  return(int + slope_all*gdata$meanModelledBa + 
           slope_cell*(gdata$ModelledBa-gdata$meanModelledBa))
}

#Likelihood function
g_ll<-function(int,slope_all,slope_cell,slope_cellmean,slope_cellsd,sigma){
  
  pred<-g_pred(int,slope_all,slope_cell[gdata$Lon_Lat])
  
  #likelihood
  loglike<-sum(dnorm(gdata$GrowthBest,pred,sigma,log=T))
  
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

fb.out.g<-filzbach(200000,200000,g_ll,nrow(gdata),fb.pars.g)

df.fb.out.g<-as.data.frame(fb.out.g)

write.table(df.fb.out.g,"fb.out.g.txt",row.names=F,quote=F,sep="\t")

#Converged?
g_llvec<-function(x) g_ll(x[1],x[2],x[3:160],x[161],x[162],x[163])
fb.out.g.ll2<-apply(fb.out.g,1,g_llvec)
plot(fb.out.g.ll2,type="l")

#Calculate goodness of fit
fb.pm.g<-colMeans(fb.out.g)
pred<-g_pred(fb.pm.g[1],fb.pm.g[2],(fb.pm.g[3:160])[gdata$Lon_Lat])
plot(pred,gdata$GrowthBest)
abline(0,1)

#Calculate credible intervals
fb.ci.g<-apply(fb.out.g,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci.g
cell.slopes.g<-fb.ci.g[2,3:159]
tot.slope.g<-fb.ci.g[,2]

#Run same model per pft

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

#Predicted potential mortality

#select age
data<-totdata2[totdata2$Age==5,]
data$Lon_Lat<-factor(data$Lon_Lat)
mdata<-data[order(data$MortBest,data$Cell,data$Pft),]

#Estimate overall slope, and slopes per grid cell
m_pred<-function(int,slope_all,slope_cell){
  return(int + slope_all*mdata$meanModelledBa + 
           slope_cell*(mdata$ModelledBa-mdata$meanModelledBa))
}

#Likelihood function
m_ll<-function(int,slope_all,slope_cell,slope_cellmean,slope_cellsd,sigma){
  
  pred<-m_pred(int,slope_all,slope_cell[mdata$Lon_Lat])
  
  #likelihood
  loglike<-sum(dnorm(mdata$MortBest,pred,sigma,log=T))
  
  #parameter hierarchy
  log_slope_hier<-sum(dnorm(slope_cell,slope_cellmean,slope_cellsd,log=T))
  
  return(loglike + log_slope_hier)
}

#Retrieve MCMC output
fb.pars.m<-list(
  int=c(-500,500,1,0,0,1),
  slope_all=c(-500,500,1,0,0,1),
  slope_cell = c(-500,500,1,0,0,1,158),
  slope_cellmean=c(-500,500,2,0,0,1),
  slope_cellsd=c(1e-6,500,2,1,0,1),
  sigma=c(1e-6,500,1,1,0,1)
)

fb.out.m<-filzbach(100000,10000,m_ll,nrow(mdata),fb.pars.m)

df.fb.out.m<-as.data.frame(fb.out.m)

write.table(df.fb.out.m,"fb.out.m.txt",row.names=F,quote=F,sep="\t")

#Converged?
m_llvec<-function(x) m_ll(x[1],x[2],x[3:160],x[161],x[162],x[163])
fb.out.m.ll2<-apply(fb.out.m,1,m_llvec)
plot(fb.out.m.ll2,type="l")

#Calculate goodness of fit
fb.pm.m<-colMeans(fb.out.m)
pred<-m_pred(fb.pm.m[1],fb.pm.m[2],(fb.pm.m[3:160])[mdata$Lon_Lat])
plot(pred,mdata$MortBest)
abline(0,1)

#Calculate credible intervals
fb.ci.m<-apply(fb.out.m,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci.m
cell.slopes.m<-fb.ci.m[2,3:159]
tot.slope.m<-fb.ci.m[,2]

#Run same model per pft

#####################################################################################
#####################################################################################

#Predicted potential recruitment

#select age
data<-totdata2[totdata2$Age==5,]
data$Lon_Lat<-factor(data$Lon_Lat)
rdata<-data[order(data$IngBest,data$Cell,data$Pft),]

#Estimate overall slope, and slopes per grid cell
r_pred<-function(int,slope_all,slope_cell){
  return(int + slope_all*rdata$meanModelledBa + 
           slope_cell*(rdata$ModelledBa-rdata$meanModelledBa))
}

#Likelihood function
r_ll<-function(int,slope_all,slope_cell,slope_cellmean,slope_cellsd,sigma){
  
  pred<-r_pred(int,slope_all,slope_cell[rdata$Lon_Lat])
  
  #likelihood
  loglike<-sum(dnorm(rdata$IngBest,pred,sigma,log=T))
  
  #parameter hierarchy
  log_slope_hier<-sum(dnorm(slope_cell,slope_cellmean,slope_cellsd,log=T))
  
  return(loglike + log_slope_hier)
}

#Retrieve MCMC output
fb.pars.r<-list(
  int=c(-100,100,1,0,0,1),
  slope_all=c(-100,100,1,0,0,1),
  slope_cell = c(-100,100,1,0,0,1,158),
  slope_cellmean=c(-100,100,1,0,0,1),
  slope_cellsd=c(1e-6,100,1,1,0,1),
  sigma=c(1e-6,100,1,1,0,1)
)

fb.out.r<-filzbach(100000,100000,r_ll,nrow(rdata),fb.pars.r)

df.fb.out.r<-as.data.frame(fb.out.r)

write.table(df.fb.out.r,"fb.out.r.txt",row.names=F,quote=F,sep="\t")

#Converged?
r_llvec<-function(x) r_ll(x[1],x[2],x[3:160],x[161],x[162],x[163])
fb.out.r.ll2<-apply(fb.out.r,1,r_llvec)
plot(fb.out.r.ll2,type="l")

#Calculate goodness of fit
fb.pm.r<-colMeans(fb.out.r)
pred<-r_pred(fb.pm.r[1],fb.pm.r[2],(fb.pm.r[3:160])[rdata$Lon_Lat])
plot(pred,rdata$IngBest)
abline(0,1)

#Calculate credible intervals
fb.ci.r<-apply(fb.out.r,2,FUN=quantile,probs=c(0.025,0.5,0.975))
fb.ci.r
cell.slopes.r<-fb.ci.r[2,3:159]
tot.slope.r<-fb.ci.r[,2]

#Run same model per pft
