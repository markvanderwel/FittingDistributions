library(fields)
library(filzbach)
library(ncf)

cols<-colorRampPalette(c("blue","dodgerblue","cyan","green","yellow","orange","red"))(1000)

harv <- read.csv("HarvestPlotsList.csv")
harv$PlotAge[harv$PlotAge<0] <- NA

#plot.effects <- read.csv("PlotEffectsOut.csv",row.names=NULL)
#names(plot.effects) <- names(plot.effects)[-1]
#plot.effects[,18] <- NULL

composition <- read.csv("plotclasses.csv")

#matches <- match(harv$PlotId,plot.effects$Cn)
matches <- match(harv$PlotId,composition$cn)
harv$PlotClass <- composition$plotclass[matches]

# harv$PlotEff <- plot.effects$ReGrMean[matches]
# quant <- quantile(harv$PlotEff,probs=c(0.33,0.67),na.rm=T)
# harv$PlotClass <- 3
# harv$PlotClass[harv$PlotEff < quant[2]] <- 2
# harv$PlotClass[harv$PlotEff < quant[1]] <- 1
# harv$PlotClass[is.na(harv$PlotEff)] <- NA

harv$Cells <- factor(trunc(harv$Lon)*100 + trunc(harv$Lat))

harv$PlotClass1 <- ifelse(harv$PlotClass==1,1,0)
harv$PlotClass2 <- ifelse(harv$PlotClass==2,1,0)
harv$PlotClass3 <- ifelse(harv$PlotClass==3,1,0)


i<-0
AnnualHarvMort <- matrix(ncol=3,nrow=length(unique(harv$Cells)))
PlotCount <- matrix(ncol=3,nrow=length(unique(harv$Cells)))
CellLonLat <- matrix(ncol=2,nrow=length(unique(harv$Cells)))
for (c in unique(harv$Cells)) {

  i <- i + 1
  
  CellLonLat[i,1] <- trunc(mean(harv$Lon[harv$Cells==c]))
  CellLonLat[i,2] <- trunc(mean(harv$Lat[harv$Cells==c]))
  
  #for each cell, distributes missing PlotClass values 
  #in the same proportion as the PlotClass distribution 
  #of harvested plots
  
  pharvs <- sweep(table(harv$PlotClass[harv$Cells==c],harv$Harv[harv$Cells==c]),1,table(harv$PlotClass[harv$Cells==c]),FUN="/")
  pharvs <- sweep(pharvs,2,colSums(pharvs),FUN="/")
  pharvs[is.nan(pharvs)] <- mean(pharvs,na.rm=T)
  
  #fill up missing row/columns with zeros
  if (!("1" %in% colnames(pharvs))) pharvs <- cbind(pharvs,"1"=0)
  
  if (!("1" %in% rownames(pharvs))) pharvs <- rbind(pharvs,"1"=0)
  if (!("2" %in% rownames(pharvs))) pharvs <- rbind(pharvs,"2"=0)
  if (!("3" %in% rownames(pharvs))) pharvs <- rbind(pharvs,"3"=0)
  
  harv$PlotClass1[harv$Cells==c & is.na(harv$PlotClass)] <- pharvs["1","1"]
  harv$PlotClass2[harv$Cells==c & is.na(harv$PlotClass)] <- pharvs["2","1"]
  harv$PlotClass3[harv$Cells==c & is.na(harv$PlotClass)] <- pharvs["3","1"]
  
  AnnualHarvMort[i,1] <- sum(harv$PlotClass1[harv$Cells==c] * harv$CutBa[harv$Cells==c]) / (sum(harv$PlotClass1[harv$Cells==c] * harv$Interval[harv$Cells==c] * harv$InitBa[harv$Cells==c]))
  AnnualHarvMort[i,2] <- sum(harv$PlotClass2[harv$Cells==c] * harv$CutBa[harv$Cells==c]) / (sum(harv$PlotClass2[harv$Cells==c] * harv$Interval[harv$Cells==c] * harv$InitBa[harv$Cells==c]))
  AnnualHarvMort[i,3] <- sum(harv$PlotClass3[harv$Cells==c] * harv$CutBa[harv$Cells==c]) / (sum(harv$PlotClass3[harv$Cells==c] * harv$Interval[harv$Cells==c] * harv$InitBa[harv$Cells==c]))
  
  PlotCount[i,1] <- sum(harv$PlotClass1[harv$Cells==c])
  PlotCount[i,2] <- sum(harv$PlotClass2[harv$Cells==c])
  PlotCount[i,3] <- sum(harv$PlotClass3[harv$Cells==c])
  
} #for

colMeans(AnnualHarvMort,na.rm=T)

stack.annharv <- c(AnnualHarvMort[,1],AnnualHarvMort[,2],AnnualHarvMort[,3])
stack.latlon <- rbind(CellLonLat,CellLonLat,CellLonLat)
stack.plots <- c(PlotCount[,1],PlotCount[,2],PlotCount[,3])
  
#remove NAs
stack.annharv2 <- stack.annharv[!is.na(stack.annharv)]
stack.annharv2[stack.annharv2==0] <- 1e-4 #replace zeros with small value
stack.latlon2 <- stack.latlon[!is.na(stack.annharv),]
stack.plots2 <- stack.plots[!is.na(stack.annharv)]

#quick correlogram shows that harvest rates are autocorrelated for up to about 5 grid cells
ncf.cor <- correlog(stack.latlon2[,1],stack.latlon2[,2],stack.annharv2,increment=1,latlon=F,resamp=F)
plot(ncf.cor)


#create indicator distance matrix for cells that are within 3 grid cells of each other
d <- as.matrix(dist(stack.latlon2,method="eucl"))
d <- (d<=3)

#calculate mean harvest rate (p) in neighbourhood of each cell
p = apply(d,1,FUN=function(x) {
  sum(stack.annharv2 * x * stack.plots2)/sum(x * stack.plots2)
})

#compare p to original data
x11()
plot(p,stack.annharv2)
abline(0,1)

#calculate logit of p
logit.p <- log(p/(1-p))


mm.sigma <- sd(logit.p - log(stack.annharv2/(1-stack.annharv2)))


#create likelihood function for random effects applied to each neighbourhood p

loglik <- function(re) {
  
  p.re <- 1/(1+exp(-(logit.p.i + re)))
  
  #work out alpha, beta of distribution based on estimated mean and # of plots
  beta.alpha = p.re * (stack.plots2.i+2)
  beta.beta = (1-p.re) * (stack.plots2.i+2)
  
  #likelihood for data
  ll = dbeta(stack.annharv2.i,beta.alpha,beta.beta,log=T) + sum(dnorm(re,0,mm.sigma,log=T))
  
  return(ll)
}

fb.pars <- list(
  re = c(-5,5,1,0,0,1)
)

re.out <- numeric(length(p))
for (i in 1:length(p)) {
  cat("\n",i,"\n")
  
  logit.p.i <- logit.p[i]
  stack.plots2.i <- stack.plots2[i]
  stack.annharv2.i <- stack.annharv2[i]
  
  fb.out.i <- filzbach(2000,5000,loglik,1,fb.pars,thinning=20)
  
  re.out[i] <- mean(fb.out.i)
  
}

#histogram of random effects
hist(re.out,breaks=50)

#new predicted harvest rates
p.re <- 1/(1+exp(-(logit.p + re.out)))
#compare to data
plot(p.re,stack.annharv2)
abline(0,1)
cor(p.re,stack.annharv2)^2

#compare mean harvest rates
sum(p.re*stack.plots2)/sum(stack.plots2)
sum(stack.annharv2*stack.plots2)/sum(stack.plots2)


#make plots of observed and predicted by plot class
x11()
layout(matrix(1:6,nrow=2))

#arrays to hold data/predictions
dat <- matrix(ncol=length(unique(CellLonLat[,2])),nrow=length(unique(CellLonLat[,1])))
colnames(dat) <- sort(unique(CellLonLat[,2]))
rownames(dat) <- sort(unique(CellLonLat[,1]))
mod <- dat

FittedHarvMort <- matrix(nrow=nrow(PlotCount),ncol=ncol(PlotCount))
RawHarvMort <- FittedHarvMort

#indices marking start of new plot class
start <- which(stack.latlon2[1,1]==stack.latlon2[,1] & stack.latlon2[1,2]==stack.latlon2[,2])
start <- c(start,nrow(stack.latlon2)+1)

#put results into 3-column output   
for (iclass in 1:3) {
  for (icell in start[iclass]:(start[iclass+1]-1)) {
    FittedHarvMort[which(stack.latlon2[icell,1]==CellLonLat[,1] & stack.latlon2[icell,2]==CellLonLat[,2]),iclass] <- p.re[icell]    
    RawHarvMort[which(stack.latlon2[icell,1]==CellLonLat[,1] & stack.latlon2[icell,2]==CellLonLat[,2]),iclass] <- stack.annharv2[icell]    
  }
}

#sort into min,med,max predictions by harvest rate
raw.seq <- t(apply(RawHarvMort,1,FUN=order)) #the sort order for each row

fit.min <- diag(FittedHarvMort[,raw.seq[,1]])
fit.med <- diag(FittedHarvMort[,raw.seq[,2]])
fit.max <- diag(FittedHarvMort[,raw.seq[,3]])
fit.ord <- cbind(fit.min,fit.med,fit.max)

raw.min <- diag(RawHarvMort[,raw.seq[,1]])
raw.med <- diag(RawHarvMort[,raw.seq[,2]])
raw.max <- diag(RawHarvMort[,raw.seq[,3]])                                      
raw.ord <- cbind(raw.min,raw.med,raw.max)

#generate maps
for (iclass in 1:3) {
  for (i in 1:nrow(fit.ord)) {
    mod[as.character(CellLonLat[i,1]),as.character(CellLonLat[i,2])] <- fit.ord[i,iclass]
    dat[as.character(CellLonLat[i,1]),as.character(CellLonLat[i,2])] <- raw.ord[i,iclass]
    
  }
  image((dat),zlim=range((raw.ord),na.rm=T),col=cols)
  image((mod),zlim=range((raw.ord),na.rm=T),col=cols)
}



#plot the trends in (min,median,max) harvest rates among plot classes by latitude
x11()
matplot(CellLonLat[,2],FittedHarvMort,col="white",ylim=c(0,0.04))

points(CellLonLat[,2],predict(loess(fit.min~CellLonLat[,2])))
points(CellLonLat[,2],predict(loess(fit.med~CellLonLat[,2])),col="red")
points(CellLonLat[,2],predict(loess(fit.max~CellLonLat[,2])),col="blue")


dat.out <- cbind(CellLonLat,PlotCount,FittedHarvMort)

dat.out[is.na(dat.out)] <- (-999)

colnames(dat.out) <- c("Lon","Lat","Plots1","Plots2","Plots3","Harv1","Harv2","Harv3")

write.csv(dat.out,"AnnualHarvOut.csv",row.names=F)


