library(ncdf4)
library(vegan)

#nc <- nc_open("D:\\code\\UF\\io\\AllEasternFiaDataNoSubplots.nc", readunlim=F)
nc <- nc_open("D:\\code\\UF\\io\\AllEasternFiaDataIncHarvesting.nc", readunlim=F)

spcd <- ncvar_get(nc,"Spcd")
status <- ncvar_get(nc,"RecentStatuscd")
dia <- ncvar_get(nc,"RecentDia")
dia2 <- ncvar_get(nc,"PrevDia")
tpa <- ncvar_get(nc,"RecentTpa_unadj")
tpa2 <- ncvar_get(nc,"PrevTpa_unadj")
treeplotcn <- ncvar_get(nc,"TreePlotCn")
subpcount <- ncvar_get(nc,"SubpCount")
plotcn <- ncvar_get(nc,"PlotCn")
lat <- ncvar_get(nc,"Lat")
lon <- ncvar_get(nc,"Lon")
age <- ncvar_get(nc,"Stdage")
moisture <- ncvar_get(nc,"Physclcd")

nc_close(nc)

treeplot <- match(treeplotcn,plotcn)

#use previous dia/tpa for dead and cut trees
dia[is.na(dia)] <- dia2[is.na(dia)]
tpa[is.na(tpa)] <- tpa2[is.na(tpa)]


treeba <- pi*(dia*2.54/200)^2 * ( 4 * tpa * 2.471 / subpcount[treeplot])

spba <- tapply(treeba,spcd,FUN=sum,na.rm=T)

#list of top 100 sp by BA
sp100 <- as.numeric(names(sort(spba,decreasing=T)[1:100]))

#keep <- spcd %in% sp100 #drop species not in the top 100
keep <- 1:length(treeplotcn)

treeba2 <- treeba[keep]
spcd2 <- spcd[keep]
treeplot2 <- treeplot[keep]

plotspba <- tapply(treeba2,list(treeplot2,spcd2),FUN=sum)
plotspba[is.na(plotspba)] <- 0
plotspba <- plotspba[rowSums(plotspba)>0,] #throw out plots with no BA

cn.ind <- as.numeric(rownames(plotspba))

lat.i <- floor(lat[cn.ind])
lon.i <- floor(lon[cn.ind])

cells <- unique(data.frame(lon.i,lat.i))

plotclasses <- integer(nrow(plotspba))

#cellsp <- matrix(nrow=nrow(cells),ncol=ncol(plotspba))
#colnames(cellsp) <- colnames(plotspba)


for (i in 1:nrow(cells)) {
  cat("\n",i,"\n")
  
  icell <- which(lon.i==cells[i,1] & lat.i==cells[i,2])
  
  if (length(icell) <= 20) { 
    plotclasses[icell] <- 2 #all get 2
    next
  }

  sumspba <- colSums(plotspba[icell,])
  cellplotspba <- plotspba[icell,sumspba>0]  
  
  cellmds <- metaMDS(cellplotspba,k=1)

  #cellage <- (age[cn.ind])[icell]
  #plot(cellage,cellmds$points,ylim=c(-0.001,0.001))
  #print(rownames(cellmds$species)[order(cellmds$species)])

  if (length(icell)>30) { #at least 31 plots
    
    q <- quantile(cellmds$points,probs=c(0.33333,0.66667))
    
    plotclasses[icell] <- 1
    plotclasses[icell[cellmds$points > q[1]]] <- 2
    plotclasses[icell[cellmds$points > q[2]]] <- 3
    
  } #if
  else { #between 21 and 30 plots
    
    q <- median(cellmds$points)
    
    plotclasses[icell] <- 1 #below median gets 1
    plotclasses[icell[cellmds$points > q]] <- 3 #above median gets 3
    
  } #else if
  
} #for (cell)

write.csv(cbind(cn=plotcn[cn.ind],plotclass=plotclasses),file="plotclasses.csv",row.names=F)





