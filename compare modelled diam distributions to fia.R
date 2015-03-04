library(ncdf4)

nc <- nc_open("D:\\code\\UF\\io\\FittingOutput.nc", readunlim=F)

fiadbh<-ncvar_get(nc,"FiaDbhDens")
modeldbh<-ncvar_get(nc,"DbhDens")
modeldbh[modeldbh==-999] <- NA
fian<-ncvar_get(nc,"FiaNPlots")
dbh<-ncvar_get(nc,"Dbh")

fiadbh2 <- matrix(nrow=10,ncol=21)
modeldbh2 <- fiadbh2

fian2 <- apply(fian,1,sum)

for (i in 1:21) {

  fiadbh2[,i] <- apply(fiadbh[,i,,,],1,FUN=weighted.mean,w=fian[i,,,],na.rm=T)
  modeldbh2[,i] <- apply(modeldbh[,i,,,],1,FUN=weighted.mean,w=fian[i,,,],na.rm=T)
}

bins <- numeric(21)
bins[1:5] <- 1
bins[6:7] <- 2
bins[8:21] <- 3

fiadbh3 <- matrix(nrow=10,ncol=4)
modeldbh3 <- fiadbh2

for (i in 1:3) {
  
  fiadbh3[,i] <- apply(fiadbh2[,bins==i],1,FUN=weighted.mean,w=fian2[bins==i],na.rm=T)
  modeldbh3[,i] <- apply(modeldbh2[,bins==i],1,FUN=weighted.mean,w=fian2[bins==i],na.rm=T)
}
fiadbh3[,4] <- apply(fiadbh2,1,FUN=weighted.mean,w=fian2,na.rm=T)
modeldbh3[,4] <- apply(modeldbh2,1,FUN=weighted.mean,w=fian2,na.rm=T)


x11(30,15)
layout(matrix(1:4,nrow=1,ncol=4,byrow=T))
par(mar=c(3,3,1,1))

binlabs <- c("Age 0-50","Age 60-70","Age 80+","All")

for (i in 1:4) {
  
  dbhcomp <- cbind(fiadbh3[,i],modeldbh3[,i])
  
  matplot(dbh-5,
          dbhcomp,
          type="b",
          pch=1,
          lty=1,
          col=c("black","red"),
          log="y",
          ylim=c(5e-2,5e3),
          xlim=c(5,95),
          xlab="DBH midpoint (cm)",
          ylab="Density (/ha)",
          main=binlabs[i]
          )

}


legend("topright",
       c("FIA","CAIN"),
       lty=1,
       pch=1,
       col=c("black","red"),
       bty="n"
)
