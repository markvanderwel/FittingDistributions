x <- read.table("D:\\code\\UF\\src\\FitCainToPftBa\\bin\\Debug\\workspace\\CainFitToPftBa_MCMC_bayes_list.txt",header=T)
x[1:10,1:10]
dim(x)
rg <- 17375:17482
colnames(x)[rg]
means <- colMeans(x[,rg])
sds <- apply(x[,rg],2,FUN=sd)

plot(means)
arrows(x0=1:108,y0=means-sds,y1=means+sds,length=0)

rT <- seq(0,1,length.out=100)
rP <- rT
TP <- expand.grid(rT,rP)

colnames(TP)<-c("T","P")
X <- cbind(1,TP$T,TP$P,TP$T^2,TP$T*TP$P,TP$P^2)

plot(NULL,NULL,xlim=c(0,1),ylim=c(0,6))
for(i in seq(1,108,by=6)) {
  y <- tapply(X%*%means[i:(i+5)],TP$P,FUN=mean)
  lines(unique(TP$P),y,type="l")
}
  

