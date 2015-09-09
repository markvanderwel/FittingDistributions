
#run code from Graphs.R up to about L156 to generate totdata

totdata$Ba2 <- totdata$ModelledBa
totdata$G2 <- totdata$GrowthModelled
totdata$M2 <- totdata$MortModelled
totdata$R2 <- totdata$IngModelled
totdata$Age2 <- totdata$Age

totdata$LogBa2 <- log(totdata$Ba2)
totdata$LogG2 <- log(totdata$G2)
totdata$LogM2 <- log(totdata$M2)
totdata$LogR2 <- log(totdata$R2)

totdata$BaGt5 <- totdata$ModelledBa > 5

totdata$LogBa2[totdata$Ba2<1] <- NA
totdata$LogG2[totdata$Ba2<1] <- NA
totdata$LogM2[totdata$Ba2<1] <- NA
totdata$LogR2[totdata$Ba2<1] <- NA

totdata$StdLogBa2 <- NA
totdata$StdLogG2 <- NA
totdata$StdLogM2 <- NA
totdata$StdLogR2 <- NA
totdata$StdAge2 <- NA

for (iPft in unique(totdata$Pft)) {
  totdata$StdLogBa2[totdata$Pft == iPft] <- scale(totdata$LogBa2[totdata$Pft == iPft])
  totdata$StdLogG2[totdata$Pft == iPft] <- scale(totdata$LogG2[totdata$Pft == iPft])
  totdata$StdLogM2[totdata$Pft == iPft] <- scale(totdata$LogM2[totdata$Pft == iPft])
  totdata$StdLogR2[totdata$Pft == iPft] <- scale(totdata$LogR2[totdata$Pft == iPft])
  totdata$StdAge2[totdata$Pft == iPft] <- scale(totdata$Age2[totdata$Pft == iPft])

#   #no scaling
#   totdata$StdLogBa2[totdata$Pft == iPft] <- (totdata$LogBa2[totdata$Pft == iPft])
#   totdata$StdLogG2[totdata$Pft == iPft] <- (totdata$LogG2[totdata$Pft == iPft])
#   totdata$StdLogM2[totdata$Pft == iPft] <- (totdata$LogM2[totdata$Pft == iPft])
#   totdata$StdLogR2[totdata$Pft == iPft] <- (totdata$LogR2[totdata$Pft == iPft])
#   totdata$StdAge2[totdata$Pft == iPft] <- (totdata$Age2[totdata$Pft == iPft])
    
}

CM.StdLogBa2 <- aggregate(totdata$StdLogBa2,by=list(totdata$Lon,totdata$Lat,totdata$Pft),FUN=mean,na.rm=T)
CM.StdLogG2 <- aggregate(totdata$StdLogG2,by=list(totdata$Lon,totdata$Lat,totdata$Pft),FUN=mean,na.rm=T)
CM.StdLogM2 <- aggregate(totdata$StdLogM2,by=list(totdata$Lon,totdata$Lat,totdata$Pft),FUN=mean,na.rm=T)
CM.StdLogR2 <- aggregate(totdata$StdLogR2,by=list(totdata$Lon,totdata$Lat,totdata$Pft),FUN=mean,na.rm=T)

totdata$StdLogBa2.CM <- CM.StdLogBa2$x[match(paste(totdata$Lon,totdata$Lat,totdata$Pft),paste(CM.StdLogG2[,1],CM.StdLogG2[,2],CM.StdLogG2[,3]))]
totdata$StdLogG2.CM <- CM.StdLogG2$x[match(paste(totdata$Lon,totdata$Lat,totdata$Pft),paste(CM.StdLogG2[,1],CM.StdLogG2[,2],CM.StdLogG2[,3]))]
totdata$StdLogM2.CM <- CM.StdLogM2$x[match(paste(totdata$Lon,totdata$Lat,totdata$Pft),paste(CM.StdLogG2[,1],CM.StdLogG2[,2],CM.StdLogG2[,3]))]
totdata$StdLogR2.CM <- CM.StdLogR2$x[match(paste(totdata$Lon,totdata$Lat,totdata$Pft),paste(CM.StdLogG2[,1],CM.StdLogG2[,2],CM.StdLogG2[,3]))]

totdata$StdLogBa2.Res <- totdata$StdLogBa2 - totdata$StdLogBa2.CM
totdata$StdLogG2.Res <- totdata$StdLogG2 - totdata$StdLogG2.CM
totdata$StdLogM2.Res <- totdata$StdLogM2 - totdata$StdLogM2.CM
totdata$StdLogR2.Res <- totdata$StdLogR2 - totdata$StdLogR2.CM

save.coef <- matrix(nrow=6,ncol=8)
rownames(save.coef) <- unique(totdata$Pft)
for (iPft in unique(totdata$Pft)) {
 print(iPft)
 lm.mod <- lm(StdLogBa2 ~ StdLogG2.CM + StdLogG2.Res + StdLogM2.CM + StdLogM2.Res + StdLogR2.CM + StdLogR2.Res + StdAge2, data=totdata, subset=(totdata$Pft==iPft))
 #logistic regression for Ba > 5 (range limit defn from earlier paper)
 #lm.mod <- glm(BaGt5 ~ StdLogG2.CM + StdLogG2.Res + StdLogM2.CM + StdLogM2.Res + StdLogR2.CM + StdLogR2.Res + StdAge2, data=totdata, subset=(totdata$Pft==iPft),family=binomial)
 print(summary(lm.mod))
 save.coef[iPft,] <- lm.mod$coefficients
}
colnames(save.coef) <- names(lm.mod$coefficients)
save.coef
barplot(t(save.coef[,-1]),beside=T,col=c("red","red","green","green","blue","blue","yellow"))

#among-cell plots
layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
for (iPft in unique(totdata$Pft)) {
  with(subset(totdata,Pft==iPft),plot(StdLogG2.CM,StdLogBa2.CM,xlim=c(-2,2)))
  title(main=iPft)
}
layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
for (iPft in unique(totdata$Pft)) {
  with(subset(totdata,Pft==iPft),plot(StdLogM2.CM,StdLogBa2.CM,xlim=c(-2,2)))
  title(main=iPft)
}
layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
for (iPft in unique(totdata$Pft)) {
  with(subset(totdata,Pft==iPft),plot(StdLogR2.CM,StdLogBa2.CM,xlim=c(-2,2)))
  title(main=iPft)
}
#within-cell plots
layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
for (iPft in unique(totdata$Pft)) {
  with(subset(totdata,Pft==iPft),plot(StdLogG2.Res,StdLogBa2.Res,xlim=c(-2,2)))
  title(main=iPft)
}
layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
for (iPft in unique(totdata$Pft)) {
  with(subset(totdata,Pft==iPft),plot(StdLogM2.Res,StdLogBa2.Res,xlim=c(-2,2)))
  title(main=iPft)
}
layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
for (iPft in unique(totdata$Pft)) {
  with(subset(totdata,Pft==iPft),plot(StdLogR2.Res,StdLogBa2.Res,xlim=c(-2,2)))
  title(main=iPft)
}



layout(1)
