# Robin Elahi
# 3 June 2013
#plot call, trochids, amp, and ton

require(vegan) 
require(BiodiversityR) # for CAPdiscrim
require(MASS)
require(labdsv)
library(plotrix)
library(sciplot)

dat<-read.csv("NEPWS_combined_130603.csv", header=TRUE)
names(dat)
nrow(dat)
colnames(dat)


xlab1 <- c("Haro", "Channel", "Rosario", "Sound", "Hood")

#practice
bargraph.CI(Region.code, CLI, data=dat, lc=TRUE, uc=TRUE, ylab="Density", xaxt='n', las=2, col=c(2,1,3,5,4), err.width=0.03, main=expression(paste(italic(Calliostoma)," ",italic(ligatum))))

# using SD
bargraph.CI(Region.code, CLI, data=dat, lc=TRUE, uc=TRUE, ylab=expression(paste("Density (no. ", m^-2, ")")), xaxt='n', las=2, col=c(2,1,3,5,4), err.width=0.03, ylim=c(0,30), ci.fun= function(x) c(mean(x)-sd(x), mean(x)+sd(x)))


# by site
bargraph.CI(Site, CLI, data=dat, lc=TRUE, uc=TRUE, ylab=expression(paste("Density (no. ", m^-2, ")")), xaxt='t', las=2, col=c(2,1,3,5,4), err.width=0.03, ylim=c(0,30))


################

par(mfrow=c(2,2), mar=c(3,5,1,1), oma=c(2,2,2,2))

bargraph.CI(Region.code, CLI, data=dat, lc=TRUE, uc=TRUE, ylab=expression(paste("Density (no. ", m^-2, ")")), xaxt='n', las=2, col=c(2,1,3,5,4), err.width=0.03, ylim=c(0,13), cex.lab=1.4, cex.axis=1.2)
title(main="Calliostoma ligatum", font.main=3)
box()

bargraph.CI(Region.code, Tro, data=dat, lc=TRUE, uc=TRUE, ylab="", xaxt='n', las=2, col=c(2,1,3,5,4), err.width=0.03, ylim=c(0,20), cex.lab=1.4, cex.axis=1.2)
title(main="Trochid snails", font.main=1)
box()

bargraph.CI(Region.code, Ton, data=dat, lc=TRUE, uc=TRUE, ylab=expression(paste("Density (no. ", m^-2, ")")), xaxt='n', las=2, col=c(2,1,3,5,4), err.width=0.03, ylim=c(0,10), cex.lab=1.4, cex.axis=1.2)
title(main="Tonicella", font.main=3)
box()

bargraph.CI(Region.code, Amp, data=dat, lc=TRUE, uc=TRUE, ylab="", xaxt='n', las=2, col=c(2,1,3,5,4), err.width=0.03, ylim=c(0,13), cex.lab=1.4, cex.axis=1.2)
title(main="Amphissa", font.main=3)
box()




