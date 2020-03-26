# Robin Elahi
# 19 June 2013

# Plot tidal flux for relevant stations 
dat2<-read.csv("noaa current3 130609.csv", header=TRUE, na.strings="NA")


xlab1 <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M','N')
ylab0 <- expression(paste("Predicted current speed (m ", s^-1, ")"))
ylab1 <- expression(paste("Average current speed (m ", s^-1, ")"))
ylab2 <- expression(paste("Maximum current speed (m ", s^-1, ")"))

#par(mfrow=c(2,1), mar=c(1,1,0.5,0.5), oma=c(4,5,0,0))

################################################

par(mfcol=c(2,1), mar=c(1,5,0.5,0.5), oma=c(3,2,0,0))

boxplot(ms.mean ~ Station.code, data=dat2, las=1, notch=FALSE, boxfill=c(2,2,2,2,0,0,0,3,3,3,5,5,4,4), ylim=c(0,2), xaxt="n", ylab="Daily average")
axis(side=1, at=c(1:14), labels=FALSE)
legend("topleft", "a", bty="n", cex=1.8, adj=2)

leg.txt <- c("Haro (A-D)", "Channel (E-G)", "Rosario (H-J)", "Sound (K,L)", "Hood (M,N)")
legend("topright", legend=leg.txt, cex=1, pch=22, col=1, pt.bg=c(2,0,3,5,4), bty='n', pt.lwd=1.5, pt.cex=1.5)

boxplot(ms.max ~ Station.code, data=dat2, las=1, notch=FALSE, boxfill=c(2,2,2,2,0,0,0,3,3,3,5,5,4,4), ylim=c(0,2), ylab="Daily maximum")
#axis(side=1, at=c(1:14), labels=xlab1)
legend("topleft", "b", bty="n", cex=1.8, adj=2)
mtext("Current station", side=1, outer=TRUE, line=1.5, cex=1.2)

mtext(ylab0, side=2, outer=TRUE, line=0, cex=1.2)

################################################

# one panel of daily averages
par(mar=c(5,5,1,1), pty="m")

boxplot(ms.mean ~ Station.code, data=dat2, las=1, notch=FALSE, boxfill=c(2,2,2,2,0,0,0,3,3,3,5,5,4,4), xaxt="n", ylab=ylab0, cex.lab=1.4)
axis(side=1, at=c(1:14), labels=xlab1)

leg.txt <- c("Haro (A-D)", "Channel (E-G)", "Rosario (H-J)", "Sound (K,L)", "Hood (M,N)")
legend("topright", legend=leg.txt, cex=1, pch=22, col=1, pt.bg=c(2,0,3,5,4), bty='n', pt.lwd=1.5, pt.cex=1.5)

mtext("Current station", side=1, outer=FALSE, line=3, cex=1.4)

################################################

#load dissolution data
clod <- read.csv("channel_sound_flow_121104.csv", header=TRUE, na.strings="NA")
names(clod)

boxplot(Diss.day.cm2*1000 ~ Region, data=clod, ylab=ylab1, notch=FALSE, boxwex=0.4, cex.lab=1.4, cex.axis=1.2, las=1, yaxt="n")
legend("topleft", "b", bty="n", cex=2, adj=c(2,-0.1))
axis(2,labels=TRUE, cex.axis=1.2, las=1)


################################################

# load correlation data
corr <- read.csv("nepws correlations.csv", header=TRUE)

# dissolution against noaa current
lm1 <- lm(diss.day.cm.mean ~ aug.mean, data=corr)
summary(lm1)
plot(lm1)
leg.txt1<- c(expression(R^{2} * " = 0.77"),expression(F[1*","*6] * " = 19.86"), expression(italic(P) * " = 0.004"))

xlab1 <- expression(paste("Predicted current speed (m ", s^-1, ")"))
ylab1 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))

mean.txt <- c("Daily average current", expression(R^{2} * " = 0.77"),expression(F[1*","*6] * " = 19.9"), expression(italic(P) * " = 0.004"))

par(mar=c(4,5,1,1))
plot(diss.day.cm.mean ~ aug.mean, data=corr, las=1, xlim=c(0.6,1.2), ylim=c(0.5, 1.5), xlab=xlab1, ylab=ylab1, cex=1.5, pch=21, lwd=1.5, bg=1)
curve(0.9292*x + 0.1038, from=0.656, to=1.174, add=TRUE, lwd=1.5, lty=1) 
legend("bottomright", leg.txt1, bty="n", cex=1.3)

################################################
# 3 panel plot
layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), heights=c(2,1), respect=TRUE)
layout.show(3)

par(mar=c(5,5,1,0.5), pty="m")
boxplot(ms.mean ~ Station.code, data=dat2, las=1, notch=FALSE, boxfill=c(2,2,2,2,0,0,0,3,3,3,5,5,4,4), xaxt="n", ylab=ylab0, cex.lab=1.4, cex.axis=1)
xlab1 <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M','N')
axis(side=1, at=c(1:14), labels=xlab1)

leg.txt <- c("Haro (A-D)", "Channel (E-G)", "Rosario (H-J)", "Sound (K,L)", "Hood (M,N)")
legend("topright", legend=leg.txt, cex=0.8, pch=22, col=1, pt.bg=c(2,0,3,5,4), bty='n', pt.lwd=1.5, pt.cex=1.5)

mtext("Current station", side=1, outer=FALSE, line=3, cex=1.2)
legend("topleft", "a", bty="n", cex=2, adj=c(2,-0.1))


par(mar=c(4,5,0,0.5))

ylab1 <- "Daily mass loss"
ylab2 <- expression(paste("(mg ", cm^-2, ")"))

boxplot(Diss.day.cm2*1000 ~ Region, data=clod, xlab="Seascape",ylab="", notch=FALSE, boxwex=0.4, cex.lab=1, cex.axis=1.05, las=1, boxfill=c(0,5), ylim=c(0.5, 2.5))
legend("topright", "b", bty="n", cex=2, adj=c(0,-0.1))
mtext(ylab1, side=2, outer=FALSE, line=4.5, cex=1)
mtext(ylab2, side=2, outer=FALSE, line=2.7, cex=1)


par(mar=c(4,5,0,0.5))
xlab2 <- expression(paste("Predicted current speed (m ", s^-1, ")"))
plot(diss.day.cm.mean ~ aug.mean, data=corr, las=1, xlim=c(0.6,1.2), ylim=c(0.5, 1.5), xlab="Predicted current speed", ylab="Daily mass loss", cex=1.2, pch=21, lwd=1.5, bg=1, cex.axis=1.05, cex.lab=1.05)
curve(0.9292*x + 0.1038, from=0.656, to=1.174, add=TRUE, lwd=1.5, lty=1) 
legend("bottomright", leg.txt1, bty="n", cex=1)
legend("topleft", "c", bty="n", cex=2, adj=c(2.5,-0.1))

################################################




