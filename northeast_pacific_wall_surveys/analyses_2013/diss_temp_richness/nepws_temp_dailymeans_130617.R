# Robin Elahi
# 17 June 2013
# Plot temperature data for NEPWS sites

library(plotrix)

dat<-read.csv("temp daily_means.csv", header=TRUE, na.strings="NA")
head(dat)
dim(dat)
colnames(dat)
summary(dat$date1)
summary(dat)

# create vector of dates in proper format - this must be labeled as a number in excel to work (i.e., not 12-jul-2008)
dates<-as.Date(dat$date1, origin="1899-12-30")
dates
summary(dates)

# merge dates with data
dat1<-cbind(dat, dates)
nrow(dat1)
ncol(dat1)
colnames(dat1)
head(dat1)
summary(dat1$dates)
summary(dat1$dates)

# load dissolution and temp data together

dat2<-read.csv("diss_temp means_130610.csv", header=TRUE, na.strings="NA")
dim(dat2)
names(dat2)
summary(dat2)
dat2$site.no

xlab1 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))

lm1 <- lm(temp.mean ~ diss.mean, data=dat2)
summary(lm1)
anova(lm1)

plot(temp.mean ~ diss.mean, data=dat2, col=c(1,5)[region], pch=c(1,5)[region], cex=3, lwd=2, las=1, ylab=my.ylab, xlab=xlab1, ylim=c(10.6, 11.5))
abline(lm1, lwd=2, col="darkgray")
text(dat2$diss.mean, dat2$temp.mean, pch=dat2$site.no)

legend(0.6,11, legend=c("y = - 1.0656x + 12.0922", expression(R^{2} * " = 0.83"),expression(F[1*","*6] * " = 29.7"), expression(italic(P) * " = 0.001")), bty = "n")


############################################################################

# plot 2panel
my.ylab = expression(paste("Temperature ", "(", degree,"C",")"))
xlab1 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))

leg.txt.Sound <- c("9. Rosario Wall", "10. Humphrey Head", "11. Frost Island", "12. Willow Island")

leg.txt.Channel <- c("5. O'Neal", "6. Shady Cove", "7. Point George", "8. Turn Island")

par(mfcol=c(2,1), mar=c(3,5,1,1), oma=c(1,0,0,0))

plot(Rosario ~ dates, data=dat1, type="l", col=5, ylab=my.ylab, cex.lab=1.2, xlab="", las=1, ylim=c(10.2,14.8), lwd=2)
axis(1, dates, labels=FALSE, lwd=0.5)

points(Humphrey ~ dates, data=dat1, type="l", col=5, lwd=2, lty=2)
points(Frost ~ dates, data=dat1, type="l", col=5, lwd=2, lty=3)
points(Willow ~ dates, data=dat1, type="l", col=5, lwd=2, lty=4)

points(ON ~ dates, data=dat1, type="l", col=1, lwd=2, lty=1)
points(SC ~ dates, data=dat1, type="l", col=1, lwd=2, lty=2)
points(PG ~ dates, data=dat1, type="l", col=1, lwd=2, lty=3)
points(Turn ~ dates, data=dat1, type="l", col=1, lwd=2, lty=4)

legend(15570, 14, legend=leg.txt.Sound, col=5, bty="n", lty=c(1,2,3,4), lwd=2, cex=0.8)
text(15575.5, 14, "Sound", font=2)

legend(15550, 14, legend=leg.txt.Channel, col=1, bty="n", lty=c(1,2,3,4), lwd=2, cex=0.8)
text(15555.5, 14, "Channel", font=2)

legend("topright", "a ", bty="n", cex=1.8, adj=0)

points(15553, 14.5, pch=21, cex=3, lwd=2) # full moon

points(15561, 14.5, pch=21, bg=1, cex=3, lwd=2) # new moon
points(15559.85, 14.5, pch=15, col="white", cex=4, lwd=2) # halve it
points(15561, 14.5, pch=21, cex=3, lwd=2) # new moon

points(15569, 14.5, pch=21, bg=1, cex=3, lwd=2) # new moon

points(15576, 14.5, pch=21, bg=1, cex=3, lwd=2) # new moon
points(15577.15, 14.5, pch=15, col="white", cex=4, lwd=2) # halve it
points(15576, 14.5, pch=21, cex=3, lwd=2) # new moon

points(15583, 14.5, pch=21, cex=3, lwd=2) # full moon


#summary(dat1$dates)

# > summary(dat$date1)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 41120   41130   41140   41140   41150   41160 

# panel B

plot(temp.mean ~ diss.mean, data=dat2, col=c(1,5)[region], pch=c(1,5)[region], cex=3, lwd=2, las=1, ylab=my.ylab, xlab="", ylim=c(10.6, 11.5), cex.lab=1.2)
abline(lm1, lwd=2, col="darkgray")
text(dat2$diss.mean, dat2$temp.mean, labels=dat2$site.no)

legend(0.6,10.9, legend=c(expression(R^{2} * " = 0.83"),expression(F[1*","*6] * " = 29.7"), expression(italic(P) * " = 0.001")), bty = "n")

legend("topright", "b ", bty="n", cex=1.8, adj=0)
mtext(side=1, xlab1, line=2.5, cex=1.2)


############################################################################
