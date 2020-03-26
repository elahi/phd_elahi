# Robin Elahi
# 25 Sept 2013
# calculate per occurrence (mean + SD) by region, using site means

dat<-read.csv("nepws correlations.csv", header=TRUE)
names(dat)


# dissolution against noaa current
lm1 <- lm(diss.day.cm.mean ~ aug.mean, data=dat)
lm2 <- lm(diss.day.cm.mean ~ aug.max, data=dat)
summary(lm1)
summary(lm2)

plot(lm1)

summary(dat)

#trying nls (power model)
lm6 <- lm(sess.q.mean ~ aug.mean, data=dat)
summary(lm6)
nls1 <- nls(sess.q.mean ~ a*aug.mean^b, data=dat, start=list(a=10, b=1))
summary(nls1)
nls1
AIC(lm6, nls1)
leg.txt6<- c(expression(R^{2} * " = 0.55"),expression(F[1*","*16] * " = 19.69"), expression(italic(P) * " < 0.001"))

newdata = new=data.frame(aug.mean=seq(0,2.15,len=200))
plot(sess.q.mean ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,20))
lines(newdata$aug.mean, predict(nls1, new=newdata, type="res"), col="red")


#############################################

# plot 2panel
xlab1 <- expression(paste("Predicted current speed (m ", s^-1, ")"))
ylab1 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))

mean.txt <- c("Daily average current", expression(R^{2} * " = 0.77"),expression(F[1*","*6] * " = 19.9"), expression(italic(P) * " = 0.004"))
max.txt <- c("Daily maximum current", expression(R^{2} * " = 0.80"),expression(F[1*","*6] * " = 23.6"), expression(italic(P) * " = 0.003"))

par(mar=c(4,5,1,1))
plot(diss.day.cm.mean ~ aug.mean, data=dat, las=1, xlim=c(0.6,1.6), ylim=c(0.6, 1.4), xlab=xlab1, ylab=ylab1, cex=1.5, pch=21, lwd=1.5)
#curve(1.0741*x + 0.001, from=0.6, to=1.2, add=TRUE, lwd=1.5, lty=2) # old reg using 1.64 for PG
curve(0.9292*x + 0.1038, from=0.656, to=1.174, add=TRUE, lwd=1.5, lty=2) 

points(diss.day.cm.mean ~ aug.max, data=dat, cex=1.5, pch=21, bg=1)
curve(0.83268*x -0.06369, from=0.906, to=1.523, add=TRUE, lwd=1.5, lty=1)

#legend("topleft", legend=c("Daily average current", "Daily maximum current"), cex=1, pch=21, col=1, pt.bg=c(1,0), bty='n', pt.lwd=1.5, pt.cex=1.5)

legend(0.6, 1.2, mean.txt, bty="n")
legend(1.2, 0.8, max.txt, bty="n")

.9292*0.62 - .1038
# 0.62 mg mass loss is equal to 0.47 mean current speed; or Harney Channel is 0.656 mean current speed
.9292*1.40 - .1038
# 1.4 mg mass loss is equal to 1.20 mean current speed; or Pear point is 1.17

#############################################
# plot A
lm3 <- lm(sess.q.mean ~ diss.day.cm.mean, data=dat)
summary(lm3) #ns
leg.txt3<- c(expression(R^{2} * " = 0.38"),expression(F[1*","*6] * " = 3.64"), expression(italic(P) * " = 0.105"))

lm4 <- lm(sess.Sobs ~ diss.day.cm.mean, data=dat)
summary(lm4) #ms
leg.txt4<- c(expression(R^{2} * " = 0.40"),expression(F[1*","*6] * " = 3.94"), expression(italic(P) * " = 0.09"))

lm5 <- lm(sess.Schao ~ diss.day.cm.mean, data=dat)
summary(lm5) #sig
leg.txt5<- c(expression(R^{2} * " = 0.77"),expression(F[1*","*6] * " = 20.66"), expression(italic(P) * " = 0.004"))

lm6 <- lm(sess.q.mean ~ aug.mean, data=dat)
summary(lm6)
leg.txt6<- c(expression(R^{2} * " = 0.55"),expression(F[1*","*16] * " = 19.69"), expression(italic(P) * " < 0.001"))

lm7 <- lm(sess.Sobs ~ aug.mean, data=dat)
summary(lm7)
leg.txt7<- c(expression(R^{2} * " = 0.51"),expression(F[1*","*16] * " = 16.89"), expression(italic(P) * " < 0.001"))

lm8 <- lm(sess.Schao ~ aug.mean, data=dat)
summary(lm8)
leg.txt8<- c(expression(R^{2} * " = 0.36"),expression(F[1*","*16] * " = 8.99"), expression(italic(P) * " = 0.008"))

#############################################
ylab1 <- expression(italic(S)[obs])
#plot(1,1, ylab=ylab1)
ylab2 <- expression(italic(S)[Chao2])
#plot(1,1, ylab=ylab2)

xlab1 <- expression(paste("Predicted current speed (m ", s^-1, ")"))
xlab2 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))

par(mfcol=c(3,2), mar=c(1,3,1,2), oma=c(5,3,3,0), pty="m")

plot(sess.q.mean ~ diss.day.cm.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,20))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Richness", side=2, line=3, cex=1.5)
legend("topright", "a", bty="n", cex=2, adj=c(0.5,-0.1))
legend("bottomright", leg.txt3, bty="n", cex=1.3)

plot(sess.Sobs ~ diss.day.cm.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,80))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext(ylab1, side=2, line=3, cex=1.5)
legend("topright", "b", bty="n", cex=2, adj=c(0.5,-0.1))
#curve(31.03*x + 11.01, from=0.6, to=1.4, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt4, bty="n", cex=1.3)

plot(sess.Schao ~ diss.day.cm.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=TRUE, cex.axis=1.5)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext(ylab2, side=2, line=3, cex=1.5)
legend("topright", "c", bty="n", cex=2, adj=c(0.5,-0.1))
mtext(xlab2, side=1, line=4, cex=1.2)
curve(55.95*x - 5.007, from=0.6, to=1.4, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt5, bty="n", cex=1.3)

plot(sess.q.mean ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,20))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels= TRUE, cex.axis=1.5, las=1)
#mtext("Richness", side=2, line=3, cex=1.5)
legend("topright", "d", bty="n", cex=2, adj=c(0.5,-0.1))
curve(6.316*x + 2.588, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt6, bty="n", cex=1.3)
abline(v=0.65, col="darkgray") # or
abline(v=1.174, col="darkgray") # 
#abline(v=0.47, col="darkgray") # 
#abline(v=1.12, col="darkgray") # 

plot(sess.Sobs ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,80))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels= TRUE, cex.axis=1.5, las=1)
#mtext(ylab1, side=2, line=3, cex=1.5)
legend("topright", "e", bty="n", cex=2, adj=c(0.5,-0.1))
curve(19.757*x + 21.097, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt7, bty="n", cex=1.3)
abline(v=0.65, col="darkgray") # or
abline(v=1.174, col="darkgray") # 

plot(sess.Schao ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels= TRUE, cex.axis=1.5)
axis(2,labels= TRUE, cex.axis=1.5, las=1)
#mtext(ylab2, side=2, line=3, cex=1.5)
legend("topright", "f", bty="n", cex=2, adj=c(0.5,-0.1))
mtext(xlab1, side=1, line=4, cex=1.2)
curve(24.856*x + 25.447, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt8, bty="n", cex=1.3)
abline(v=0.65, col="darkgray") # or
abline(v=1.174, col="darkgray") # 

#############################################
# Stats for solitary inverts

lm9 <- lm(TETR*100 ~ aug.mean, data=dat)
summary(lm9) # sig
leg.txt9<- c(expression(R^{2} * " = 0.64"),expression(F[1*","*16] * " = 28.34"), expression(italic(P) * " < 0.001"))

lm10 <- lm(Metr*100 ~ aug.mean, data=dat)
summary(lm10) # sig
leg.txt10<- c(expression(R^{2} * " = 0.43"),expression(F[1*","*16] * " = 12.07"), expression(italic(P) * " = 0.003"))

lm11 <- lm(PSCH*100 ~ aug.mean, data=dat)
summary(lm11) # sig
leg.txt11<- c(expression(R^{2} * " = 0.61"),expression(F[1*","*16] * " = 24.85"), expression(italic(P) * " < 0.001"))

lm12 <- lm(BAEL*100 ~ aug.mean, data=dat)
summary(lm12) # sig
leg.txt12<- c(expression(R^{2} * " = 0.41"),expression(F[1*","*16] * " = 11.16"), expression(italic(P) * " = 0.004"))

lm13 <- lm(POMA*100 ~ aug.mean, data=dat)
summary(lm13) # NS
leg.txt13<- c(expression(R^{2} * " = 0.15"),expression(F[1*","*16] * " = 2.77"), expression(italic(P) * " = 0.115"))

lm14 <- lm(BACR*100 ~ aug.mean, data=dat)
summary(lm14) # NS
leg.txt14<- c(expression(R^{2} * " = 0.14"),expression(F[1*","*16] * " = 2.60"), expression(italic(P) * " = 0.127"))

lm15 <- lm(SCUN*100 ~ aug.mean, data=dat)
summary(lm15) # sig
leg.txt15<- c(expression(R^{2} * " = 0.49"),expression(F[1*","*16] * " = 15.45"), expression(italic(P) * " = 0.001"))

lm16 <- lm(CNFI*100 ~ aug.mean, data=dat)
summary(lm16) # sig
leg.txt16<- c(expression(R^{2} * " = 0.15"),expression(F[1*","*16] * " = 2.72"), expression(italic(P) * " = 0.118"))



#############################################
# original
par(mfcol=c(3,2), mar=c(1,3,1,3), oma=c(5,5,1,0), pty="m")

plot(TETR*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Terebretalia", side=2, line=3, cex=1.1, font=3)
legend("topright", "a", bty="n", cex=2, adj=c(0.5,-0.1))
curve(46.110*x - 15.007, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt9, bty="n", cex=1.3)

plot(Metr*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Metridium", side=2, line=3, cex=1.1, font=3)
legend("topright", "b", bty="n", cex=2, adj=c(0.5,-0.1))
curve(40.04*x - 17.34, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt10, bty="n", cex=1.3)

plot(PSCH*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=TRUE, cex.axis=1.5)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Psolus", side=2, line=3, cex=1.1, font=3)
legend("topright", "c", bty="n", cex=2, adj=c(0.5,-0.1))
curve(46.752*x - 21.126, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt11, bty="n", cex=1.3)

plot(BAEL*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Balanophyllia", side=2, line=3, cex=1.1, font=3)
legend("topright", "d", bty="n", cex=2, adj=c(0.5,-0.1))
curve(43.185*x - 5.486, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt12, bty="n", cex=1.3)

plot(POMA*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Pododesmus", side=2, line=3, cex=1.1, font=3)
legend("topright", "e", bty="n", cex=2, adj=c(0.5,-0.1))
legend("topright", leg.txt13, bty="n", cex=1.3, adj=0.5)

plot(BACR*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=TRUE, cex.axis=1.5)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Balanus", side=2, line=3, cex=1.1, font=3)
legend("topright", "f", bty="n", cex=2, adj=c(0.5,-0.1))
#mtext(xlab1, side=1, line=4, cex=1.2)
legend("topright", leg.txt14, bty="n", cex=1.3, adj=0.5)

mtext(xlab1, side=1, line=2.5, cex=1.2, outer=TRUE)
mtext("Occurence in quadrats (%)", side=2, line=2.5, cex=1.2, outer=TRUE)

#############################################

# 2nd version
par(mfcol=c(3,2), mar=c(1,3,1,3), oma=c(5,5,1,0), pty="m")

plot(BAEL*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Balanophyllia", side=2, line=3, cex=1.1, font=3)
legend("topright", "a", bty="n", cex=2, adj=c(0.5,-0.1))
curve(43.185*x - 5.486, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt12, bty="n", cex=1.3)

plot(TETR*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Terebretalia", side=2, line=3, cex=1.1, font=3)
legend("topright", "b", bty="n", cex=2, adj=c(0.5,-0.1))
curve(46.110*x - 15.007, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt9, bty="n", cex=1.3)

plot(Metr*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=TRUE, cex.axis=1.5, las=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Metridium", side=2, line=3, cex=1.1, font=3)
legend("topright", "c", bty="n", cex=2, adj=c(0.5,-0.1))
curve(40.04*x - 17.34, from=0.45, to=2.15, add=TRUE, lwd=1.5, lty=1)
legend("bottomright", leg.txt10, bty="n", cex=1.3)

plot(POMA*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Pododesmus", side=2, line=3, cex=1.1, font=3)
legend("topright", "d", bty="n", cex=2, adj=c(0.5,-0.1))
legend("topright", leg.txt13, bty="n", cex=1.3, adj=0.5)

plot(BACR*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=FALSE, cex.axis=1.5)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Balanus", side=2, line=3, cex=1.1, font=3)
legend("topright", "e", bty="n", cex=2, adj=c(0.5,-0.1))
#mtext(xlab1, side=1, line=4, cex=1.2)
legend("topright", leg.txt14, bty="n", cex=1.3, adj=0.5)

plot(CNFI*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))
axis(1, labels=TRUE, cex.axis=1.5)
axis(2,labels=TRUE, cex.axis=1.5, las=1)
mtext("Cnemidocarpa", side=2, line=3, cex=1.1, font=3)
legend("topright", "f", bty="n", cex=2, adj=c(0.5,-0.1))
#mtext(xlab1, side=1, line=4, cex=1.2)
legend("topright", leg.txt16, bty="n", cex=1.3, adj=0.5)

mtext(xlab1, side=1, line=2.5, cex=1.2, outer=TRUE)
mtext("Occurence in quadrats (%)", side=2, line=2.5, cex=1.2, outer=TRUE)


#############################################

require(graphics)
with(cars, scatter.smooth(speed, dist))
cars.lo <- loess(dist ~ speed, cars)
predict(cars.lo, data.frame(speed = seq(5, 30, 1)), se = TRUE)
# to allow extrapolation
cars.lo2 <- loess(dist ~ speed, cars,
  control = loess.control(surface = "direct"))
predict(cars.lo2, data.frame(speed = seq(5, 30, 1)), se = TRUE)

#############################################
with(dat, scatter.smooth(aug.max, POMA*100, degree=1))

scatter.smooth(dat$TETR*100 ~ dat$aug.mean, degree=1)
scatter.smooth(dat$sedim.mean ~ dat$aug.mean, degree=1)

l1 <- loess(TETR*100 ~ aug.mean, data=dat)
summary(l1)

plot(TETR*100 ~ aug.mean, data=dat, las=1, xlab="", ylab="", cex=1.5, pch=21, lwd=1.5, cex.lab=1.5, cex.axis=1, xaxt='n', yaxt='n', bg=1, ylim=c(0,100))

lines(dat$aug.mean, predict(l1))

###
names(dat)

plot(aug.mean ~ ann.curr.mean, data=dat)
summary(lm(ann.curr.mean ~ aug.mean, data=dat))

summary(lm(ENCA ~ aug.mean, data=dat))
summary(lm(ENNA ~ aug.mean, data=dat))
summary(lm(BREN ~ aug.mean, data=dat))
summary(lm(Rfil ~ aug.mean, data=dat))

par(mfrow=c(4,4), mar=c(4,5,0.5,0.5))

plot(diss.day.cm.mean ~ aug.mean, data=dat)
plot(sess.Sobs ~ aug.mean, data=dat)
plot(sess.Schao ~ aug.mean, data=dat)
plot(sess.q.mean ~ aug.mean, data=dat)
plot(mob.Sobs ~ aug.mean, data=dat)
plot(mob.Schao ~ aug.mean, data=dat)
plot(mob.q.mean ~ aug.mean, data=dat)
plot(sedim.mean ~ aug.mean, data=dat)
plot(Metr ~ aug.mean, data=dat)
plot(BAEL ~ aug.mean, data=dat)
plot(HYSP ~ aug.mean, data=dat)
plot(Agla ~ aug.mean, data=dat)
plot(Abie ~ aug.mean, data=dat)
plot(META ~ aug.mean, data=dat)
plot(POMA ~ aug.mean, data=dat)
plot(BACR ~ aug.mean, data=dat)
plot(PSCH ~ aug.mean, data=dat)
plot(TETR ~ aug.mean, data=dat)
plot(SCUN ~ aug.mean, data=dat)
plot(CNFI ~ aug.mean, data=dat)



plot(ENCA ~ aug.mean, data=dat)
plot(ENNA ~ aug.mean, data=dat)
plot(BREN ~ aug.mean, data=dat)
abline(lm(BREN ~ aug.mean, data=dat))


