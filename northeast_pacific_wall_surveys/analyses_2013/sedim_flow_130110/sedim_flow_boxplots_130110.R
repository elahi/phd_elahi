# Robin Elahi
# 10 June 2013
# Sediment and dissolution boxplots, removing one PG block outlier

# quad level data - sediment
dat<-read.csv("NEPWS_rich_120624.csv", header=TRUE)
names(dat)
par(mfrow=c(2,1))
boxplot(richness ~ Region, data=dat)
# reorder by median richness
rich.med <- with(dat, reorder(Region, -richness, median))
rich.med

boxplot(sedim.per ~ rich.med, data=dat)

#load dissolution data
clod <- read.csv("channel_sound_flow_121104.csv", header=TRUE, na.strings="NA")
names(clod)

# get means + SD
library(plyr)
diss_summary <- ddply(clod, c("Site"), function(df)return(c(Average=mean(df$Diss.day.cm2*1000), SD=sd(df$Diss.day.cm2*1000))))
diss_summary
write.csv(diss_summary, "diss_summary.csv")

###############
# Final 2 panel plot of sediment cover and dissolution data

xlab1 <- c("Haro", "Channel", "Rosario", "Sound", "Hood")
ylab1 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))
plot(1,1, ylab=ylab1)

par(mfcol=c(2,1), mar=c(3,5,1,1), pty="m")

boxplot(sedim.per ~ rich.med, data=dat, ylab="Sediment cover (%)", notch=TRUE, cex.lab=1.4, cex.axis=1.2, las=1, xaxt="n",yaxt="n")
legend("topleft", "a", bty="n", cex=2, adj=c(2,-0.1))
axis(1, labels=xlab1, at=c(1:5), cex.axis=1.2)
axis(2,labels=TRUE, cex.axis=1.2, las=1)

boxplot(Diss.day.cm2*1000 ~ Region, data=clod, ylab=ylab1, notch=FALSE, boxwex=0.4, cex.lab=1.4, cex.axis=1.2, las=1, yaxt="n")
legend("topleft", "b", bty="n", cex=2, adj=c(2,-0.1))
axis(2,labels=TRUE, cex.axis=1.2, las=1)


# LMER because I will remove one value

names(clod)
summary(clod)
loss <- clod$Diss.day.cm2*1000
hist(loss)# p = 0.003
loss.ln <- log(loss)
hist(loss.ln)# p = 0.23

library(lme4)
library(languageR)
library(memisc)

## Inference, region as fixed
mod1<-lmer(loss.ln ~ Region + (1|Site), data=clod, REML=TRUE)
mod1

table1 <- mtable(mod1, summary.stats=FALSE, coef.style="horizontal")
table1
write.mtable(mtable(mod1, summary.stats=FALSE, coef.style="horizontal"))

# check normality
qqnorm(resid(mod1)) # normal quantile-quantile plot of the values
qqline(resid(mod1)) 
# check residuals against fitted values
plot(resid(mod1) ~ fitted(mod1)); abline(h=0) 
shapiro.test(residuals(mod1)) 
# significance tests
mod1.pvals<-pvals.fnc(mod1, nsim=10000)
mod1.pvals$fixed 

# LMER because I will remove one value
loss
clod2 <- clod[-9,]
dim(clod2)
names(clod2)
loss2 <- clod2$Diss.day.cm2*1000
hist(loss2)

## Inference, region as fixed
mod1 <- lmer(loss2 ~ Region + (1|Site), data=clod2, REML=TRUE)
mod1
# check normality
qqnorm(resid(mod1)) # normal quantile-quantile plot of the values
qqline(resid(mod1)) 
# check residuals against fitted values
plot(resid(mod1) ~ fitted(mod1)); abline(h=0) 
shapiro.test(residuals(mod1)) 
# significance tests
mod1.pvals<-pvals.fnc(mod1, nsim=10000)
mod1.pvals$fixed 

# get means + SD
library(plyr)
reg_summary <- ddply(clod2, c("Region"), function(df)return(c(Average=mean(df$Diss.day.cm2*1000), SD=sd(df$Diss.day.cm2*1000))))
reg_summary



# Variance calculation
vc1<-lmer(loss2 ~ (1|Region) + (1|Site), data=clod2, REML=TRUE)
vc1
# partition the variance
VC = as.numeric(VarCorr(vc1))
names(VC) = names(VarCorr(vc1))
VC = c(VC, Resid=attr(VarCorr(vc1), 'sc')^2)
VC

GetPercent = function(x)
{
	rslt = 0
	for(i in 1:length(x))
	{
		rslt[i] = x[i] / sum(x)
	}
	return(rslt)
}


data.frame(Variance=round(VC,3),PercVar= round(GetPercent(VC),5)*100, row.names=names(VC))

       # Variance PercVar
# Site      0.025  25.393
# Region    0.067  67.444
# Resid     0.007   7.163

# significance test for random effects
vc1.pvals<-pvals.fnc(vc1, nsim=10000)
vc1.pvals$random


###############
# New 2 panel plot of sediment cover and dissolution data

xlab1 <- c("Haro", "Channel", "Rosario", "Sound", "Hood")
ylab1 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))
plot(1,1, ylab=ylab1)

par(mfcol=c(2,1), mar=c(3,5,1,1), pty="m")

boxplot(sedim.per ~ rich.med, data=dat, ylab="Sediment cover (%)", notch=TRUE, cex.lab=1.4, cex.axis=1.2, las=1, xaxt="n",yaxt="n")
legend("topleft", "a", bty="n", cex=2, adj=c(2,-0.1))
axis(1, labels=xlab1, at=c(1:5), cex.axis=1.2)
axis(2,labels=TRUE, cex.axis=1.2, las=1)

boxplot(Diss.day.cm2*1000 ~ Region, data=clod2, ylab=ylab1, notch=FALSE, boxwex=0.4, cex.lab=1.4, cex.axis=1.2, las=1, yaxt="n")
legend("topleft", "b", bty="n", cex=2, adj=c(2,-0.1))
axis(2,labels=TRUE, cex.axis=1.2, las=1)


# how many sd's away is the outlier?
clod.mean <- mean(clod$Diss.day.cm2*1000)
clod.sd <- sd(clod$Diss.day.cm2*1000)

clod.mean # 0.984
clod.sd # 0.340
# outlier is 2.121
2.121 - 0.984 #1.137 is the difference
1.137/0.34 # 3.344 SD's away from the mean

# get new means + SD
library(plyr)
diss_summary <- ddply(clod2, c("Site"), function(df)return(c(Average=mean(df$Diss.day.cm2*1000), SD=sd(df$Diss.day.cm2*1000))))
diss_summary
write.csv(diss_summary, "diss_summary_130610.csv")


