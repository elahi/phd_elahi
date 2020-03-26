# Robin Elahi
# 19 June 2013
# Sediment  boxplots

# quad level data - sediment
dat<-read.csv("NEPWS_rich_121114.csv", header=TRUE)
names(dat)
par(mfrow=c(2,1))
boxplot(richness ~ Region, data=dat)
# reorder by median richness
rich.med <- with(dat, reorder(Region, -richness, median))
rich.med

boxplot(sedim.per ~ rich.med, data=dat)

###############
# Final plot of sediment cover 

xlab1 <- c("Haro", "Channel", "Rosario", "Sound", "Hood")
ylab1 <- expression(paste("Daily mass loss (mg ", cm^-2, ")"))
plot(1,1, ylab=ylab1)

par(mfcol=c(1,1), mar=c(3,5,1,1), pty="m")

boxplot(sedim.per ~ rich.med, data=dat, ylab="Sediment cover (%)", notch=TRUE, cex.lab=1.4, cex.axis=1.2, las=1, xaxt="n",yaxt="n")
#legend("topleft", "a", bty="n", cex=2, adj=c(2,-0.1))
axis(1, labels=xlab1, at=c(1:5), cex.axis=1.2)
axis(2,labels=TRUE, cex.axis=1.2, las=1)

#############################################

# LMER for sediment

names(dat)
summary(dat)
library(lme4)
library(languageR)
library(memisc)

library(car)
sed.log <- logit(dat$sedim)
hist(sed.log)
hist(dat$sedim)
sed.int <- round(dat$sedim.per, digits=0)
sed.int

## Inference, region as fixed
mod1<-lmer(sedim.per ~ Region + (1|Site) + (1|Transect), data=dat, REML=TRUE)
mod1<-lmer(sedim ~ Region + (1|Site) + (1|Transect), data=dat, REML=TRUE, family=binomial(link="logit"))
mod1<-lmer(sed.log ~ Reg + (1|Site) + (1|Transect), data=dat, REML=TRUE) # use this for analysis

# check residuals against fitted values
plot(resid(mod1) ~ fitted(mod1)); abline(h=0) 
summary(mod1)

table1 <- mtable(mod1, summary.stats=FALSE, coef.style="horizontal")
table1
write.mtable(mtable(mod1, summary.stats=FALSE, coef.style="horizontal"))

Vcov <- vcov(mod1, useScale=FALSE)
Vcov
betas <- fixef(mod1)
betas
se <- sqrt(diag(Vcov))
se
zval <- betas/se
zval
pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
cbind(betas, se, zval, pval)

# check normality
qqnorm(resid(mod1)) # normal quantile-quantile plot of the values
qqline(resid(mod1)) 
# check residuals against fitted values
plot(resid(mod1) ~ fitted(mod1)); abline(h=0) 
shapiro.test(residuals(mod1)) 
# significance tests
mod1.pvals<-pvals.fnc(mod1, nsim=10000)
mod1.pvals$fixed 




# Variance calculation
vc1<-lmer(sed.log ~ (1|Region) + (1|Site) + (1|Transect), data=dat, REML=TRUE)
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
# Transect    0.332   6.188
# Site        0.986  18.389
# Region      3.761  70.169
# Resid       0.282   5.254

# significance test for random effects
vc1.pvals<-pvals.fnc(vc1, nsim=10000)
vc1.pvals$random




