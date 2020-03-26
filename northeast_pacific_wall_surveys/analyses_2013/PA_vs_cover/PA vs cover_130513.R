# Robin Elahi
## 13 May 2013
## For NEPWS paper, how important is percent cover vs presence absence data?
## only using RE data

# Load the necessary libraries
library(vegan)
library(MASS)
library(labdsv)

dat<-read.csv("sjc_sessile_130513_mds.csv", header=TRUE, na.strings="NA")
dim(dat)
summary(dat)
colnames(dat)
str(dat)

subset 1009 data (to best match NEPWS sampling)
dat1 <- subset(dat, dat$month=="1009")
dim(dat1)
colnames(dat1)

# create spp matrix (including ON and PG)
abun <- dat1[,13:75]
summary(abun)

# transform matrix
abun.sq <- sqrt(abun)

# presence/absence
abun.pa <- decostand(abun, method="pa")

# create group matrix
group <- dat1[,1:12]

# Perform the NMDS in 2 dimensions on raw data
abun.MDS2 = metaMDS(abun, distance="bray", k=2, engine = 'monoMDS', autotransform=FALSE, noshare=0.1, trymax=40, zerodist='add')
abun.MDS2$stress # stress = 0.225

# Perform the NMDS in 2 dimensions on sqrt data
abun.sq.MDS2 = metaMDS(abun.sq, distance="bray", k=2, engine = 'monoMDS', autotransform=FALSE, noshare=0.1, trymax=40, zerodist='add')
abun.sq.MDS2$stress # stress = 0.236


# Perform the NMDS in 2 dimensions on sqrt data
abun.pa.MDS2 = metaMDS(abun.pa, distance="bray", k=2, engine = 'monoMDS', autotransform=FALSE, noshare=0.1, trymax=40, zerodist='add')
abun.pa.MDS2$stress # stress = 0.241



# plot 3 panels
par(mfrow=c(1,3), pty="s", mar=c(1,1,1,1), oma=c(2,4,2,1))

plot(abun.MDS2$points, col=c(group$site), xaxt='n', yaxt='n', xlab='', ylab='', main="Raw cover")
plot(abun.sq.MDS2$points, col=c(group$site), xaxt='n', yaxt='n', xlab='', ylab='', main="Square root")
plot(abun.pa.MDS2$points, col=c(group$site), xaxt='n', yaxt='n', xlab='', ylab='', main="Presence-Absence")

legend("topright", legend=levels(group$site), pch=1, col=c(1:3), bty='y')


# now do the same for the NEPWS data (pa only)
dat2<-read.csv("NEPWS_channel_130513.csv", header=TRUE, na.strings="NA")
colnames(dat2)

# create spp matrix (including ON and PG)
abun2 <- dat2[,15:87]
summary(abun2)

# create group matrix
group2 <- dat2[,1:14]
group2

# Perform the NMDS in 2 dimensions on raw data
nepws.MDS2 = metaMDS(abun2, distance="bray", k=2, engine = 'monoMDS', autotransform=FALSE, noshare=0.1, trymax=40, zerodist='add')
nepws.MDS2$stress # stress = 0.264

# compare both PA plots

par(mfrow=c(1,2), pty="s", mar=c(1,1,1,1), oma=c(2,4,2,1))

plot(nepws.MDS2$points, col=c(group$site), xaxt='n', yaxt='n', xlab='', ylab='', main="NEPWS")
legend("topright", legend=levels(group2$site), pch=1, col=c(1:3), bty='y')


plot(abun.pa.MDS2$points, col=c(group$site), xaxt='n', yaxt='n', xlab='', ylab='', main="Presence-Absence_%Cover")


