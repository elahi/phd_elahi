# Robin Elahi
# 31 May 2013
# plot centroids of nmds results of mobile quadrat data
# use square root transform (?)

require(vegan) 
require(BiodiversityR) # for CAPdiscrim
require(MASS)
require(labdsv)

dat<-read.csv("nepws_mob5MDS_130531.csv", header=TRUE)
names(dat)
nrow(dat)

colnames(dat) # 10:41
summary(dat)
quad.env <- dat[,1:14]

#create matrix for mds, for mobile data only
mat<-dat[,15:46]
mat
colnames(mat)
ncol(mat) # 32 spp matrix, with mobile abundances
nrow(mat)
summary(mat)
head(mat)

# sqrt transform matrix
mat.sq <- sqrt(mat)
mat.sq

# presence/absence
mat.pa <- decostand(mat, "pa")
mat.pa

# get dissim matrix for untransformed data
mat.dist <- vegdist(mat)
mat.dist
# get dissim matrix for sqrt data
mat.sq.dist <- vegdist(mat.sq)
# get dissimilarity matrix, use bray curtis
mat.pa.dist <- vegdist(mat.pa)


###############################################################

### Do site-level mds for Square root data
### get principal coordinates; use classic (metric) MDS scaling
### what k should i use??? arbitrary choice - 100?
mds.sq <- cmdscale(mat.sq.dist, k=100)
mds.sq

summary(mds.sq)
plot(mds.sq, col=c(quad.env$Region))
names(quad.env)

# one outlier - identify
ordiplot(mds.sq, type="t") # more than one outlier - they are overlapping.  probably because they don't have any mobile taxa (and share the dummy taxon)

### bind
pc.df.sq <- cbind(quad.env, mds.sq)
names(pc.df.sq)
nrow(pc.df.sq)

### now take arithmetic averages of pc's for each site
site.mean.sq <- aggregate(pc.df.sq[,15:114], list(pc.df.sq$Site), mean)
site.mean.sq
names(site.mean.sq)

### create matrix of pc's
mat.pc.sq <- site.mean.sq[,2:101]
mat.pc.sq
### run nMDS using euclidean distance on this matrix, reduce to two dimensions
mds.pc.sq <- metaMDS(mat.pc.sq,distance="euclidean", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
mds.pc.sq$stress # 0.07
mds.pc.sq

### load env data at site level
site.env <- read.csv("site.dat2_edit_120712.csv", header=TRUE)
names(site.env)
site.env <- site.env[,1:12]
mds.pc.sq$points
plot(mds.pc.sq$points, col=c(site.env$region))
ordisurf(mds.pc.sq, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
surf.sq <- ordisurf(mds.pc.sq, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
summary(surf.sq)

###############################################################


### try cluster analysis, perhaps ordicluster?
### cluster analysis using hclust
# set up distance matrix from principle coordinates using euclidean distance (to match mds)
mat.pc.dist <- vegdist(mat.pc.sq, 'euclidean')
clust1<-hclust(mat.pc.dist, 'ward')
clust1
summary(clust1)
plot(clust1, labels=site.env$site)
rect.hclust(clust1, k=2, border='red')

# what is the optimum number of classes?
dat.km <- cascadeKM(mat.pc.sq, 2, 16, criterion="calinski")
dat.km
dat.km$results
dat.km$partition
plot(dat.km, sortq=TRUE)

###############################################################
sq.scrs <- scores(mds.pc.sq, display=c("sites"))
sq.scrs
xlimit.sq <- range(sq.scrs[,1])
ylimit.sq <- range(sq.scrs[,2])

sq.scrs <- cbind(site.env, sq.scrs)

# subset the data for sq
Haro.sq <- (subset(sq.scrs, region == "Haro"))[,13:14]
Channel.sq <- (subset(sq.scrs, region == "Channel"))[,13:14]
Georgia.sq <- (subset(sq.scrs, region == "Georgia"))[,13:14]
Sound.sq <- (subset(sq.scrs, region == "Sound"))[,13:14]
Hood.sq <- (subset(sq.scrs, region == "Hood"))[,13:14]


###############################################################


par(mfrow=c(2,1), pty="s", mar=c(0.5, 0.5, 0.5, 0.5), oma=c(2,1,2,1))

plot(mds.pc.sq$points, type="n", xlab="", ylab="", xaxt='n', yaxt='n')

ordisurf(mds.pc.sq$points, site.env$sedim.per, add=TRUE, labcex=0.8, col="gray40", lty=2)

points(Channel.sq, pch=22, col=1, cex=2.1, lwd=2, bg="white")
points(Haro.sq, pch=24, col=2, cex=2.1, lwd=2, bg="white")
points(Georgia.sq, pch=25, col=3, cex=2.1, lwd=2, bg="white")
points(Sound.sq, pch=23, col=5, cex=2.1, lwd=2, bg="white")
points(Hood.sq, pch=21, col=4, cex=2.1, lwd=2, bg="white")

text(mds.pc.sq$points, labels=site.env$site.no, cex=0.7, font=2)

legend("bottomleft", "a", adj=c(2,1), bty='n', cex=1.5)

leg.txt <- c("Haro", "Channel", "Rosario", "Sound", "Hood")

#legend("bottomright", legend=leg.txt, cex=0.8, pch=c(24,22,25,23,21), col=c(2,1,3,5,4), bty='n', pt.lwd=2)

#text(.4, -0.2, '2D Stress = 0.07', cex=0.8)


# plot cluster analysis
plot(clust1, labels=site.env$site.no, font=1, cex=0.8, sub='', main='', xlab='', hang=-1, las=1)

legend("topright", "b", adj=c(2,0), bty='n', cex=1.5)

points(1,0, pch=21, cex=1.5, col=4, lwd=2)
points(2,0, pch=23, cex=1.5, col=5, lwd=2)
points(3,0, pch=21, cex=1.5, col=4, lwd=2)
points(4,0, pch=21, cex=1.5, col=4, lwd=2)
points(5,0, pch=23, cex=1.5, col=5, lwd=2)
points(6,0, pch=23, cex=1.5, col=5, lwd=2)
points(7,0, pch=22, cex=1.5, col=1, lwd=2)
points(8,0, pch=24, cex=1.5, col=2, lwd=2)
points(9,0, pch=23, cex=1.5, col=5, lwd=2)
points(10,0, pch=22, cex=1.5, col=1, lwd=2)
points(11,0, pch=24, cex=1.5, col=2, lwd=2)
points(12,0, pch=22, cex=1.5, col=1, lwd=2)
points(13,0, pch=25, cex=1.5, col=3, lwd=2)
points(14,0, pch=25, cex=1.5, col=3, lwd=2)
points(15,0, pch=25, cex=1.5, col=3, lwd=2)
points(16,0, pch=24, cex=1.5, col=2, lwd=2)
points(17,0, pch=24, cex=1.5, col=2, lwd=2)
points(18,0, pch=22, cex=1.5, col=1, lwd=2)

box()

rect.hclust(clust1, k=2, border="gray")


###############################################################
# simper
# at the quad level
simp1 <- simper(mat.sq, group=quad.env$Region)
simp1
# what to plot: Call, Troc, Ton, Amp

