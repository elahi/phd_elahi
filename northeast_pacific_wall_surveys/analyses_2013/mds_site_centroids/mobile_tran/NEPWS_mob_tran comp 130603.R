# Robin Elahi
# 3 June 203
# plot centroids of nmds results of mobile transect
# compare abundances, sqrt, P/A

require(vegan) 
require(BiodiversityR) # for CAPdiscrim
require(MASS)
require(labdsv)

dat <- read.csv("NEPWS_mob_tran_MDS_120823.csv", header=TRUE)
colnames(dat) 
summary(dat)

# create envt file
tran.env <- dat[,1:6]
names(tran.env)

#create matrix for mds, for mobile data only
mat<-dat[,7:14]
mat
colnames(mat)
ncol(mat) # 8 taxa
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
### Do site-level mds for untransformed data
### get principal coordinates; use classic (metric) MDS scaling
### what k should i use??? arbitrary choice - 100?
mds.un <- cmdscale(mat.dist, k=75)
mds.un

summary(mds.un)
plot(mds.un, col=c(tran.env$region))
names(tran.env)

# one outlier - identify
ordiplot(mds.un, type="t") # more than one outlier - they are overlapping.  probably because they don't have any mobile taxa (and share the dummy taxon)

### bind
pc.df.un <- cbind(tran.env, mds.un)
names(pc.df.un)
nrow(pc.df.un)

### now take arithmetic averages of pc's for each site
site.mean.un <- aggregate(pc.df.un[,7:40], list(pc.df.un$site), mean)
site.mean.un
dim(site.mean.un)
names(site.mean.un)

### create matrix of pc's
mat.pc.un <- site.mean.un[,2:34]
mat.pc.un
### run nMDS using euclidean distance on this matrix, reduce to two dimensions
mds.pc.un <- metaMDS(mat.pc.un,distance="euclidean", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
mds.pc.un$stress # 0.10
mds.pc.un

### load env data at site level
site.env <- read.csv("site.dat2_edit_120712.csv", header=TRUE)
names(site.env)
site.env <- site.env[,1:12]
mds.pc.un$points
plot(mds.pc.un$points, col=c(site.env$region))
ordisurf(mds.pc.un, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
surf.un <- ordisurf(mds.pc.un, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
summary(surf.un)

###############################################################

### Do site-level mds for Square root data
### get principal coordinates; use classic (metric) MDS scaling
### what k should i use??? arbitrary choice - 100?
mds.sq <- cmdscale(mat.sq.dist, k=75)
mds.sq

summary(mds.sq)
plot(mds.sq, col=c(tran.env$region))
names(tran.env)

# one outlier - identify
ordiplot(mds.sq, type="t") # more than one outlier - they are overlapping.  probably because they don't have any mobile taxa (and share the dummy taxon)

### bind
pc.df.sq <- cbind(tran.env, mds.sq)
names(pc.df.sq)
nrow(pc.df.sq)

### now take arithmetic averages of pc's for each site
site.mean.sq <- aggregate(pc.df.sq[,7:38], list(pc.df.sq$site), mean)
site.mean.sq
names(site.mean.sq)

### create matrix of pc's
mat.pc.sq <- site.mean.sq[,2:33]
mat.pc.sq
### run nMDS using euclidean distance on this matrix, reduce to two dimensions
mds.pc.sq <- metaMDS(mat.pc.sq,distance="euclidean", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
mds.pc.sq$stress # 0.09
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

### Do site-level mds for PA data
### get principal coordinates; use classic (metric) MDS scaling
### what k should i use??? arbitrary choice - 100?
mds.pa <- cmdscale(mat.pa.dist, k=75)
mds.pa

summary(mds.pa)
plot(mds.pa, col=c(tran.env$region))
names(tran.env)

# one outlier - identify
ordiplot(mds.pa, type="t") # more than one outlier - they are overlapping.  probably because they don't have any mobile taxa (and share the dummy taxon)

### bind
pc.df.pa <- cbind(tran.env, mds.pa)
names(pc.df.pa)
dim(pc.df.pa)
nrow(pc.df.pa)

### now take arithmetic averages of pc's for each site
site.mean.pa <- aggregate(pc.df.pa[,7:39], list(pc.df.pa$site), mean)
site.mean.pa
names(site.mean.pa)

### create matrix of pc's
mat.pc.pa <- site.mean.pa[,2:33]
mat.pc.pa
### run nMDS using euclidean distance on this matrix, reduce to two dimensions
mds.pc.pa <- metaMDS(mat.pc.pa,distance="euclidean", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
mds.pc.pa$stress # 0.09
mds.pc.pa

### load env data at site level
site.env <- read.csv("site.dat2_edit_120712.csv", header=TRUE)
names(site.env)
site.env <- site.env[,1:12]
mds.pc.pa$points
plot(mds.pc.pa$points, col=c(site.env$region))
ordisurf(mds.pc.pa, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
surf.pa <- ordisurf(mds.pc.pa, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
summary(surf.pa)


### 3 panel plot

un.scrs <- scores(mds.pc.un, display=c("sites"))
un.scrs
xlimit.un <- range(un.scrs[,1])
ylimit.un <- range(un.scrs[,2])

sq.scrs <- scores(mds.pc.sq, display=c("sites"))
sq.scrs
xlimit.sq <- range(sq.scrs[,1])
ylimit.sq <- range(sq.scrs[,2])

pa.scrs <- scores(mds.pc.pa, display=c("sites"))
pa.scrs
xlimit.pa <- range(pa.scrs[,1])
ylimit.pa <- range(pa.scrs[,2])

# cbind the site.env columns with mds scores

un.scrs <- cbind(site.env, un.scrs)
sq.scrs <- cbind(site.env, sq.scrs)
pa.scrs <- cbind(site.env, pa.scrs)

colnames(un.scrs)

# subset the data for un
Haro.un <- (subset(un.scrs, region == "Haro"))[,13:14]
Channel.un <- (subset(un.scrs, region == "Channel"))[,13:14]
Georgia.un <- (subset(un.scrs, region == "Georgia"))[,13:14]
Sound.un <- (subset(un.scrs, region == "Sound"))[,13:14]
Hood.un <- (subset(un.scrs, region == "Hood"))[,13:14]

# subset the data for sq
Haro.sq <- (subset(sq.scrs, region == "Haro"))[,13:14]
Channel.sq <- (subset(sq.scrs, region == "Channel"))[,13:14]
Georgia.sq <- (subset(sq.scrs, region == "Georgia"))[,13:14]
Sound.sq <- (subset(sq.scrs, region == "Sound"))[,13:14]
Hood.sq <- (subset(sq.scrs, region == "Hood"))[,13:14]

# subset the data for pa
Haro.pa <- (subset(pa.scrs, region == "Haro"))[,13:14]
Channel.pa <- (subset(pa.scrs, region == "Channel"))[,13:14]
Georgia.pa <- (subset(pa.scrs, region == "Georgia"))[,13:14]
Sound.pa <- (subset(pa.scrs, region == "Sound"))[,13:14]
Hood.pa <- (subset(pa.scrs, region == "Hood"))[,13:14]


### 

par(mfrow=c(1,3), pty="s", mar=c(0.5, 0.5, 0.5, 0.5), oma=c(2,1,2,1))

# panel A - untranformed data
plot(mds.pc.un$points, type="n", xlab="", ylab="", xaxt='n', yaxt='n')

ordisurf(mds.pc.un$points, site.env$sedim.per, add=TRUE, labcex=0.8, col="gray40", lty=2)

points(Channel.un, pch=22, col=1, cex=2.1, lwd=2, bg="white")
points(Haro.un, pch=24, col=2, cex=2.1, lwd=2, bg="white")
points(Georgia.un, pch=25, col=3, cex=2.1, lwd=2, bg="white")
points(Sound.un, pch=23, col=5, cex=2.1, lwd=2, bg="white")
points(Hood.un, pch=21, col=4, cex=2.1, lwd=2, bg="white")

text(mds.pc.un$points, labels=site.env$site.no, cex=0.7, font=2)

title(main="Untransformed", line=0.5)

# panel B - sqrt transformed data
plot(mds.pc.sq$points, type="n", xlab="", ylab="", xaxt='n', yaxt='n')

ordisurf(mds.pc.sq$points, site.env$sedim.per, add=TRUE, labcex=0.8, col="gray40", lty=2)

points(Channel.sq, pch=22, col=1, cex=2.1, lwd=2, bg="white")
points(Haro.sq, pch=24, col=2, cex=2.1, lwd=2, bg="white")
points(Georgia.sq, pch=25, col=3, cex=2.1, lwd=2, bg="white")
points(Sound.sq, pch=23, col=5, cex=2.1, lwd=2, bg="white")
points(Hood.sq, pch=21, col=4, cex=2.1, lwd=2, bg="white")

text(mds.pc.sq$points, labels=site.env$site.no, cex=0.7, font=2)

title(main="Square root transformed", line=0.5)

# panel C - presence absence
plot(mds.pc.pa$points, type="n", xlab="", ylab="", xaxt='n', yaxt='n')

ordisurf(mds.pc.pa$points, site.env$sedim.per, add=TRUE, labcex=0.8, col="gray40", lty=2)

points(Channel.pa, pch=22, col=1, cex=2.1, lwd=2, bg="white")
points(Haro.pa, pch=24, col=2, cex=2.1, lwd=2, bg="white")
points(Georgia.pa, pch=25, col=3, cex=2.1, lwd=2, bg="white")
points(Sound.pa, pch=23, col=5, cex=2.1, lwd=2, bg="white")
points(Hood.pa, pch=21, col=4, cex=2.1, lwd=2, bg="white")

text(mds.pc.pa$points, labels=site.env$site.no, cex=0.7, font=2)

title(main="Presence-Absence", line=0.5)


#legend("bottomleft", "a", adj=c(2,1), bty='n', cex=1.5)

#leg.txt <- c("Haro", "Channel", "Rosario", "Sound", "Hood")

#legend("bottomright", legend=leg.txt, cex=0.8, pch=c(24,22,25,23,21), col=c(2,1,3,5,4), bty='n', pt.lwd=2)

#text(.4, -0.2, '2D Stress = 0.07', cex=0.8)










### try cluster analysis, perhaps ordicluster?
### cluster analysis using hclust
# set up distance matrix from principle coordinates using euclidean distance (to match mds)
mat.pc.dist <- vegdist(mat.pc, 'euclidean')
clust1<-hclust(mat.pc.dist, 'ward')
clust1
summary(clust1)
plot(clust1, labels=site.env$site)
rect.hclust(clust1, k=2, border='red')

clust2<-hclust(mat.pc.dist, 'single')
clust2
summary(clust2)
plot(clust2, labels=site.env$site)
rect.hclust(clust1, k=2, border='red')

# what is the optimum number of classes?
dat.km <- cascadeKM(mat.pc, 2, 16, criterion="calinski")
dat.km
dat.km$results
dat.km$partition
plot(dat.km, sortq=TRUE)

         # 2 groups 3 groups 4 groups  5 groups  6 groups  7 groups
# SSE      1.648616 1.279437 1.058067 0.9023451 0.7824762 0.6644489
# calinski 5.561725 5.523432 5.132216 4.7519123 4.4143298 4.2967008

### final 2 panel plot
scrs <- scores(mds.pc, display=c("sites"))
scrs
xlimit <- range(scrs[,1])
ylimit <- range(scrs[,2])
group <- site.env[,1:12]
names(site.env)
# cbind the group columns with mds scores
group.scrs <- cbind(site.env, scrs)
colnames(group.scrs)
summary(group.scrs)

#write.csv(group.scrs, "group.scrs.csv")

# subset the data
Haro <- (subset(group.scrs, region == "Haro"))[,13:14]
Channel <- (subset(group.scrs, region == "Channel"))[,13:14]
Georgia <- (subset(group.scrs, region == "Georgia"))[,13:14]
Sound <- (subset(group.scrs, region == "Sound"))[,13:14]
Hood <- (subset(group.scrs, region == "Hood"))[,13:14]

#Haro: green triangles (pch 24)
#Channel: black squares (pch 22)
#Georgia: red inverse triangles (pch 25)
#Sound: light blue diamonds (pch 23)
#Hood: blue circles (pch 21)


###############################################################
# plot 2 panels
par(mfrow=c(2,1), pty="s", mar=c(0.5, 0.5, 0.5, 0.5), oma=c(2,1,2,1))

# panel A
#plot.new()
#plot.window(xlim = xlimit, ylim = ylimit, asp=1)

plot(mds.pc$points, type="n", xlab="", ylab="", xaxt='n', yaxt='n')

ordisurf(mds.pc$points, group$sedim.per, add=TRUE, labcex=0.8, col="gray40", lty=2)

points(Channel, pch=22, col=1, cex=2.1, lwd=2, bg="white")
points(Haro, pch=24, col=2, cex=2.1, lwd=2, bg="white")
points(Georgia, pch=25, col=3, cex=2.1, lwd=2, bg="white")
points(Sound, pch=23, col=5, cex=2.1, lwd=2, bg="white")
points(Hood, pch=21, col=4, cex=2.1, lwd=2, bg="white")

# points(Haro, pch=24, col="gray40", bg="red", cex=2.1)
# points(Channel, pch=22, col="gray40", bg="yellow", cex=2.1)
# points(Georgia, pch=25, col="gray40", bg="green",cex=2.1)
# points(Sound, pch=23, col="gray40", bg="lightblue1", cex=2.1)
# points(Hood, pch=21, col="blue", bg="white", cex=2.1)

text(mds.pc$points, labels=group$site.no, cex=0.7, font=2)

legend("bottomleft", "a", adj=c(2,1), bty='n', cex=1.5)

leg.txt <- c("Haro", "Channel", "Rosario", "Sound", "Hood")

legend("bottomright", legend=leg.txt, cex=0.8, pch=c(24,22,25,23,21), col=c(2,1,3,5,4), bty='n', pt.lwd=2)

text(.4, -0.2, '2D Stress = 0.07', cex=0.8)

# plot cluster analysis
plot(clust1, labels=group$site.no, font=1, cex=0.8, sub='', main='', xlab='', hang=-1, las=1)

legend("topleft", "b", adj=c(2,0), bty='n', cex=1.5)

points(1,0, pch=25, cex=1.5, col=3, lwd=2)
points(2,0, pch=25, cex=1.5, col=3, lwd=2)
points(3,0, pch=25, cex=1.5, col=3, lwd=2)
points(4,0, pch=24, cex=1.5, col=2, lwd=2)
points(5,0, pch=24, cex=1.5, col=2, lwd=2)
points(6,0, pch=22, cex=1.5, col=1, lwd=2)
points(7,0, pch=24, cex=1.5, col=2, lwd=2)
points(8,0, pch=22, cex=1.5, col=1, lwd=2)
points(9,0, pch=22, cex=1.5, col=1, lwd=2)
points(10,0, pch=24, cex=1.5, col=2, lwd=2)
points(11,0, pch=21, cex=1.5, col=4, lwd=2)
points(12,0, pch=21, cex=1.5, col=4, lwd=2)
points(13,0, pch=23, cex=1.5, col=5, lwd=2)
points(14,0, pch=21, cex=1.5, col=4, lwd=2)
points(15,0, pch=23, cex=1.5, col=5, lwd=2)
points(16,0, pch=23, cex=1.5, col=5, lwd=2)
points(17,0, pch=22, cex=1.5, col=1, lwd=2)
points(18,0, pch=23, cex=1.5, col=5, lwd=2)

box()

rect.hclust(clust1, k=2, border="gray")
###############################################################


### Significance testing
### 1.  Permdisp
### 2.  Permanova
