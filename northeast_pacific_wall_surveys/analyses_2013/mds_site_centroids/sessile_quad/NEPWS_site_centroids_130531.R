# Robin Elahi
# 31 May 2013
# plot centroids of nmds results of sess data only
# P/A only!  

require(vegan) 
require(BiodiversityR) # for CAPdiscrim
require(MASS)
require(labdsv)

dat0<-read.csv("NEPWS_combined_120626.csv", header=TRUE)
head(dat0)
names(dat0)
colnames(dat0) # 15:87 all data
summary(dat0)
quad.env <- dat0[,1:14]

#create matrix for mds, does not include space or sediment cover
mat0<-dat0[,15:87]
colnames(mat0)
ncol(mat0) # 132 spp mat0rix, with mobile abundances
nrow(mat0)
summary(mat0)
head(mat0)

# convert to presence absence
mat0.pa <- decostand(mat0, "pa")
# get dissimilarity matrix, use bray curtis
mat0.pa.dist <- vegdist(mat0.pa)


### get principal coordinates; use classic (metric) MDS scaling
### what k should i use??? arbitrary choice - 100?
mds1 <- cmdscale(mat0.pa.dist, k=100)
mds1
summary(mds1)

### bind
pc.df <- cbind(quad.env, mds1)
names(pc.df)
nrow(pc.df)

### now take arithmetic averages of pc's for each site
site.mean.pc <- aggregate(pc.df[,15:86], list(pc.df$Site), mean)
site.mean.pc
names(site.mean.pc)

### create matrix of pc's
mat.pc <- site.mean.pc[,2:73]
mat.pc
### run nMDS using euclidean distance on this matrix, reduce to two dimensions
mds.pc <- metaMDS(mat.pc,distance="euclidean", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
mds.pc$stress # 0.11
mds.pc

### load env data at site level

site.env <- read.csv("site.dat2_edit_120712.csv", header=TRUE)
names(site.env)
site.env <- site.env[,1:12]
mds.pc$points
plot(mds.pc$points, col=c(site.env$region))
ordisurf(mds.pc, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
surf1 <- ordisurf(mds.pc, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
summary(surf1)


### try cluster analysis, perhaps ordicluster?
### cluster analysis using hclust
# set up distance matrix from principle coordinates using euclidean distance (to match mds)
mat.pc.dist <- vegdist(mat.pc, 'euclidean')
clust1<-hclust(mat.pc.dist, 'ward')
clust1
summary(clust1)
plot(clust1, labels=group$site)
rect.hclust(clust1, k=2, border='red')

clust2<-hclust(mat.pc.dist, 'single')
clust2
summary(clust2)
plot(clust2, labels=group$site)
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
names(group)
# cbind the group columns with mds scores
group.scrs <- cbind(group, scrs)
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

text(.4, -0.2, '2D Stress = 0.11', cex=0.8)

# plot cluster analysis
plot(clust1, labels=group$site.no, font=1, cex=0.8, sub='', main='', xlab='', hang=-1, las=1)

legend("topleft", "b", adj=c(2,0), bty='n', cex=1.5)

# 22square, 25downtri; 24uptri; 

points(1,0, pch=22, cex=1.5, col=1, lwd=2)
points(2,0, pch=24, cex=1.5, col=2, lwd=2)
points(3,0, pch=24, cex=1.5, col=2, lwd=2)
points(4,0, pch=24, cex=1.5, col=2, lwd=2)
points(5,0, pch=22, cex=1.5, col=1, lwd=2)
points(6,0, pch=22, cex=1.5, col=1, lwd=2)
points(7,0, pch=24, cex=1.5, col=2, lwd=2)
points(8,0, pch=22, cex=1.5, col=1, lwd=2)
points(9,0, pch=25, cex=1.5, col=3, lwd=2)
points(10,0, pch=25, cex=1.5, col=3, lwd=2)
points(11,0, pch=25, cex=1.5, col=3, lwd=2)
points(12,0, pch=21, cex=1.5, col=4, lwd=2)
points(13,0, pch=21, cex=1.5, col=4, lwd=2)
points(14,0, pch=23, cex=1.5, col=5, lwd=2)
points(15,0, pch=21, cex=1.5, col=4, lwd=2)
points(16,0, pch=23, cex=1.5, col=5, lwd=2)
points(17,0, pch=23, cex=1.5, col=5, lwd=2)
points(18,0, pch=23, cex=1.5, col=5, lwd=2)

box()

rect.hclust(clust1, k=2, border="gray")
###############################################################


### Significance testing
### 1.  Permdisp
### 2.  Permanova
