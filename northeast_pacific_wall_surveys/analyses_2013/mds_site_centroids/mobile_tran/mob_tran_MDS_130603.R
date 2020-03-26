# Robin Elahi
# 3 June 2013
# plot centroids of nmds results of mobile transect data abundance
# removed rows without any species (22 rows)
# removed species occuring in less than 5% of plots (3/11 species)

require(vegan) 
require(BiodiversityR) # for CAPdiscrim
require(MASS)
require(labdsv)

dat <- read.csv("NEPWS_mob_tran_MDS_120823.csv", header=TRUE)
names(dat)
nrow(dat)

# create envt file
env <- dat[,1:6]
#create spp matrix
spp <- dat[,7:14]
# 5% minimum occurence for inclusion in MDS
spp1 <- vegtab(spp, minval=5)
ncol(spp1)
names(spp1) # the same
# get dissimilarity matrix, use bray curtis
spp.dist <- vegdist(spp)
spp.dist

##########################################################################
# metaMDS of transects
meta1 <- metaMDS(spp,distance="bray", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
meta1 #stress = 0.18
plot(meta1, type="t")  
meta1$points # this contains the range needed for plotting both species and sites

# extact species values
plot(meta1, display="species", type="text")
spp.val <- plot(meta1, display="species", type="n")  # extract the species values
spp.val
plot(spp.val$species)
spp.points <- spp.val$species
spp.points
spp.names <- row.names(spp.points)
spp.names

# extract transect values
scrs <- scores(meta1, display=c("sites"))
scrs
xlimit <- range(scrs[,1])
ylimit <- range(scrs[,2])
names(env)

# cbind the group columns with mds scores
group.scrs <- cbind(env, scrs)
colnames(group.scrs)
summary(group.scrs)

#write.csv(group.scrs, "group.scrs.csv")

# subset the data

Haro <- (subset(group.scrs, region == "Haro"))[,7:8]
Channel <- (subset(group.scrs, region == "Channel"))[,7:8]
Rosario <- (subset(group.scrs, region == "Georgia"))[,7:8]
Sound <- (subset(group.scrs, region == "Sound"))[,7:8]
Hood <- (subset(group.scrs, region == "Hood"))[,7:8]

#Haro: green triangles (pch 24)
#Channel: black squares (pch 22)
#Georgia: red inverse triangles (pch 25)
#Sound: light blue diamonds (pch 23)
#Hood: blue circles (pch 21)

# plot species
plot(spp.val$species, xlim=range(meta1$points[,1]), ylim=range(meta1$points[,2]), type='n', xlab="", ylab="", xaxt='n', yaxt='n')
text(spp.points[,1], spp.points[,2], labels=spp.names, cex=0.5)

text(-1.1, -1, 'Stress = 0.18')

# plot transects
points(Channel, pch=22, col=1, cex=1, lwd=1, bg="white")
points(Haro, pch=24, col=2, cex=1, lwd=1, bg="white")
points(Rosario, pch=25, col=3, cex=1, lwd=1, bg="white")
points(Sound, pch=23, col=5, cex=1, lwd=1, bg="white")
points(Hood, pch=21, col=4, cex=1, lwd=1, bg="white")

text(mds.pc$points, labels=group$site.no, cex=0.7, font=2)

legend("bottomleft", "a", adj=c(2,1), bty='n', cex=1.5)

leg.txt <- c("Haro", "Channel", "Rosario", "Sound", "Hood")

legend("bottomright", legend=leg.txt, cex=0.8, pch=c(24,22,25,23,21), col=c(2,1,3,5,4), bty='n', pt.lwd=2)

text(.27, -0.2, '2D Stress = 0.11', cex=0.8)

##########################################################################
# CENTROIDS
### get principal coordinates; use classic (metric) MDS scaling
### what k should i use??? maximum = 8, because 8 species
mds1 <- cmdscale(spp.dist, k=8)
mds1
summary(mds1)
plot(mds1, col=c(env$region))

### bind
pc.df <- cbind(env, mds1)
names(pc.df)
nrow(pc.df)

### now take arithmetic averages of pc's for each site
site.mean.pc <- aggregate(pc.df[,7:14], list(pc.df$site), mean)
site.mean.pc
names(site.mean.pc)
# write.csv(site.mean.pc, "site.mean.pc.csv")

### create matrix of pc's
mat.pc <- site.mean.pc[,2:9]
mat.pc
### run nMDS using euclidean distance on this matrix, reduce to two dimensions
mds.pc <- metaMDS(mat.pc,distance="euclidean", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
mds.pc$stress # 0.11
mds.pc
plot(mds.pc)

### load env data at site level

site.env <- read.csv("site.dat2_edit_120712.csv", header=TRUE)
names(site.env)
site.env <- site.env[,1:12]
mds.pc$points
plot(mds.pc$points, col=c(site.env$region))
ordisurf(mds.pc, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
summary(ordisurf(mds.pc, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2))

### try cluster analysis, perhaps ordicluster?
### cluster analysis using hclust
# set up distance matrix from principle coordinates using euclidean distance (to match mds)
mat.pc.dist <- vegdist(mat.pc, 'euclidean')
clust1<-hclust(mat.pc.dist, 'ward')
clust1
summary(clust1)
plot(clust1, labels=site.env$site)
rect.hclust(clust1, k=3, border='red')


clust2<-hclust(mat.pc.dist, 'single')
clust2
summary(clust2)
plot(clust2, labels=site.env$site)
rect.hclust(clust1, k=2, border='red')

# what is the optimum number of classes?
dat.km <- cascadeKM(mat.pc, 2, 16, criterion="calinski")
dat.km
plot(dat.km, sortq=TRUE)


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

write.csv(group.scrs, "group.scrs.csv")
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

text(mds.pc$points, labels=group$site.no, cex=0.7, font=2)

legend("bottomleft", "a", adj=c(2,1), bty='n', cex=1.5)

leg.txt <- c("Haro", "Channel", "Rosario", "Sound", "Hood")

legend("bottomright", legend=leg.txt, cex=0.8, pch=c(24,22,25,23,21), col=c(2,1,3,5,4), bty='n', pt.lwd=2)

text(.27, -0.2, '2D Stress = 0.11', cex=0.8)

# plot cluster analysis
plot(clust1, labels=group$site.no, font=1, cex=0.8, sub='', main='', xlab='', hang=-1, las=1)

legend("topleft", "b", adj=c(2,0), bty='n', cex=1.5)

#Haro: red triangles (pch 24)
#Channel: black squares (pch 22)
#Georgia: green inverse triangles (pch 25)
#Sound: light blue diamonds (pch 23)
#Hood: blue circles (pch 21)

points(1,0, pch=22, cex=1.5, col=1, lwd=2) # 
points(2,0, pch=24, cex=1.5, col=2, lwd=2)
points(3,0, pch=25, cex=1.5, col=3, lwd=2)
points(4,0, pch=24, cex=1.5, col=2, lwd=2)
points(5,0, pch=25, cex=1.5, col=3, lwd=2)
points(6,0, pch=22, cex=1.5, col=1, lwd=2)
points(7,0, pch=24, cex=1.5, col=2, lwd=2)
points(8,0, pch=22, cex=1.5, col=1, lwd=2)
points(9,0, pch=24, cex=1.5, col=2, lwd=2)
points(10,0, pch=22, cex=1.5, col=1, lwd=2)
points(11,0, pch=21, cex=1.5, col=4, lwd=2)
points(12,0, pch=25, cex=1.5, col=3, lwd=2)
points(13,0, pch=21, cex=1.5, col=4, lwd=2)
points(14,0, pch=21, cex=1.5, col=4, lwd=2)
points(15,0, pch=23, cex=1.5, col=5, lwd=2)
points(16,0, pch=23, cex=1.5, col=5, lwd=2)
points(17,0, pch=23, cex=1.5, col=5, lwd=2)
points(18,0, pch=23, cex=1.5, col=5, lwd=2)

box()

rect.hclust(clust1, k=3, border="gray")
###############################################################
# simper
# at the tran level
names(env)
simp1 <- simper(spp, group=env$region)
simp1
#Hen, PCA, SFR, PHE, DIM

simp2 <- simper(spp, group=env$site)
simp2
Hen, PCA, SFR, PHE, DIM
