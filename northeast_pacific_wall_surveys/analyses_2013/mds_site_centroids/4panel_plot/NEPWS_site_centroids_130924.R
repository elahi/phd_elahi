# Robin Elahi
# 24 sept 2013, modified for a 2 panel plot
# plot centroids of nmds results of sess and mob quadrat data
# P/A only!  

require(vegan) 
require(BiodiversityR) # for CAPdiscrim
require(MASS)
require(labdsv)

dat0<-read.csv("NEPWS_combined_120626.csv", header=TRUE)
head(dat0)
names(dat0)
colnames(dat0) # 15:146 all data
summary(dat0)
quad.env <- dat0[,1:14]

#create matrix for mds, does not include space or sediment cover
mat0<-dat0[,15:146]
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
plot(mds1, col=c(region0))

### bind
pc.df <- cbind(quad.env, mds1)
names(pc.df)
nrow(pc.df)

### now take arithmetic averages of pc's for each site
site.mean.pc <- aggregate(pc.df[,15:112], list(pc.df$Site), mean)
site.mean.pc
names(site.mean.pc)
# write.csv(site.mean.pc, "site.mean.pc.csv")

### create matrix of pc's
mat.pc <- site.mean.pc[,2:99]
mat.pc
### run nMDS using euclidean distance on this matrix, reduce to two dimensions
mds.pc <- metaMDS(mat.pc,distance="euclidean", engine='monoMDS', autotransform=FALSE,noshare=0.1,k=2,tol=1e-5,trymax=40, zerodist='add')
mds.pc$stress # 0.10
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

write.csv(group.scrs, "quad.group.scrs.csv")

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
par(mfrow=c(1,2), pty="s", mar=c(0.5, 0.5, 0.5, 0.5), oma=c(0,0,2,0))

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

legend("bottomright", legend=leg.txt, cex=1, pch=c(24,22,25,23,21), col=c(2,1,3,5,4), bty='n', pt.lwd=2)

text(.38, -0.2, '2D Stress = 0.10', cex=1)

mtext("Quadrat (sessile and small mobile taxa)", side=3, line=0.5)

###############################################################
# transect data

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

# CENTROIDS
### get principal coordinates; use classic (metric) MDS scaling
### what k should i use??? maximum = 8, because 8 species
mds1 <- cmdscale(spp.dist, k=8)
mds1
summary(mds1)
#plot(mds1, col=c(env$region))

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
#plot(mds.pc)

### load env data at site level

site.env <- read.csv("site.dat2_edit_120712.csv", header=TRUE)
names(site.env)
site.env <- site.env[,1:12]
#mds.pc$points
#plot(mds.pc$points, col=c(site.env$region))
#ordisurf(mds.pc, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2)
#summary(ordisurf(mds.pc, site.env$sedim.per, add=TRUE, col="gray40", labcex=1, lty=2))

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

write.csv(group.scrs, "tran.group.scrs.csv")

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
# panel B
plot(mds.pc$points, type="n", xlab="", ylab="", xaxt='n', yaxt='n')

ordisurf(mds.pc$points, group$sedim.per, add=TRUE, labcex=0.8, col="gray40", lty=2)

points(Channel, pch=22, col=1, cex=2.1, lwd=2, bg="white")
points(Haro, pch=24, col=2, cex=2.1, lwd=2, bg="white")
points(Georgia, pch=25, col=3, cex=2.1, lwd=2, bg="white")
points(Sound, pch=23, col=5, cex=2.1, lwd=2, bg="white")
points(Hood, pch=21, col=4, cex=2.1, lwd=2, bg="white")

text(mds.pc$points, labels=group$site.no, cex=0.7, font=2)

legend("bottomleft", "b", adj=c(2,1), bty='n', cex=1.5)

leg.txt <- c("Haro", "Channel", "Rosario", "Sound", "Hood")

#legend("bottomright", legend=leg.txt, cex=1, pch=c(24,22,25,23,21), col=c(2,1,3,5,4), bty='n', pt.lwd=2)

text(.25, .28, '2D Stress = 0.11', cex=1)
mtext("Transect (large mobile taxa)", side=3, line=0.5)

###############################################################




