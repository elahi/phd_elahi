# Robin Elahi
# 25 Sept 2013
# use ISA to compare clusters 1 and 2, low retention vs high retention (ON is no longer in cluster 2 with Sound and Hood sites)

# Appendix S1.  Indicator.Value function to calculate DufrÍne & Legendreís (1997) Indicator
# Value for all species in two groups in R (R Development Core Team 2008).
# Text following a ë#í annotates the code.
# J.D. Bakker, January 11 2007

# Data must be in txt format, with treatment in first column (coded as 1 or 2)
# Must have only 2 groups, and must have >1 species (can use a dummy species)
# No spaces in column headings (use '.' instead)
# All data values >= 0 (no blank cells)
# Number of permutations set by numitr argument in function


### INDICATOR.VALUE FUNCTION
Indicator.Value <- function( sppmatrix, numitr) {

## Set initial values for variables
grps <- sppmatrix[,1]
fulldata <- sppmatrix[,2:dim(sppmatrix)[2]]
data <- fulldata[,colSums(fulldata)!=0]  # removes species absent from both groups
numspp <- ncol(data)
relabu <- c()
relfrq <- c()
indval <- c()
rand <- c()
pval <- c()

## Calculate mean and presence/absence in each group
grp.1 <- subset(data,grps==1)
count1 <- nrow(grp.1) # number of lines in group 1
count3=(grp.1>0) # TRUE if >0, else FALSE
x.1 <- apply(grp.1,2,mean) # mean abundance of species in group 1
grp.2 <- subset(data,grps==2)
count2 <- nrow(grp.2) # number of lines in group 2
count4=(grp.2>0) # TRUE if >0, else FALSE
x.2 <- apply(grp.2,2,mean) # mean abundance of species in group 2

## Calculate Actual Specificity
a.1 <- x.1/(x.1 + x.2)
a.2 <- x.2/(x.1 + x.2)
A <- c(a.1, a.2)  # vector of relative abundances
relabu <- cbind(relabu,A)
relabu.actual <- relabu[,1]

## Calculate Actual Fidelity
b.1 <- apply(count3, 2, sum)
b.2 <- apply(count4, 2, sum)
b1 <- b.1/count1
b2 <- b.2/count2
B <- c(b1,b2)
relfrq <- cbind(relfrq,B)
relfrq.actual <- relfrq[,1]

## Calculate Actual IV
IV <- A * B * 100
indval <- cbind(indval,IV)
indval.actual <- indval[,1]

## BEGIN PERMUTATIONS
# Calculate Permutational IVs
for(i in 1:numitr) {
cat("\n","Iteration ", i)

rand <- matrix(sample(grps))  	# Permutes group identities
data.perm <- cbind(rand,data)

# Calculate mean and presence/absence in each Permutation Group
grp.1.perm <- subset(data.perm,rand==1,select=2:dim(data.perm)[2])
# Don't need to recalculate count1 as it's unchanged from actual data
count3 = (grp.1.perm>0) # TRUE if >0, else FALSE
x.1.perm <- apply(grp.1.perm,2,mean) # mean abundance of species in group 1
grp.2.perm <- subset(data.perm,rand==2,select=2:dim(data.perm)[2])
# Don't need to recalculate count2 as it's unchanged from actual data
count4 = (grp.2.perm>0) # TRUE if >0, else FALSE
x.2.perm <- apply(grp.2.perm,2,mean) # mean abundance of species in group 2

# Calculate Permutation Specificity
a.1.perm <- x.1.perm/(x.1.perm + x.2.perm)
a.2.perm <- x.2.perm/(x.1.perm + x.2.perm)
A.perm <- c(a.1.perm,a.2.perm) # vector of relative abundances
relabu <- cbind(relabu,A.perm)

# Calculate Permutation Fidelity
b.1.perm <- apply(count3, 2, sum)
b.2.perm <- apply(count4, 2, sum)
b1.perm <- b.1.perm/count1
b2.perm <- b.2.perm/count2
B.perm <- c(b1.perm,b2.perm)
relfrq <- cbind(relfrq,B.perm)

# Calculate Permutation IV
IV.perm <- A.perm * B.perm * 100
indval <- cbind(indval,IV.perm)
}  # END OF PERMUTATION i

## Calculate p-value for each species in each group
indval.perm <- indval[,2:dim(indval)[2]]
count5 = (indval.perm >= indval.actual)
gt.IV <- apply(count5, 1, sum)
pval=(gt.IV+1)/(numitr+1)

## Output results
out <- cbind(rep(1:2,each=numspp), relabu.actual, relfrq.actual, indval.actual, pval)
colnames(out) <- c("Grp","A","B","IV","pval")
rownames(out) <- rep(colnames(data),2)
out

## Save result files to R working directory
write.csv(out, file="IV.results.csv", append=FALSE)  # Summary data (Species, Group, A, B, IV, p-value)
write.csv(indval, file="IV.indval.csv", append=FALSE)  # All IVs (actual and permutations)

} # End Indicator.Value function


### READ DATA FILE AND EXECUTE FUNCTION
dat <- as.matrix(read.table(file.choose(),header=TRUE))
dim(dat)
# input file must be a text file!!!!
# convert to presence absence
library(vegan)
spp.matrix <- decostand(dat[,2:133], "pa")
dim(spp.matrix)
# cbind cluster column(1) with spp.matrix
spp.df <- cbind(dat[,1],spp.matrix)
dim(dat)
summary(spp.df)
# this analysis uses presence/absence data
Indicator.Value(spp.df, numitr = 9999)


#References
#DufrÍne, M. & Legendre, P. (1997) Species assemblages and indicator species: the need for a
#  flexible asymmetrical approach. Ecological Monographs, 67, 345-366.
#R Development Core Team (2008) R: A Language and Environment for Statistical Computing.
#  R Foundation for Statistical Computing, Vienna, Austria.

