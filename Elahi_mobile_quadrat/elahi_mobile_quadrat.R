#################################################
# Tidy up PhD data
# Mobile taxa in quadrats
# Date: 151218
#################################################
# Datasets: Robin Elahi
# Analyses of subtidal rock wall community structure

# Mobile taxa in quadrats

rm(list=ls(all=TRUE))

##### LOAD LIBRARIES #####

library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(tidyr)

##### DATASET 1 #####
### TIME-SERIES OF PERMANENT QUADRATS AND BIRKELAND DATA TO MATCH
### IN SAN JUAN CHANNEL
### (NOTE THAT THE ELAHI ET AL. 2013 PUB SIMPLIFIED A FEW CATEGORIES DIFFERENTLY)

# load data
dat <- read.csv("./data/mobq_master_130725.csv")
dat <- read.csv("./data/birke_elahi_mob_130725.csv")

# categorical factors
names(dat)
names(be)

# Remove functional groups
df1 <- select(dat, quadrat:depth)
mat1 <- select(dat, SFR:Gunnel)

# Add in a removals column
df1$removals <- with(df1, ifelse(time < 905 | time > 1006, 
                                 "no", "yes"))

# Reassemble dataframe
dat1 <- cbind(df1, mat1)
names(dat1)

##' NAMING RULES:
##' Tonicella insignis = TIN (full genus and species)
##' Tonicella species = Ton (genus only)
##' dorid species = dorid (family or functional group name)
##' 

# Collapse some columns
dat2 <- dat1 %>% 
  mutate(Gastropod = GA1 + GA2 + GA5 + Gas, 
         Tonicella_spp = Ton + Ton3) %>%
  select(-GA1, -GA2, -GA5, - Gas, -ANA, -AHU, 
         -Aca, -FVE, -FTR, -Fla, -Ton, -Ton3) %>%
  rename(Aca = Acanthodoris_spp, chiton = Chiton_UNID, 
         crab = Crab_UNID, Den = Dendronotus, 
         Dio = Diodora, dorid = Dorid, Fla = Flabellina_spp, 
         flatworm = Flatworm, gastropod = Gastropod, 
         gunnel = Gunnel, hermit = Her, limpet = Limpet_UNID, 
         trochid = Marcal, nudi = Nudi_UNID, ophioroid = Oph, 
         Pug = Pugettia, sculpin = Sculpin_UNID, shrimp = Shr, 
         Ton = Tonicella_spp, Cal = CAL, CPA = cuke
         ) %>% 
  select(-c(transect.no, treat, time))
names(dat2)

# Order the species columns alphabetically
df2 <- dat2 %>% select(quadrat:removals)
mat2 <- dat2 %>% select(SFR:Ton)

mat3 <- mat2[, order(colnames(mat2), decreasing = FALSE)]
df3 <- cbind(df2, mat3)
names(df3)

sp1 <- df3 %>% select(Aca:trochid) %>% names() # %>% sort()

# Write this csv file
write.csv(df3, "./output/elahi_mobile_quad_data.csv")

##### DATASET 2 #####
### URCHIN ADDITION EXPERIMENT
### IN SAN JUAN CHANNEL
###

urch <- read.csv("./data/urch_add_mobile master_160502.csv")

# Collapse some columns
urch2 <- urch %>% 
  mutate(Ton_spp = Ton + TLI + Ton2 + Ton3, 
         Mop_spp = MKE + Mop, 
         Marcal = MAR + GA3_MAR,
         Amphissa_spp = Amp + GA4_AMP, 
         Gastropod = GA1 + GA2 + GA5 + GAS, 
         Acanthodoris_spp = ACA, 
         Flabellina_spp = FVE + FTR, 
         Her = HER, 
         Shr = SHR) %>%
  select(-c(Ton, TLI, Ton2, Ton3, MKE, Mop, MAR, GA3_MAR, GA4_AMP, Amp, 
            GA1, GA2, GA5, GAS, ACA, FVE, FTR, HER, SHR, WNEM, eggs))
names(urch2)

sp2 <- urch2 %>% select(SFR:Shr) %>% names() %>% sort()
sp2

# Rename columns to follow rules
urch3 <- urch2 %>% 
  rename(Aca = Acanthodoris_spp, Amp = Amphissa_spp, 
         Cal = CAL, chiton = chiton_other, 
         dorid = Dorid, Fla = Flabellina_spp, 
         gastropod = Gastropod, hermit = Her, 
         trochid = Marcal, Mop = Mop_spp, Nuc = NUC, 
         ophioroid = Oph, sculpin = Sculpin, 
         shrimp = Shr, Ton = Ton_spp)

names(urch3)
sp3 <- urch3 %>% select(SFR:shrimp) %>% names() %>% sort()
sp3

sp3unique <- which(!sp3 %in% sp1)
sp1unique <- which(!sp1 %in% sp3)

sp3[sp3unique] # no unique values in sp3 (urchin addition experiment)
sp1[sp1unique] # 16 unique species in sp1 (permanent quads)

# Order the species columns alphabetically
names(urch3)
df <- urch3 %>% select(Photo:Time)
mat1 <- urch3 %>% select(SFR:shrimp)

mat2 <- mat1[, order(colnames(mat1), decreasing = FALSE)]
urch4 <- cbind(df, mat2)
names(urch4)

# Write this csv file
write.csv(urch4, "./output/urch_add_mobile_data.csv")


