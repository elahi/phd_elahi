#################################################
# Prepare raw data for biodiversity change database
# Author: Robin Elahi
# Date: 151215
#################################################
# Datasets: Robin Elahi
# Analyses of subtidal rock wall community structure

# Sessile taxa

rm(list=ls(all=TRUE))

##### LOAD LIBRARIES #####

library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(tidyr)

##### DATASET 1 #####
### TIME-SERIES OF PERMANENT QUADRATS
### IN SAN JUAN CHANNEL
###

##### LOAD ORIGINAL RAW DATA #####
#Load in the master data set
dat <- read.csv("./data/sjc_sessile_130513.csv")

# categorical factors
names(dat)

# Remove functional groups
df1 <- select(dat, quad:treat)
mat1 <- select(dat, BASP:RUBB)

# Reassemble dataframe
dat1 <- cbind(df1, mat1)
names(dat1)

# Write this csv file
write.csv(dat1, "./output/sjc_data.csv")
# I modified this file in two ways:
# 1. Added a 'removals' column - was the urchin-chiton removal experiment happening?
# 2. Added a 2nd header row with full species names

##### COLLAPSE SOME OF THE UNIDENTIFED SPECIES INTO 'OTHER' CATEGORIES #####
# Read in the modified data file
dat <- read.csv("./output/sjc_data_modified.csv", skip = 1)

names(dat)

# Collapse SPS3, SPS8, MSP1, and sponge_other
dat3 <- dat %>% mutate(sponge_other = SPS3 + SPS8 + MSP1 + sponge_other) %>%
  select(-SPS3, -SPS8, -MSP1)

# Collapse other clonal tunicates
dat4 <- dat3 %>% mutate(tunicate_clonal_other = DITR + EUPU + SYPA + 
                          APPO + TCM2 + tunicate_clonal_other) %>%
  select(-DITR, -EUPU, -SYPA, -APPO, -TCM2)

# Collapse other solitary tunicates
dat5 <- dat4 %>% mutate(tunicate_solitary_other = 
                          TSV1 + TSV2 + tunicate_solitary_other) %>%
  select(-TSV1, -TSV2)
  
names(dat5)

datW <- dat5 %>% select(-X)
write.csv(datW, "./output/sjc_sessile_data_final_wide.csv")

##### WIDE TO LONG DATA #####
names(datW)
str(datW)

datW$dateR <- as.Date(datW$comm.date, origin = "1900-01-01")


##### DATASET 2 #####
### URCHIN ADDITION EXPERIMENT
### IN SAN JUAN CHANNEL
###

##### LOAD ORIGINAL RAW DATA #####
#Load in the master data set
dat <- read.csv("./data/urch_add_sessile_master_6_111202_edit_pg395_spp_mat.csv")

names(dat)

# Remove functional groups
df1 <- select(dat, Photo:Time)
mat1 <- select(dat, BASP:RUBB)

# Reassemble dataframe
dat1 <- cbind(df1, mat1)
names(dat1)

# Write this csv file
write.csv(dat1, "./output/urch_add_data.csv")
# I modified this file in two ways:
# 1. Added a 'removals' column - was the urchin-chiton removal experiment happening?
# 2. Added a 2nd header row with full species names

##### COLLAPSE SOME OF THE UNIDENTIFED SPECIES INTO 'OTHER' CATEGORIES #####
# Read in the modified data file
dat <- read.csv("./output/urch_add_data_modified.csv", skip = 1)

names(dat)

# Collapse SPS3, SPS8, and sponge_other
dat2 <- dat %>% mutate(sponge_other = SPS3 + SPS8 + MSP1 + sponge_other) %>%
  select(-SPS3, -SPS8, -MSP1)
names(dat2)

# Collapse Metridium columns
dat3 <- dat2 %>% mutate(Metridium_spp = 
                          ANE1 + Metridium_spp) %>%
  select(-ANE1)
names(dat3)

# Collapse other clonal tunicates
dat4 <- dat3 %>% mutate(tunicate_clonal_other = DITR + EUPU + SYPA + 
                          APPO + TCM2 + tunicate_clonal_other) %>%
  select(-DITR, -EUPU, -SYPA, -APPO, -TCM2)
names(dat4)

# Collapse other solitary tunicates
dat5 <- dat4 %>% mutate(tunicate_solitary_other = 
                          TSV1 + TSV2 + tunicate_solitary_other) %>%
  select(-TSV1, -TSV2)

names(dat5)

# Final dataset
datW <- dat5 %>% select(-X, -EALB)
write.csv(datW, "./output/urch_add_sessile_data_final_wide.csv")

##### WIDE TO LONG DATA #####
names(datW)
str(datW)

datW$dateR <- as.Date(datW$comm.date, origin = "1900-01-01")


##### DATASET 3 #####
### ANALYSIS OF BIRKELAND PHOTOGRAPHS 
### AT SHADY COVE (1969-1974)
###

# See 'birkeland_sessileQ.R' script


##### COMPILE ALL THREE DATASETS #####
### INTO A SINGLE DATAFRAME
### 
###

##### LOAD THREE FINAL MODIFIED SPREADSHEETS #####
# 
birke <- read.csv("./output/sjc_master_compiled/birkeland_sessile_data_final_wide_modified.csv")

elahi <- read.csv("./output/sjc_master_compiled/sjc_sessile_data_final_wide_modified.csv")

urchAdd <- read.csv("./output/sjc_master_compiled/urch_add_sessile_data_final_wide_modified.csv")

sessile_comp <- rbind(birke, elahi, urchAdd)

head(sessile_comp)
sessile_comp$dateR <- as.Date(sessile_comp$comm.date, format = "%m/%d/%Y")

write.csv(sessile_comp, "./output/sjc_master_compiled/sessile_compiled_final_wide.csv")

library(ggplot2)
ggplot(data = sessile_comp, aes(dateR, Metandrocarpa_taylori, color = site)) + 
  geom_point(alpha = 0.5) + facet_wrap(~ dataset, scales = "free", nrow = 3) + 
  geom_smooth()

library(dplyr)
sessile_comp %>% filter(dataset == "birkeland_permanent_quads") %>% 
  select(comm.date, dateR)

crap <- read.csv("./output/sjc_master_compiled/sessile_compiled_final_wide.csv")
crap %>% select(real.date, comm.date, dateR) %>% head()

