#################################################
# Prepare raw data for biodiversity change database
# Author: Robin Elahi
# Date: 151215
#################################################
# Dataset: Robin Elahi
# Subtidal rock wall community structure at Shady Cove, WA 1969-1974
# Data from the analysis of photographs taken by Charles Birkeland

# Sessile taxa
# Top left corner of quadrats were assessed for percent cover
# These data published in Elahi et al. 2013 Marine Biology

rm(list=ls(all=TRUE))

##### LOAD LIBRARIES #####

library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(tidyr)

##### LOAD ORIGINAL RAW DATA #####
#Load in the master data set
dat <- read.csv("./data/be SC_130719.csv")

# categorical factors
names(dat)

# Remove functional groups
df1 <- select(dat, quadrat:treat)
mat1 <- select(dat, BASP:RUBB)

# Reassemble dataframe
dat1 <- cbind(df1, mat1)
names(dat1)

# Write this csv file
write.csv(dat1, "./output/birkeland_sessile_data.csv")
# I modified this file in two ways:
# 1. Added a 'removals' column - was the urchin-chiton removal experiment happening?
# 2. Added a 2nd header row with full species names

##### COLLAPSE SOME OF THE UNIDENTIFED SPECIES INTO 'OTHER' CATEGORIES #####
# Read in the modified data file
dat <- read.csv("./output/birkeland_sessile_data_modified.csv", skip = 1)

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

# Collapse Metridium columns
dat6 <- dat5 %>% mutate(Metridium_spp = 
                          MEFA + ANE1 + Metridium_spp) %>%
  select(-MEFA, -ANE1)

names(dat6)

datW <- dat6 %>% select(-X)
write.csv(datW, "./output/birkeland_sessile_data_final_wide.csv")

##### WIDE TO LONG DATA #####
names(datW)
str(datW)

datW$dateR <- as.Date(datW$comm.date, origin = "1900-01-01")

