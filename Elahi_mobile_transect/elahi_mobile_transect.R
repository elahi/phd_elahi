#################################################
# Tidy up PhD data
# Mobile taxa on transects
# Author: Robin Elahi
# Date: 160503
# NEED TO FIX DATES - IN SEPTEMBER 2011 (?)
#################################################

rm(list=ls(all=TRUE))

##### LOAD LIBRARIES AND DATA #####

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

dat <- read_csv("Elahi_mobile_transect/data/mobiledata_160503.csv")
glimpse(dat)

# Examine data
range(dat$length) # note that 0.01 indicates that size data were not available, but the animal was observed

# Deal with species codes
sppL <- sort(unique(dat$code))

# Dor becomes dorid
# ForEM lumped with FORem
# Her and HER becomes hermit
# Lit becomes lithode
# MKE and MOP become Mop
# Oel becomes OEL
# Per becomes perch
# Scu becomes sculpin
# SHR becomes shrimp
# TCA becomes TCAn
# 

code2 <- dat$code

code2 <- gsub('Dor', 'dorid', code2)
code2 <- gsub('ForEM', 'FORem', code2)
code2 <- gsub('Her', 'hermit', code2)
code2 <- gsub('HER', 'hermit', code2)
code2 <- gsub('Lit', 'lithode', code2)
code2 <- gsub('MKE', 'Mop', code2)
code2 <- gsub('MOP', 'Mop', code2)
code2 <- gsub('Oel', 'OEL', code2)
code2 <- gsub('Per', 'perch', code2)
code2 <- gsub('Scu', 'sculpin', code2)
code2 <- gsub('SHR', 'shrimp', code2)
code2 <- gsub('TCA', 'TCAn', code2)
code2 <- gsub('TCAnt', 'TCAt', code2)
code2 <- gsub('Tfe', 'TFE', code2)

sort(unique(code2))
sort(unique(dat$code))

sppMetaData <- data.frame("code" = sort(unique(code2)), 
                       "species" = NA)

# write.csv(sppMetaData, 'output/sppMetaData.csv')

dat$code <- code2
head(dat)

# Check dates
with(dat, table(comm.date, transect))

write.csv(dat, 'output/elahi_mobile_tran_data.csv')
unique(dat$surveyor)

