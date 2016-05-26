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

##### FIX SPECIES CODES #####

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

##### ADD INFO ON URCHIN REMOVALS #####

mobDat <- read_csv("Elahi_mobile_quadrat/output/elahi_mobile_quad_data.csv")
glimpse(mobDat)

mobDat2 <- mobDat %>% filter(era == "Modern.R") %>%
  select(transect, comm.date, urch, removals) 
head(mobDat2)

urchinInfo <- mobDat %>% filter(era == "Modern.R") %>%
  select(transect, urch) %>% 
  distinct(transect)

# Bind treatment info to transect data
dat2 <- left_join(dat, urchinInfo, by = "transect")

# Add info for urchin removal dates
mobDat %>% filter(era == "Modern.R") %>%
  select(comm.date, removals) %>% 
  distinct(comm.date) 

dat2 %>% filter(urch == "u_rem" & code == "SFR") %>% 
  group_by(site, transect, code, comm.date) %>%
  summarise(n = n()) %>% ungroup() %>% 
  group_by(comm.date) %>% 
  summarise(urchins = mean(n, na.rm = TRUE)) %>% 
  ggplot(., aes(comm.date, urchins)) + geom_line()

# Grazer removals began on 18 April 2009 and 
# continued every two weeks until 24 January 2010. 
# (PLOS1 paper)
sort(unique(dat2$comm.date))

dat3 <- dat2 %>% 
  mutate(comm.date = ymd(comm.date), 
         removals = ifelse(
    comm.date > ymd(20090417) & comm.date < ymd(20100124), 
    "yes", "no"))

with(dat3, table(removals))

write.csv(dat3, 'Elahi_mobile_transect/output/elahi_mobile_tran_data.csv')
unique(dat$surveyor)

