################################################################################
##' @title Process raw survey data
##' @author Robin Elahi
##' @contact elahi.robin@gmail.com
##' @date 2019-01-23
##' @log Add a log here
################################################################################

##' Reproduce the data used in 2012 Elahi & Sebens MEPS
##' Interested in quadrat level estimates of:
##' available space, sessile species richness, chiton density
##' Also transect level estimates of urchin density

#### LOAD PACKAGES, DATA ####

# Tidyverse
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(vegan)

#### SESSILE DATA FROM QUADRATS ####

# Sessile data
ses <- read_csv("Elahi_sessile_quadrat/output/sjc_sessile_data_final_wide.csv") %>% 
  select(-X1)
unique(ses$quad)

# First calculate available space
ses <- ses %>% 
  mutate(space = bare_rock + encrusting_coralline_algae + 
           encrusting_non_calcified_algae)

# Create species matrix
names(ses)
spec_mat <- ses %>% 
  select(red_foliose_algae:invert_other)

ses$richness <- specnumber(spec_mat)
ses$simpsondiv <- diversity(spec_mat, index = "simpson")

ses2 <- ses %>% select(quad:removals, space, richness, simpsondiv)

## Select appropriate sampling dates - July 2008
ses2 <- ses2 %>% filter(month == 0807)


#### MOBILE DATA FROM QUADRATS ####
## Need to fix quadrat SC5_50 SC5_50z

mobq <- read_csv("Elahi_mobile_quadrat/output/elahi_mobile_quad_data.csv") %>% 
  select(-X1)

unique(mobq$quadrat)
mobq$quadrat <- gsub(pattern = "z", replacement = "",  
                     x = mobq$quadrat, fixed = TRUE)
names(mobq)

# Sum abundances of all Tonicella species
mobq <- mobq %>% 
  mutate(tonicella = TIN + Ton)

mobq2 <- mobq %>% 
  select(quadrat:removals, tonicella) %>% 
  filter(comm.date < "2008-08-01") %>% 
  filter(comm.date > "2006-12-01")

mobq2 %>% count(removals)

#### MOBILE DATA FROM TRANSECTS ####

mobt <- read_csv("Elahi_mobile_transect/output/elahi_mobile_tran_data.csv") 
names(mobt)
mobt %>% count(transect)

# Filter to relevant time period
mobt2 <- mobt %>% 
  filter(comm.date < "2008-08-01") %>% 
  select(id:length, comm.date)

## Get a list of all the transects
mobt_transects <- mobt2 %>% count(comm.date, site, transect) %>% select(-n)

## Summarise urchin abundances per transect
mobt3 <- mobt2 %>% 
  filter(code == "SFR") %>%
  count(comm.date, site, transect)
  
mobt3 %>% count(transect)

## join with complete transects
mobt3 <- left_join(mobt_transects, mobt3, 
                   by = c("comm.date", "transect", "site")) %>% 
  mutate(n = replace_na(n, replace = 0))


#### LINK THE DATASETS ####

mobt3 <- mobt3 %>% 
  group_by(site, transect) %>% 
  summarise(urchin_mean = mean(n)) %>% 
  ungroup()

mobq3 <- mobq2 %>% 
  group_by(site, transect, quadrat) %>% 
  summarise(chiton_mean = mean(tonicella)) %>% 
  ungroup()

ses3 <- ses2 %>% 
  select(quad, site, space:richness) %>% 
  rename(quadrat = quad)

dat <- left_join(mobq3, ses3, by = c("quadrat", "site"))
dat <- left_join(dat, mobt3, by = c("site", "transect"))
write.csv(dat, "processed_data/meps_2012_survey_data.csv")

dat %>% select(chiton_mean:urchin_mean) %>% pairs()
