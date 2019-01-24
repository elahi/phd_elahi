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

# First calculate available space
ses <- ses %>% 
  mutate(space = bare_rock + encrusting_coralline_algae + encrusting_non_calcified_algae)

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

mobq <- read_csv("Elahi_mobile_quadrat/output/elahi_mobile_quad_data.csv") %>% 
  select(-X1)
names(mobq)
mobq

# Sum abundances of all Tonicella species
mobq <- mobq %>% 
  mutate(tonicella = TIN + Ton)

mobq2 <- mobq %>% 
  select(quadrat:removals, tonicella) %>% 
  filter(comm.date < "2008-08-01") %>% 
  filter(comm.date > "2006-12-01")

mobq2 %>% count(removals)

#### MOBILE DATA FROM TRANSECTS ####


