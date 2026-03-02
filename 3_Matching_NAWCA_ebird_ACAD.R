###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_ _###_###_###_###_###_###_###_###_ 
## checking for data availability for teh species of interest
## ##
## Objective: To check of the species of interest, for which ones we have ebird data 
## This includes making sure that the taxonomy is properly matched 
## 
## 1) Check the access code to the ebird data  
## 2) Check the ebird list of species
## 3) 
## 4) 
## i) 
## ii) 
## iii) 
##
## Updated and annotated by Jenny Munoz
## Last updated: February 2026
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###__###_###_###_###_###_###_###_###_ 

# ================================================================
# 0) SETUP
# ================================================================

# Install packages (run once only if needed)
install.packages(c(
  "tidyverse", "janitor", "glue", "fs", "png",
  "viridis", "scales", "fields", "readr",
  "rnaturalearth", "sf", "raster", "ebirdst",
  "rmapshaper", "terra"
))

# Load libraries
library(dplyr)
library(janitor)
library(rnaturalearth)
library(sf)
library(raster)
library(ebirdst)
library(rmapshaper)
library(terra)
library(ggplot2)

# ================================
# 1)  DATA 
# ================================

# The ebird data-species list 
ebirdst::ebirdst_runs

nawca_list<- read.csv("data/nawca_acad_species_match.csv", stringsAsFactors = FALSE) %>%
  as_tibble() %>% 
  mutate(common_name=NAWCA_species) %>% 
  filter(ACAD=="Yes")

unique(nawca_list$ACAD)

# ================================
# 2) See FOR WHICH SPECIES WE HAVE EBIRD DATA 
# ================================

ebirdst_runs_selected<-ebirdst_runs %>% 
  dplyr::select (species_code,scientific_name,common_name, breeding_start, breeding_end, breeding_quality)
#View(ebirdst_runs_selected) # the overall list of ebird species 


# list of Bc birds
nawca_species<-nawca_list %>% dplyr::select(common_name,global_or_american_predominantly, ACAD, PIF )

match_nawca_species_ebird <- nawca_species %>%
  left_join(ebirdst_runs_selected, by = "common_name")
View(match_nawca_species_ebird  )

# # for which species we dont have it
missing_species_nawca_ebird <- match_nawca_species_ebird %>%
  filter(is.na(species_code)) %>%
  dplyr::select(common_name, scientific_name) %>%
  mutate(ebird_data="no")


# --#---#-------------------------------------------

# This is the list of species for which we dont have the data 

#  1 Ancient Murrelet          NA              no        
# 2 Black-crowned Night-Heron NA              no        
# 3 Cassin's Auklet           NA              no        
#  4 Common Murre              NA              no        
#  5 Hudsonian Whimbrel        NA              no        
#  6 Masked Duck               NA              no        
#  7 Northern Goshawk          NA              no        
#  8 Red-necked Phalarope      NA              no        
#  9 White-rumped Sandpiper    NA              no        
# 10 Yellow Rail               NA              no        
# > 

# Check manually for taxonomy or other causes of mistmach

 View(ebirdst::ebirdst_runs)

#Ancient Murrelet # no data 
#Black-crowned Night-Heron # difference in name, i just remove the "-" from Night-Heron from the NAWCA list manually 
#Cassin's Auklet  # no data 
#Common Murre # no data 
#   Hudsonian Whimbrel           
#   Masked Duck    
# Northern Goshawk  
# Red-necked Phalarope 
#White-rumped Sandpiper
#Yellow Rail     

# ------------------------------------------------



