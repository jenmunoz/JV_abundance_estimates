###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
## Relative to absolute abundance
## ##
## Objective: Translate relative abundance to absolute abundance using ebird and estimates of continental popultaion from PIF
## This code does the following:
## 0) Download e-bird 3*3 km rasters for species of interest
## 
## The subsequent scripts do the other steps 
##i) 
## ii) 
## iii) 
##
## Updated and annotated by Jenny Munoz
## Last updated: February 26
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

# ================================
# 0) SETUP
# ================================

# Install required Libraries 
# To access ebird data 

# Data Manipulation
install.packages("tidyverse")
install.packages("janitor")
install.packages("glue") # String Manipulation
install.packages("fs") # File Operations
install.packages("png") # Image Handling
# Data Visualization
install.packages("viridis")
install.packages("scales")
install.packages("fields")
install.packages("readr") # Data Input/Output
# Geospatial Data
install.packages("rnaturalearth")
install.packages("sf")
install.packages("raster")
install.packages("ebirdst")
install.packages("rmapshaper")
install.packages("terra")


# Load necessary libraries
library(dplyr) # Essential for data manipulation with functions to filter, arrange, summarize, etc.
library(janitor) # Functions for simple data cleaning
library(glue) # Useful for data-driven string interpolations
library(fs) # A cross-platform interface for file system operations
library(png) # Allows reading and writing PNG images
library(viridis) # Provides color scales for improving data visualization accessibility
library(scales) # Graphical scales for mapping data to aesthetics in visualizations
library(fields) # Tools for spatial data
library(readr) # Fast and friendly way to read rectangular data like CSV files
library(rnaturalearth) # Provides map data for creating high-quality maps
library(sf) # Used for handling simple features to work with geographic data
library(raster) # For raster data manipulation and analysis
library(ebirdst) # Tools for accessing and analyzing eBird Status and Trends data
library(rmapshaper) # Simplifies shapes for data visualizations
library(terra) # For working with raster and vector data in spatial analysis


library(ebirdst) # ebird data 
library(ggplot2) # fpr plots


# ================================================================
# 0) SETUP & ACCESS KEY
# === === === == === === == == ==  === === === == === === == == == 
#  eBird S&T access key 
# An access key is required to download eBird Status & Trends data.
# 1) Request a key here or look at your key here  https://ebird.org/st/request
# 2) Save the key for this session with set_ebirdst_access_key().
#    

set_ebirdst_access_key("f6me7thr51ul")  # <- replace for local testing only
# Where am I running this from (useful for path debugging)?
getwd()
# ebirdst package version (useful for reproducibility)
ebirdst_version() # Uisng version year 2022

# ================================================================
# 1) CHOOSING THE DATA DIRECTORY
# === === === == === === == == ==  === === === == === === == == == 

# ebirdst_data_dir() resolves the download directory using:
# 1) EBIRDST_DATA_DIR env var if set, otherwise
# 2) tools::R_user_dir("ebirdst", which = "data")
ebirdst_data_dir()

# If you want to override the default for THIS SESSION ONLY: this will be the folder where data is downloaded 
# (Pick a fast local SSD or a managed project folder.)
ebirdst::ebirdst_data_dir()
Sys.setenv(EBIRDST_DATA_DIR = "C:/Users/jmunoz/Local_BirdsCanada/1_JV_science_coordinator_role_local/1_Projects/14_relativetoabsolute_abundance/data/species_layers")

#### CHECk for which species we have data, extract the breeding season dates and teh quality of teh data rating == 2 and 3 consider acceptable

ebirdst_runs$breeding_start
ebirdst_runs$breeding_end
ebirdst_runs$breeding_quality

# ================================
# 2) GETTING THE EBIRD DATA -DOWNLOAD DATA FOR ONE SPECIES *** Practice
# === === === == === === == == ==  === === === == === === == == == 
# The function below checks what exists (dry_run = TRUE) before downloading.
# pattern = "full-year_max_3km" targets the 3-km full-year max product.
# Set download_occurrence=TRUE if you also want occurrence rasters.
# force = TRUE avoids prompts and overwrites if needed.

# DRY RUN: list what would be downloaded for Cinnamon Teal
#ebirdst_download_status( "cinnamon teal",pattern = "full-year_max_3km",download_occurrence = TRUE,dry_run = TRUE)

ebirdst_download_status( "cinnamon teal",pattern = c("_3km","seasonal"),download_abundance= TRUE, dry_run = TRUE)
ebirdst_download_status( "Cinnamon teal",download_regional= TRUE, pattern = "regional", dry_run = TRUE)

ebirdst_download()

stats<-load_regional_stats("Cinnamon Teal")

names(stats)

# Making sure the species modeled population is just North America ( if not the percentage of population might be for more continnent than America)
#total_pop_percent: the proportion of the seasonal modeled population within the region.
#continent_pop_percent: the proportion of the seasonal modeled population for the continent within the region. The continent_name column identifies the continent that the region falls within. Note that Yellow-bellied Sapsucker only occurs in North American so the total and continental proportions are identical.

regional <- load_regional_stats("Cinnamon teal")





# ========STEP ONE-=======================
#  ESTIMATING ON OR NORTH AMERICAN POPULATION PER PIXEL PER SEASONS

# ================================
# 3a)PROPORTION OR NORTH AMERICAN POPULATION PER PIXEL [BREEDING SEASON-ALL SPECIES]
# === === === == === === == == == 
# Here we are instead of using the global population we are calculating the percentage of population for each pixel in USA and Canada
# In other words we are cropping the ebird data only to NorthAMerica as we assume this species will only breed there and we have teh estimates specific for those areas
# seasonal relative abundance
abd_seasonal <- load_raster("cintea", 
                            product = "abundance",
                            period = "seasonal",
                            resolution="3km")

# load country polygon, union into a single polygon, and project

northamerica<- ne_countries(country = c("United States of America", 
                                        "Canada" )) %>% #"Mexico"
  st_union() %>% 
  st_transform(crs = st_crs(abd_seasonal)) %>% 
  # vect converts an sf object to terra format for mask()
  vect()


plot(northamerica)
# mask seasonal abundance ( cookie cutter of the abundance map to Nort America)
abd_seasonal_northamerica <- mask(abd_seasonal, northamerica)

# total north american relative abundance for each season
abd_noram_total <- global(abd_seasonal_northamerica , fun = "sum", na.rm = TRUE)

# proportion of north american population
prop_pop_northamerica  <- abd_seasonal_northamerica  / abd_noram_total$sum

plot(prop_pop_noram)

# ================================
# 3b)PROPORTION OR  AMERICAN POPULATION PER PIXEL [MIGRATORY AND NON-BREEDING SEASON- SPECIES WITH DISTRIBUTIONS ONLY IN CONTINENTAL AMERICA]
# === === === == === === == == == 
# Here we are instead of using the global population we are calculating the percentage of population for each pixel in Continental America
# seasonal relative abundance
abd_seasonal <- load_raster("cintea", 
                            product = "abundance",
                            period = "seasonal",
                            resolution="3km")

# load country polygon, union into a single polygon, and project
#america<- rnaturalearth::ne_countries(continent = "North America" )
#plot(america$geometry) # only works before transforming to a vector

america<- rnaturalearth::ne_countries(continent = c("North America", "South America" )) %>% 
  st_union() %>% 
  st_transform(crs = st_crs(abd_seasonal)) %>% 
  # vect converts an sf object to terra format for mask()
  vect()

plot(america)
# mask seasonal abundance ( cookie cutter of the abundance map to Nort America)
abd_seasonal_america <- mask(abd_seasonal, america)

# total north american relative abundance for each season
abd_america_total <- global(abd_seasonal_america , fun = "sum", na.rm = TRUE)

# proportion of north american population
prop_pop_america  <- abd_seasonal_america  / abd_america_total$sum

plot(prop_pop_america)


# ================================
# 3c)PROPORTION OR GLOBAL POPULATION PER PIXEL [MIGRATORY AND NON-BREEDING SEASON- SPECIES WITH DISTRIBUTIONS THAT EXPANDS OUTSIDE CONTINENTAL AMERICA]
# === === === == === === == == == 
# Here we are instead of using the global population we are calculating the percentage of population for each pixel globally
# Consider removing areas with no ebird coverage (Rusia )
# seasonal relative abundance
abd_seasonal <- load_raster("cintea", 
                            product = "abundance",
                            period = "seasonal",
                            resolution="3km")

# load country polygon, union into a single polygon, and project # Consider removing areas with no ebird coverage (Rusia )

# mask seasonal abundance ( cookie cutter of the abundance map to Nort America)

# total north american relative abundance for each season
abd_global_total <- global(abd_seasonal , fun = "sum", na.rm = TRUE)
# proportion of global population
prop_pop_global <- abd_seasonal / abd_global_total$sum
plot(prop_pop_global )

# ================================
# 4) GEOGRAPHIC CONSERVATION AREA OF INTEREST POLYGON
# === === === == === === == == == 
# Here you can bring the polygon of interest ()
# check teh projections 
crs(prop_pop_northamerica )
crs(abd_seasonal )
crs(bc_boundary)
# Boundary of interest 
conservation_area <- sf::st_read("data/conservation_polygon/BC_boundary_layer.shp") # vector file 

# transform to match projection of raster data ( the raster data uses teh ebird projection)
conservation_area_proj<-conservation_area %>% 
  st_transform(crs = st_crs(prop_pop_northamerica)) %>% 
  vect() # # vect converts an sf object to terra format for mask()

crs(bc_boundary_proj)

# ================================
# 4b) GEOGRAPHIC CONSERVATION AREA OF INTEREST POLYGON
# === === === == === === == == == 
# Here you can bring the polygon of interest ()
# check teh projections 
crs(prop_pop_northamerica )
crs(abd_seasonal )
crs(bc_boundary)
# Boundary of interest 
conservation_area <- sf::st_read("data/conservation_polygon/ColumbiaWetland/Columbia Wetland Corridor area FINAL.shp") # vector file 
plot(conservation_area)

# transform to match projection of raster data ( the raster data uses teh ebird projection)
conservation_area_proj<-conservation_area %>% 
  st_transform(crs = st_crs(prop_pop_northamerica)) %>% 
  vect() # # vect converts an sf object to terra format for mask()

crs(bc_boundary_proj)
# ================================
# 5) CROP THE LAYERS OF PROPORTION OF POPULATION TO CONSERVATION AREA
# === === === == === === == == == 
# Note that there are 3 scales 
# NorthAmerica
# mask seasonal abundance ( cookie cutter of the abundance map to the conservation Area)
prop_pop_northamerica_conservation_area <- mask(prop_pop_northamerica, conservation_area_proj)
# total percentage of population for each season for the conservation area
prop_pop_northamerica_conservation_area_seasonal <- global(prop_pop_northamerica_conservation_area , fun = "sum", na.rm = TRUE)
print (prop_pop_northamerica_conservation_area )

# Continental America ( This piece might be unnecesary!!!! )
# mask seasonal abundance ( cookie cutter of the abundance map to the conservation Area)
prop_pop_america_conservation_area <- mask(prop_pop_america, conservation_area_proj)
# total percentage of population for each season for the conservation area
prop_pop_america_conservation_area_seasonal <- global(prop_pop_america_conservation_area , fun = "sum", na.rm = TRUE)
print (prop_pop_america_conservation_area )


# Global 
# mask seasonal abundance ( cookie cutter of the abundance map to the conservation Area)
prop_pop_global_conservation_area <- mask(prop_pop_global, conservation_area_proj)
# total percentage of population for each season for the conservation area
prop_pop_global_conservation_area_seasonal <- global(prop_pop_global_conservation_area , fun = "sum", na.rm = TRUE)
print (prop_pop_global_conservation_area )

# For species that have a distribution only in Americas, the global and America proportion of population is the same. 
# =========STEPTWO-=======================
# === === === == === === == == == 
# 8) EXTRACT THE ACAD-PIF ESTIMATE OF SPECIES POPULATIONS
# === === === == === === == == == 
# INCLUDE THE GLOBAL AND NORTHAMERICA POPULATIONS ( CANADA AND USA)
# If downloading the ACAD file form "https://pif.birdconservancy.org/avian-conservation-assessment-database-scores/",
# Make sure that you change the number of the columns to something easy to read for R, for instance need to rename manually in teh csv file the column Pop.Size_US.Ca_num and Global.Pop.Size_num
options(scipen = 999)

ACAD_estimates<-read.csv("data/ACAD_global_2024_05_23.csv")
names(ACAD_estimates)

ACAD_estimates_selected <- ACAD_estimates %>%
  mutate(common_name = Common.Name) %>% #   mutate(common_name = gsub("\\s+", "_", tolower(Common.Name))) %>% 
  mutate(pop_size_us_ca_numb = parse_number(Pop.Size_US.Ca_num)) %>% 
  mutate(global_pop_size_num=parse_number(Global.Pop.Size_num))%>% 
  dplyr::select(common_name, pop_size_us_ca_numb,global_pop_size_num)

# ================================
# 9) ESTIMATE THE ABSOLUTE ABUNDANCE IN THE AREA OF CONSERVATION INTEREST [BREEDING SEASON]
# === === === == === === == == == 
# prop_pop_northamerica_conservation_area 
# prop_pop_america_conservation_area 
# prop_pop_global_conservation_area

ACAD_estimates_northamerica_cintea <- ACAD_estimates_selected %>%
  dplyr::filter(common_name =="Cinnamon Teal") %>%
  dplyr::select(pop_size_us_ca_numb) 

#Multiply the estimated proportion of population the conservation area, by the the estimates numbers of breeding birds in Canada and teh USA
est_abundance_conservation_area_seasonal<-prop_pop_northamerica_conservation_area_seasonal*(ACAD_estimates_northamerica_cintea$pop_size_us_ca_numb)

# Select breeding only
est_abundance_conservation_area_breeding<-est_abundance_conservation_area_seasonal["breeding", ] # select breeding only 

# ================================
# 10) ESTIMATE THE ABSOLUTE ABUNDANCE IN THE AREA OF CONSERVATION INTEREST [NON-BREEDING AND MIGRATORY SEASON - SPECIES WITH DISTRIBUTIONS ONLY IN CONTINENTAL AMERICA]
# === === === == === === == == == 
# USE THE GLOBAL POPULAION ESTIMATE FROM ACAD-PIF
# USE THE TOTAL RELATIVE ABUNDANCE SUMMED ACROSS CONTINENTAL AMERICA ( )
# prop_pop_northamerica_conservation_area 
# prop_pop_america_conservation_area 
# prop_pop_global_conservation_area

ACAD_estimates_global_cintea <- ACAD_estimates_selected %>%
  dplyr::filter(common_name =="Cinnamon Teal") %>%
  dplyr::select(global_pop_size_num) 

#Multiply the estimated proportion of population the conservation area, by the the estimates numbers of breeding birds Continental America
est_abundance_conservation_area_american_sp<-prop_pop_america_conservation_area_seasonal*(ACAD_estimates_global_cintea$global_pop_size_num) 

# Select all the other seasons ( except breeding)

est_abundance_conservation_area_seasons_other_american_sp <-
  est_abundance_conservation_area_american_sp[ c("nonbreeding","prebreeding_migration","postbreeding_migration"), ,  drop = FALSE ]

# ================================
# 11) ESTIMATE THE ABSOLUTE ABUNDANCE IN THE AREA OF CONSERVATION INTEREST [NON-BREEDING AND MIGRATORY SEASON - SPECIES WITH DISTRIBUTIONS OUTSIDE CONTINENTAL AMERICA]
# === === === == === === == == == 
# USE THE GLOBAL POPULAION ESTIMATE FROM ACAD-PIF
# USE THE TOTAL RELATIVE ABUNDANCE SUMMED ACROSS CONTINENTAL AMERICA ( )

# Use teh proportion of gloal population per pixel
# prop_pop_northamerica_conservation_area 
# prop_pop_america_conservation_area 
prop_pop_global_conservation_area

ACAD_estimates_global_cintea <- ACAD_estimates_selected %>%
  dplyr::filter(common_name =="Cinnamon Teal") %>%
  dplyr::select(global_pop_size_num) 

#Multiply the estimated proportion of population the conservation area, by the the estimates numbers of breeding birds globally 
est_abundance_conservation_area_global_sp<-prop_pop_global_conservation_area_seasonal*(ACAD_estimates_global_cintea$global_pop_size_num) 
est_abundance_conservation_area_seasons_other_global_sp <-est_abundance_conservation_area_global_sp[ c("nonbreeding","prebreeding_migration","postbreeding_migration"), ,  drop = FALSE ]

# the abundance from ACAD # breeding population estimates 
ACAD_estimates_northamerica_cintea
ACAD_estimates_america_cintea
ACAD_estimates_global_cintea

# the proportion of population 
prop_pop_northamerica_conservation_area_seasonal
prop_pop_america_conservation_area_seasonal
prop_pop_global_conservation_area_seasonal

# the estimated values using different scales, north america and america gives us different values, part of it is becase the estimates are different so the NAmerican only adds more precision
# would it be simpler to just use the global estimate of population size?
est_abundance_conservation_area_seasonal # this estimate shoudl be more precise for breeding population, but it is not appropiate for other seaosn as it assumes that  the species does migration only within North america,

est_abundance_conservation_area_american_sp # intead when using the global estimates  if include the otehr areas where it moves and the global estimate
est_abundance_conservation_area_global_sp


############################################################
# POPULATION SCALING FRAMEWORK
############################################################
# We use ACAD–PIF (Partners in Flight) as the trusted source of population estimates.
# IMPORTANT: PIF provides BREEDING population estimates only.

# Population estimates are available at two spatial scales:
#   1) Global
#   2) North America (USA + Canada)

############################################################
# BREEDING SEASON ASSUMPTION
############################################################
# For the breeding season:
# - We assume North American estimates are more precise (finer spatial scale).
# - Therefore, we use the North American total population  to scale breeding relative abundance.
# This is defensible because:
# - PIF estimates are explicitly breeding-season estimates for north america

############################################################
# NONBREEDING & MIGRATION SEASONS ISSUE
############################################################
# For nonbreeding and migration seasons:
# Using the North American total population would imply:
#   - The entire breeding population remains in North America.
#   - Individuals only redistribute within North America.
# This creates a biological inconsistency because:
#   - Many species migrate outside North America.
#   - The total NA breeding population would be forced into a smaller seasonal distribution area durig the otehr seasons 
#
# Consequence:
#   -> Pixel-level abundance becomes inflated.
#   -> Seasonal densities are artificially high.
#
# This explains why abundance values are larger when using
# North American totals compared to global totals.

############################################################
# SOLUTION FOR NONBREEDING SEASONS
############################################################
# For nonbreeding and migration seasons we Use GLOBAL population estimate for scaling.
#
# Rationale:
# - We dont have information for the rest of the seasons,
# - Total number of individuals remains constant globally during the year 
# - Individuals redistribute across the globe seasonally.
# - Only a fraction of the global population is present in North America in nonbreeding/migration seasons.
#
# This avoids inflating abundance within North America.
############################################################
# OPEN QUESTION / CONCEPTUAL CHECK
############################################################
# Are we decreasing the uncertainty by using the NAmerican estimates instead if the global total for the breeding season?
#
############################################################
# NOTES
# For Cinnamon teal, teh estimated number of individuals is as expected different for North America than for all America because there is a large population in south Amecica, 
# but it this is just estimating the breeding population what does it mean
# ASK Barry what does that mean, does it mean that there are some Cinnamon teal breeding outside Canada and USA. 
# ANSWER YES Cinnamon teal breeds in south america too. so when we are using north america only, we are estimating using a finer scale calculation, 

# What does it means that using different  scales we have different values? Does it mean
## 


