###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
## Relative to absolute abundance
## ##
## Objective: Translate relative abundance to absolute abundance using ebird and estimates of continental population from PIF
## This code does the following:
## 
## 1) Check the access code to the ebird data  
## 2) Download e-bird 3*3 km rasters for species of interest
## 3) Enter the input data, including the species name and conservation area of interest
## 4) Choose the relevant function for your estimates 
## i) Breeding seasons all species
## ii) Species distributed in the Americas: Migratory and non-breeding seasons
## iii) Species distributed globally: Migratory and non-breeding seasons
## Updated and annotated by Jenny Munoz
## Last updated: February 2026
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

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
library(glue)
library(fs)
library(png)
library(viridis)
library(scales)
library(fields)
library(readr)
library(rnaturalearth)
library(sf)
library(raster)
library(ebirdst)
library(rmapshaper)
library(terra)
library(ggplot2)

# ================================================================
# 1) EBIRD STATUS & TRENDS ACCESS
# ================================================================
# This step allow you to download data from ebird 

# Request key: https://ebird.org/st/request
set_ebirdst_access_key("f6me7thr51ul")  # replace locally

getwd()              # Working directory
ebirdst_version()    # Data version (e.g., 2022)
ebirdst_data_dir()   # List the directory where you will be saving the ebird rasters

# Set data download directory (session only), if you want to have it on a speciefic folder
Sys.setenv(
  EBIRDST_DATA_DIR =
    "C:/Users/jmunoz/Local_BirdsCanada/1_JV_science_coordinator_role_local/1_Projects/14_relativetoabsolute_abundance/JV_birds_abundance_estimates/data/species_layers")

ebirdst_data_dir()

# ================================================================
# 2) DOWNLOAD SPECIES DATA
# ================================================================

# Dry run example, check what is available for a given species 

ebirdst_download_status( "Pied-billed grebe", pattern = "_3km", download_occurrence = TRUE, dry_run = TRUE)

# Download example
ebirdst_download_status( "American Coot", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Bank Swallow", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Sora", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "American Bittern", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Eared Grebe", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Marsh Wren", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Pied-billed grebe", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Virginia Rail", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Belted Kingfisher", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Common nighthawk", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)
ebirdst_download_status( "Common Loon", pattern = "_3km",download_abundance = TRUE,dry_run = FALSE)

load_raster()

# ================================================================
# 3) INPUT DATA
# ================================================================

# Create a folder named data and place there:
# 1)Polygon of conservation interest shape file 
# 2)ACAD file, can be downloaded from here :  "https://pif.birdconservancy.org/avian-conservation-assessment-database-scores/", SAVE IT as a csv file 
# Make sure that you change the naming of the columns is readable to R

# Check that there is information on ebird about your species and use the name as presented there 
# View(ebirdst_runs) # Available modeled species

#species<-"Sora"
#species<-"American Bittern"
#species<-"Marsh Wren"
#species<-"Eared Grebe"

# READ YOUR CONSERVATION POLYGON
# Conservation polygon (sf object)
# conservation_polygon <- 
#   st_read("data/conservation_polygon/BC/BC_boundary_layer.shp")

conservation_polygon <-
  st_read("data/conservation_polygon/ColumbiaWetland/Columbia Wetland Corridor area FINAL.shp")

#plot(conservation_polygon )
# ACAD population data
#ACAD_raw <- read.csv("data/ACAD_global_2024_05_23.csv")
ACAD_raw <- read.csv("data/ACAD Global 2024.05.23.csv")

ACAD_clean <- ACAD_raw %>%
  mutate(
    common_name = Common.Name,
    pop_us_ca   = parse_number(Pop.Size_US.Ca.),
    pop_global  = parse_number(Global.Pop.Size.)
  ) %>%
  dplyr::select(common_name, pop_us_ca, pop_global)

# ================================================================
# 4) FUNCTIONS – RELATIVE → ABSOLUTE ESTIMATES
# ================================================================
# Run once
# ------------------------------------------------
# BREEDING SEASON (North America only)
# Uses PIF US+Canada estimate
# Should be used for all species 
# ------------------------------------------------
estimate_pop_conservation_area_breeding <- function(species,
                                      conservation_polygon,
                                      ACAD_clean) {
  abd <- load_raster(
    species,
    product = "abundance",
    period = "seasonal",
    resolution = "3km")
  
  northamerica <- ne_countries(
    country = c("United States of America", "Canada") ) %>%
    st_union() %>%
    st_transform(crs = st_crs(abd)) %>%
    vect()
  
  abd_na <- mask(abd, northamerica)
  total_na <- global(abd_na, fun = "sum", na.rm = TRUE)
  
  prop_na <- abd_na / total_na$sum
  
  conservation_area <- conservation_polygon %>%
    st_transform(crs = st_crs(prop_na)) %>%
    vect()
  
  prop_area <- mask(prop_na, conservation_area)
  prop_area_sum <- global(prop_area, fun = "sum", na.rm = TRUE)
  
  pop_size <- ACAD_clean %>%
    filter(common_name == species) %>%
    pull(pop_us_ca)
  
  if (length(pop_size) == 0)
    stop("Species not found in ACAD.")
  
  abundance_est <- prop_area_sum * pop_size  # Compute absolute abundance
  
  note<-"Important, use the breeding season estimate only"
  note1<-"Uses ACAD US-CAD estimates and ebird cropped to US-CAD"
  

# Return data frame
abundance_estimate_northamerica_scaled <- data.frame(abundance = abundance_est, 
  species = species,
  note=note,
  framework=note1
)

  return(abundance_estimate_northamerica_scaled)
}


# ------------------------------------------------
# NON-BREEDING and MIGRATION – Species restricted to Americas
# Uses GLOBAL population estimate
# ------------------------------------------------

estimate_pop_conservation_area_americas <- function(species,
                                 conservation_polygon,
                                 ACAD_clean) {
  abd <- load_raster(
    species,
    product = "abundance",
    period = "seasonal",
    resolution = "3km")
  
  america <- ne_countries(
    continent = c("North America", "South America")) %>%
    st_union() %>%
    st_transform(crs = st_crs(abd )) %>% 
    vect() 
    
  abd_am <- mask(abd, america)
  total_am <- global(abd_am, fun = "sum", na.rm = TRUE)
  
  prop_am <- abd_am / total_am$sum
  
  conservation_area <- conservation_polygon %>%
    st_transform(crs = st_crs(abd)) %>%
    vect()
  
  prop_area <- mask(prop_am, conservation_area)
  prop_area_sum <- global(prop_area, fun = "sum", na.rm = TRUE)
  
  pop_size <- ACAD_clean %>%
    filter(common_name == species) %>%
    pull(pop_global)
  
  if (length(pop_size) == 0)
    stop("Species not found in ACAD.")
  
  abundance_est <- prop_area_sum * pop_size  # Compute absolute abundance
  
  note<-"Use for species distributed in the Americas. Extract the non-breeding and migratory seasons estimates"
  note1<-"Uses ACAD global estimates and ebird cropped to the americas"
  
  # Return data frame
  abundance_estimate_americas_scaled <- data.frame(abundance = abundance_est, 
                                    species = species,
                                   framework=note1,
                                    note=note )
  
  return(abundance_estimate_americas_scaled)
  
}


# ------------------------------------------------
# NON-BREEDING – Globally distributed species
# Uses GLOBAL population estimate
# ------------------------------------------------
estimate_pop_conservation_area_global <- function(species,
                                conservation_polygon,
                                ACAD_clean) {
  abd <- load_raster(
    species,
    product = "abundance",
    period = "seasonal",
    resolution = "3km")
  
  total_global <- global(abd, fun = "sum", na.rm = TRUE)
  prop_global <- abd / total_global$sum
  
  conservation_area <- conservation_polygon %>%
    st_transform(crs = st_crs(abd)) %>%
    vect()
  
  prop_area <- mask(prop_global, conservation_area)
  prop_area_sum <- global(prop_area, fun = "sum", na.rm = TRUE)
  
  pop_size <- ACAD_clean %>%
    filter(common_name == species) %>%
    pull(pop_global)
  
  if (length(pop_size) == 0)
    stop("Species not found in ACAD.")
  
  note<-"Use these estimates for species distributed Globally. Use the non-breeding and migratory seasons estimates, for breeding season use the breeding function"
  note1<-"Uses ACAD global estimates and global ebird"
    
  abundance_est <- prop_area_sum * pop_size  # Compute absolute abundance
  
  # Return data frame
  abundance_estimate_global_scaled <- data.frame(abundance = abundance_est, 
                                    species = species,
                                    method=note1,
                                    note=note )
  
  return(abundance_estimate_global_scaled)
  
}

# ================================================================
# 5) EXAMPLE CALL
# ================================================================

estimate_pop_conservation_area_breeding( species,conservation_polygon,ACAD_clean)

estimate_pop_conservation_area_breeding("American Coot",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_breeding("Sora",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_breeding("Marsh Wren",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_breeding("Eared Grebe",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_breeding("Virginia Rail",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_breeding("Bank Swallow",conservation_polygon,ACAD_clean)


estimate_pop_conservation_area_americas(species,conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_americas("American Coot",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_americas("Marsh Wren",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_americas("Eared Grebe",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_americas("American Coot",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_americas("Belted Kingfisher",conservation_polygon,ACAD_clean)


estimate_pop_conservation_area_global( species,conservation_polygon, ACAD_clean)
estimate_pop_conservation_area_global( "American Coot", conservation_polygon, ACAD_clean)
estimate_pop_conservation_area_global("American Bittern",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_global( "Sora",conservation_polygon, ACAD_clean)
estimate_pop_conservation_area_global("Eared Grebe",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_global("Bank Swallow",conservation_polygon,ACAD_clean)
estimate_pop_conservation_area_global("Belted Kingfisher",conservation_polygon,ACAD_clean)

# Here are some problem that we will encounter

# Species that have have small numbers ans are very secretive
estimate_pop_conservation_area_breeding("American Bittern",conservation_polygon,ACAD_clean) #


# Sandbox -----------------------------------------------------------------
# Here is a rational, we are using data from breeding population so the breeding estimate is restricted to the more precise information we have which it is US and Canada 
# For the other seasons we will use the one that bets describe the species distribution, America or Global for computational process

conservation_polygon <- 
  st_read("data/conservation_polygon/BurnsBog/Burns Bog Ecological Conservancy Area.kmz")

#### Rain check
#Calculate the percentage of population in Northamerica based in ebird4
#Compare it with the percentage of population form ACAD = North american pop/ Global population 
# This helps us with the argument that this values are comparable 

# I want to know the sum (total proportion ) of population that leaves in NA 

species<-"American Coot"

species<-"Bank swallow"

species<-"Common nighthawk"

species<-"Common Loon"



  
  

abd_seasonal<- load_raster(
  species,
  product = "abundance",
  period = "seasonal",
  resolution = "3km")


abd_global_total <- global(
  abd_seasonal,
  fun = "sum",
  na.rm = TRUE
)

prop_pop_global <- abd_seasonal /
  abd_global_total$sum


northamerica <- ne_countries(
  country = c("United States of America", "Canada")
) %>%
  st_union() %>%
  st_transform(crs = st_crs(abd_seasonal)) %>%
  vect()

prop_northamerica <- mask(prop_pop_global, northamerica)

layer1<-prop_northamerica[[1]]
total <- global(layer1, fun = "sum", na.rm = TRUE)


# == = == === ==== === == = === == === == === = ==
#### The end
# == = == === ==== === == = === == === == === = ==

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ================================================================
# THE WHOLE CODE – STEP BY STEP
# ================================================================
# NO NEED TO RUN THIS PART OF THE CODE
# THIS IS IF YOU ARE CURIOS ON INDIVIDUAL DECISIONS TAKEN AND STEP BY STEP
# ================================================================
# 0) SETUP
# ================================================================

# Install required libraries (run once if needed)
install.packages(c(
  "tidyverse", "janitor", "glue", "fs", "png",
  "viridis", "scales", "fields", "readr",
  "rnaturalearth", "sf", "raster",
  "ebirdst", "rmapshaper", "terra"
))

# Load libraries
library(dplyr)
library(janitor)
library(glue)
library(fs)
library(png)
library(viridis)
library(scales)
library(fields)
library(readr)
library(rnaturalearth)
library(sf)
library(raster)
library(ebirdst)
library(rmapshaper)
library(terra)
library(ggplot2)


# ================================================================
# 1) EBIRD STATUS & TRENDS ACCESS
# ================================================================

# Request key: https://ebird.org/st/request
set_ebirdst_access_key("f6me7thr51ul")  # replace locally

getwd()              # Working directory
ebirdst_version()    # Model version (e.g., 2022)

# Set download directory (session only)
Sys.setenv(
  EBIRDST_DATA_DIR =
    "C:/Users/jmunoz/Local_BirdsCanada/1_JV_science_coordinator_role_local/1_Projects/14_relativetoabsolute_abundance/data/species_layers"
)

ebirdst::ebirdst_data_dir()
ebirdst_runs


# ================================================================
# 2) DOWNLOAD EBIRD DATA
# ================================================================

# Example download
ebirdst_download_status(
  "Bank Swallow",
  pattern = "_3km",
  download_abundance = TRUE,
  dry_run = FALSE
)


# ================================================================
# STEP 1 – PROPORTION OF POPULATION
# ================================================================

# ------------------------------------------------
# 3a) NORTH AMERICAN PROPORTION (Breeding season)
# ------------------------------------------------

abd_seasonal <- load_raster(
  "banswa",
  product = "abundance",
  period = "seasonal",
  resolution = "3km"
)

northamerica <- ne_countries(
  country = c("United States of America", "Canada")
) %>%
  st_union() %>%
  st_transform(crs = st_crs(abd_seasonal)) %>%
  vect()

abd_seasonal_northamerica <- mask(abd_seasonal, northamerica)

abd_noram_total <- global(
  abd_seasonal_northamerica,
  fun = "sum",
  na.rm = TRUE
)

prop_pop_northamerica <- abd_seasonal_northamerica /
  abd_noram_total$sum

plot(prop_pop_northamerica)


# ------------------------------------------------
# 3b) CONTINENTAL AMERICA PROPORTION
# (Migratory / Non-breeding species restricted to Americas)
# ------------------------------------------------

america <- rnaturalearth::ne_countries(
  continent = c("North America", "South America")
) %>%
  st_union() %>%
  st_transform(crs = st_crs(abd_seasonal)) %>%
  vect()

abd_seasonal_america <- mask(abd_seasonal, america)

abd_america_total <- global(
  abd_seasonal_america,
  fun = "sum",
  na.rm = TRUE
)

prop_pop_america <- abd_seasonal_america /
  abd_america_total$sum

plot(prop_pop_america)


# ------------------------------------------------
# 3c) GLOBAL PROPORTION
# (Species extending outside continental America)
# ------------------------------------------------

abd_global_total <- global(
  abd_seasonal,
  fun = "sum",
  na.rm = TRUE
)

prop_pop_global <- abd_seasonal /
  abd_global_total$sum

plot(prop_pop_global)




# ================================================================
# 4) CONSERVATION AREA (BC)
# ================================================================

# Check projections
crs(prop_pop_northamerica)
crs(abd_seasonal)

# Load boundary
conservation_area <- sf::st_read(
  "data/conservation_polygon/BC_boundary_layer.shp"
)

# Transform to raster CRS
conservation_area_proj <- conservation_area %>%
  st_transform(crs = st_crs(prop_pop_northamerica)) %>%
  vect()

crs(conservation_area_proj)


# ================================================================
# 5) CROP PROPORTION LAYERS TO CONSERVATION AREA
# ================================================================

# --- North America
prop_pop_northamerica_conservation_area <-
  mask(prop_pop_northamerica, conservation_area_proj)

prop_pop_northamerica_conservation_area_seasonal <-
  global(prop_pop_northamerica_conservation_area,
         fun = "sum", na.rm = TRUE)


# --- Continental America
prop_pop_america_conservation_area <-
  mask(prop_pop_america, conservation_area_proj)

prop_pop_america_conservation_area_seasonal <-
  global(prop_pop_america_conservation_area,
         fun = "sum", na.rm = TRUE)


# --- Global
prop_pop_global_conservation_area <-
  mask(prop_pop_global, conservation_area_proj)

prop_pop_global_conservation_area_seasonal <-
  global(prop_pop_global_conservation_area,
         fun = "sum", na.rm = TRUE)

print(prop_pop_global_conservation_area)

# Note:
# If species occur only in the Americas,
# global and continental America proportions will be identical.


# ================================================================
# STEP 2 – ACAD POPULATION ESTIMATES
# ================================================================

options(scipen = 999)

ACAD_estimates <- read.csv(
  "data/ACAD_global_2024_05_23.csv"
)

species <- "Bank Swallow"

ACAD_estimates_selected <- ACAD_estimates %>%
  mutate(common_name = Common.Name) %>%
  mutate(pop_size_us_ca_numb =
           parse_number(Pop.Size_US.Ca_num)) %>%
  mutate(global_pop_size_num =
           parse_number(Global.Pop.Size_num)) %>%
  dplyr::select(common_name,
                pop_size_us_ca_numb,
                global_pop_size_num)


ACAD_estimates_northamerica_species <-
  ACAD_estimates_selected %>%
  dplyr::filter(common_name == species) %>%
  dplyr::select(pop_size_us_ca_numb)


# ================================================================
# 7) ESTIMATE ABSOLUTE ABUNDANCE
# ================================================================

# ----------------------------
# BREEDING (North America)
# ----------------------------

est_abundance_conservation_area_seasonal <-
  prop_pop_northamerica_conservation_area_seasonal *
  (ACAD_estimates_northamerica_species$pop_size_us_ca_numb)

est_abundance_conservation_area_breeding <-
  est_abundance_conservation_area_seasonal["breeding", ]


# ----------------------------
# NON-BREEDING – Americas only
# ----------------------------

ACAD_estimates_america_species <-
  ACAD_estimates_selected %>%
  dplyr::filter(common_name == "Bank Swallow") %>%
  dplyr::select(global_pop_size_num)

est_abundance_conservation_area_american_sp <-
  prop_pop_america_conservation_area_seasonal *
  (ACAD_estimates_america_species$global_pop_size_num)

est_abundance_conservation_area_seasons_other_american_sp <-
  est_abundance_conservation_area_american_sp[
    c("nonbreeding",
      "prebreeding_migration",
      "postbreeding_migration"),
    ,
    drop = FALSE
  ]


# ----------------------------
# NON-BREEDING – Global species
# ----------------------------

ACAD_estimates_global_species <-
  ACAD_estimates_selected %>%
  dplyr::filter(common_name == "Bank Swallow") %>%
  dplyr::select(global_pop_size_num)

est_abundance_conservation_area_global_sp <-
  prop_pop_global_conservation_area_seasonal *
  (ACAD_estimates_global_species$global_pop_size_num)

est_abundance_conservation_area_seasons_other_global_sp <-
  est_abundance_conservation_area_global_sp[
    c("nonbreeding",
      "prebreeding_migration",
      "postbreeding_migration"),
    ,
    drop = FALSE
  ]


# ================================================================
# 8) RESULTS
# ================================================================

# ACAD population values
ACAD_estimates_northamerica_species
ACAD_estimates_global_species
ACAD_estimates_america_species


# --- BC example outputs

northamerica_bc_bankswallow <-
  prop_pop_northamerica_conservation_area_seasonal

america_bc_bankswallow <-
  prop_pop_america_conservation_area_seasonal

global_bc_bankswallow <-
  prop_pop_global_conservation_area_seasonal


northamerica_bc_bankswallow_est <-
  est_abundance_conservation_area_seasonal

america_bc_bankswallow_est <-
  est_abundance_conservation_area_american_sp

global_bc_bankswallow_est <-
  est_abundance_conservation_area_global_sp








