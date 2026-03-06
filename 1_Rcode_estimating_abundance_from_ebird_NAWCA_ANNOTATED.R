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
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# POPULATION SCALING FRAMEWORK
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# We use ACAD–PIF (Partners in Flight) as the trusted source of population estimates.
# IMPORTANT: PIF provides BREEDING population estimates only.

# Population estimates are available at two spatial scales:
#   1) Global
#   2) North America (USA + Canada)

###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# BREEDING SEASON ASSUMPTION
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# For the breeding season:
# - We assume North American estimates are more precise (finer spatial scale).
# - Therefore, we use the North American total population  to scale breeding relative abundance.
# This is defensible because:
# - PIF estimates are explicitly breeding-season estimates for north america

###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# NONBREEDING & MIGRATION SEASONS ISSUE
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
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

###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# SOLUTION FOR NONBREEDING SEASONS
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# For nonbreeding and migration seasons we Use GLOBAL population estimate for scaling.
#
# Rationale:
# - We dont have information for the rest of the seasons,
# - Total number of individuals remains constant globally during the year 
# - Individuals redistribute across the globe seasonally.
# - Only a fraction of the global population is present in North America in nonbreeding/migration seasons.
#
# This avoids inflating abundance within North America.
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# OPEN QUESTION / CONCEPTUAL CHECK
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# Are we decreasing the uncertainty by using the NAmerican estimates instead if the global total for the breeding season?
#
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# NOTES
#

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
# ================================================================
# eBird Status & Trends access key
# An access key is required to download eBird Status & Trends data.
#
# 1) Request or view your key here:
#    https://ebird.org/st/request
#
# 2) Save the key for this session with set_ebirdst_access_key()

set_ebirdst_access_key("f6me7thr51ul")  # <- replace for local testing only

# Where am I running this from (useful for path debugging)?
getwd()

# ebirdst package version (useful for reproducibility)
ebirdst_version()   # Using version year 2022


# ================================================================
# 1) CHOOSING THE DATA DIRECTORY
# ================================================================
# ebirdst_data_dir() resolves the download directory using:
# 1) EBIRDST_DATA_DIR environment variable if set, otherwise
# 2) tools::R_user_dir("ebirdst", which = "data")

ebirdst_data_dir()

# If you want to override the default for THIS SESSION ONLY
# This will be the folder where data is downloaded.
# (Pick a fast local SSD or a managed project folder.)

ebirdst::ebirdst_data_dir()

Sys.setenv(
  EBIRDST_DATA_DIR = "C:/Users/jmunoz/Local_BirdsCanada/1_JV_science_coordinator_role_local/1_Projects/14_relativetoabsolute_abundance/data/species_layers"
)


# ================================================================
# 2) CHECK SPECIES DATA AVAILABLE
# ================================================================
# Check breeding season dates and data quality.
# Quality ratings of 2 and 3 are considered acceptable.

ebirdst_runs$breeding_start
ebirdst_runs$breeding_end
ebirdst_runs$breeding_quality


# ================================================================
# 3) DOWNLOAD EBIRD STATUS & TRENDS DATA
# ================================================================
# The function below checks what exists (dry_run = TRUE) before downloading.
#
# pattern = "full-year_max_3km" targets the 3-km full-year max product.
# download_occurrence = TRUE if occurrence rasters are also required.
# force = TRUE avoids prompts and overwrites if needed.

# ---- DRY RUN EXAMPLES ----

# List seasonal 3-km abundance rasters
ebirdst_download_status(
  "cinnamon teal",
  pattern = c("_3km", "seasonal"),
  download_abundance = TRUE,
  dry_run = TRUE
)

# List regional statistics
ebirdst_download_status(
  "Cinnamon teal",
  download_regional = TRUE,
  pattern = "regional",
  dry_run = TRUE
)

# ---- DOWNLOAD ----
ebirdst_download()

# Load regional statistics
stats <- load_regional_stats("Cinnamon Teal")

names(stats)

# ---------------------------------------------------------------
# Check modeled population coverage
# ---------------------------------------------------------------
# total_pop_percent: proportion of the seasonal modeled population
# within the region.
#
# continent_pop_percent: proportion of the seasonal modeled
# population for the continent within the region.
#
# continent_name identifies the continent the region falls within.
#
# Example:
# Yellow-bellied Sapsucker only occurs in North America,
# therefore total and continental proportions are identical.

regional <- load_regional_stats("Cinnamon teal")


# ================================================================
# STEP ONE
# ESTIMATING PROPORTION OF POPULATION PER PIXEL
# ================================================================


# ================================================================
# 3a) PROPORTION OF NORTH AMERICAN POPULATION PER PIXEL
# (BREEDING SEASON – ALL SPECIES)
# ================================================================
# Instead of using the global population, we calculate the percentage
# of population for each pixel within USA and Canada.

# Seasonal relative abundance
abd_seasonal <- load_raster(
  "cintea",
  product = "abundance",
  period = "seasonal",
  resolution = "3km"
)

# Load country polygons, union them, and project
northamerica <- ne_countries(
  country = c("United States of America", "Canada")
) %>%
  st_union() %>%
  st_transform(crs = st_crs(abd_seasonal)) %>%
  vect()

plot(northamerica)

# Mask seasonal abundance to North America
abd_seasonal_northamerica <- mask(abd_seasonal, northamerica)

# Total North American relative abundance for each season
abd_noram_total <- global(
  abd_seasonal_northamerica,
  fun = "sum",
  na.rm = TRUE
)

# Proportion of North American population
prop_pop_northamerica <- abd_seasonal_northamerica / abd_noram_total$sum

plot(prop_pop_noram)


# ================================================================
# 3b) PROPORTION OF AMERICAN POPULATION PER PIXEL
# (MIGRATORY & NON-BREEDING – SPECIES ONLY IN THE AMERICAS)
# ================================================================

abd_seasonal <- load_raster(
  "cintea",
  product = "abundance",
  period = "seasonal",
  resolution = "3km"
)

america <- rnaturalearth::ne_countries(
  continent = c("North America", "South America")
) %>%
  st_union() %>%
  st_transform(crs = st_crs(abd_seasonal)) %>%
  vect()

plot(america)

abd_seasonal_america <- mask(abd_seasonal, america)

abd_america_total <- global(
  abd_seasonal_america,
  fun = "sum",
  na.rm = TRUE
)

prop_pop_america <- abd_seasonal_america / abd_america_total$sum

plot(prop_pop_america)


# ================================================================
# 3c) PROPORTION OF GLOBAL POPULATION PER PIXEL
# (SPECIES DISTRIBUTED OUTSIDE THE AMERICAS)
# ================================================================

abd_seasonal <- load_raster(
  "cintea",
  product = "abundance",
  period = "seasonal",
  resolution = "3km"
)

# Total global relative abundance
abd_global_total <- global(
  abd_seasonal,
  fun = "sum",
  na.rm = TRUE
)

# Proportion of global population
prop_pop_global <- abd_seasonal / abd_global_total$sum

plot(prop_pop_global)


# ================================================================
# 4) GEOGRAPHIC CONSERVATION AREA OF INTEREST
# ================================================================
# Import polygon of conservation interest.

crs(prop_pop_northamerica)
crs(abd_seasonal)
crs(bc_boundary)

conservation_area <- sf::st_read(
  "data/conservation_polygon/BC_boundary_layer.shp"
)

# Transform to match raster CRS
conservation_area_proj <- conservation_area %>%
  st_transform(crs = st_crs(prop_pop_northamerica)) %>%
  vect()

crs(bc_boundary_proj)


# ================================================================
# 4b) ALTERNATIVE CONSERVATION AREA
# ================================================================

conservation_area <- sf::st_read(
  "data/conservation_polygon/ColumbiaWetland/Columbia Wetland Corridor area FINAL.shp"
)

plot(conservation_area)

conservation_area_proj <- conservation_area %>%
  st_transform(crs = st_crs(prop_pop_northamerica)) %>%
  vect()

crs(bc_boundary_proj)


# ================================================================
# 5) CROP POPULATION PROPORTION TO CONSERVATION AREA
# ================================================================

# ---- NORTH AMERICA ----
prop_pop_northamerica_conservation_area <-
  mask(prop_pop_northamerica, conservation_area_proj)

prop_pop_northamerica_conservation_area_seasonal <-
  global(prop_pop_northamerica_conservation_area, fun = "sum", na.rm = TRUE)

print(prop_pop_northamerica_conservation_area)


# ---- CONTINENTAL AMERICA ----
prop_pop_america_conservation_area <-
  mask(prop_pop_america, conservation_area_proj)

prop_pop_america_conservation_area_seasonal <-
  global(prop_pop_america_conservation_area, fun = "sum", na.rm = TRUE)

print(prop_pop_america_conservation_area)


# ---- GLOBAL ----
prop_pop_global_conservation_area <-
  mask(prop_pop_global, conservation_area_proj)

prop_pop_global_conservation_area_seasonal <-
  global(prop_pop_global_conservation_area, fun = "sum", na.rm = TRUE)

print(prop_pop_global_conservation_area)

# For species distributed only in the Americas,
# global and America proportions will be identical.


# ================================================================
# STEP TWO
# EXTRACT ACAD-PIF POPULATION ESTIMATES
# ================================================================

options(scipen = 999)

ACAD_estimates <- read.csv("data/ACAD_global_2024_05_23.csv")

names(ACAD_estimates)

ACAD_estimates_selected <- ACAD_estimates %>%
  mutate(common_name = Common.Name) %>%
  mutate(pop_size_us_ca_numb = parse_number(Pop.Size_US.Ca_num)) %>%
  mutate(global_pop_size_num = parse_number(Global.Pop.Size_num)) %>%
  dplyr::select(
    common_name,
    pop_size_us_ca_numb,
    global_pop_size_num
  )


# ================================================================
# 9) ABSOLUTE ABUNDANCE (BREEDING SEASON)
# ================================================================

ACAD_estimates_northamerica_cintea <- ACAD_estimates_selected %>%
  dplyr::filter(common_name == "Cinnamon Teal") %>%
  dplyr::select(pop_size_us_ca_numb)

# Multiply proportion by ACAD estimate
est_abundance_conservation_area_seasonal <-
  prop_pop_northamerica_conservation_area_seasonal *
  (ACAD_estimates_northamerica_cintea$pop_size_us_ca_numb)

# Select breeding season
est_abundance_conservation_area_breeding <-
  est_abundance_conservation_area_seasonal["breeding", ]


# ================================================================
# 10) ABSOLUTE ABUNDANCE
# NON-BREEDING & MIGRATION (SPECIES IN THE AMERICAS)
# ================================================================

ACAD_estimates_global_cintea <- ACAD_estimates_selected %>%
  dplyr::filter(common_name == "Cinnamon Teal") %>%
  dplyr::select(global_pop_size_num)

est_abundance_conservation_area_american_sp <-
  prop_pop_america_conservation_area_seasonal *
  (ACAD_estimates_global_cintea$global_pop_size_num)

est_abundance_conservation_area_seasons_other_american_sp <-
  est_abundance_conservation_area_american_sp[
    c("nonbreeding", "prebreeding_migration", "postbreeding_migration"),
    ,
    drop = FALSE
  ]


# ================================================================
# 11) ABSOLUTE ABUNDANCE
# NON-BREEDING & MIGRATION (GLOBAL DISTRIBUTION SPECIES)
# ================================================================

ACAD_estimates_global_cintea <- ACAD_estimates_selected %>%
  dplyr::filter(common_name == "Cinnamon Teal") %>%
  dplyr::select(global_pop_size_num)

est_abundance_conservation_area_global_sp <-
  prop_pop_global_conservation_area_seasonal *
  (ACAD_estimates_global_cintea$global_pop_size_num)

est_abundance_conservation_area_seasons_other_global_sp <-
  est_abundance_conservation_area_global_sp[
    c("nonbreeding", "prebreeding_migration", "postbreeding_migration"),
    ,
    drop = FALSE
  ]


# ================================================================
# SUMMARY OUTPUTS
# ================================================================

# ACAD population estimates
ACAD_estimates_northamerica_cintea
ACAD_estimates_america_cintea
ACAD_estimates_global_cintea

# Proportion of population in conservation area
prop_pop_northamerica_conservation_area_seasonal
prop_pop_america_conservation_area_seasonal
prop_pop_global_conservation_area_seasonal

# Estimated abundance using different spatial scales
est_abundance_conservation_area_seasonal
est_abundance_conservation_area_american_sp
est_abundance_conservation_area_global_sp
