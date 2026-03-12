###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
### CHECKING DATA AVAILABILITY FOR SPECIES OF INTEREST
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
##
## Objective:
## Check for which species of interest we have (1) eBird Status & Trends data, (2) BAM
## density models, (3) CGAM density models, (4) DUC density models.
## This includes verifying that taxonomy between datasets is properly matched.
##
## Workflow:
## 1) Check the access code to the eBird data
## 2) Retrieve the eBird list of species available in ebirdst
## 3) Compare with the NAWCA species list
## 4) Identify missing species and investigate taxonomy mismatches
##
## Updated and annotated by: Jenny Munoz
## Last updated: February 2026
###
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

# ================================================================
# 0) SETUP
# ================================================================

# Install packages (run once only if needed)
packages <- c(
  "tidyverse", "janitor", "glue", "fs", "png",
  "viridis", "scales", "fields", "readr",
  "rnaturalearth", "sf", "raster", "ebirdst",
  "rmapshaper", "terra", "BAMexploreR"
)

# Identify which ones are not yet installed
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]

# Install only the missing ones
if(length(new_packages) >0) {
  install.packages(new_packages)
}

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
library(BAMexploreR)

# ================================================================
# 1) DATA
# ================================================================

# eBird Status & Trends species list
ebirdst::ebirdst_runs

# NAWCA species list matched with ACAD
nawca_list <- read.csv(
  "data/nawca_acad_species_match.csv",
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(common_name = NAWCA_species) %>%
  filter(ACAD == "Yes")

# Check unique ACAD values
unique(nawca_list$ACAD)

# ================================================================
# 2) CHECK FOR WHICH SPECIES WE HAVE EBIRD, BAM, CGAM, AND DUC MODELS FOR
# ================================================================

# Select key fields from eBird runs table
ebirdst_runs_selected <- ebirdst_runs %>%
  dplyr::select(
    species_code,
    scientific_name,
    common_name,
    breeding_start,
    breeding_end,
    breeding_quality
  )

#create simplified eBird table to combine with NAWCA species table below
eBirdspp <- data.frame(common_name = ebirdst_runs_selected$common_name,
                       eBird = "Yes")

#get list of species for BAM v4 and v5
bam <- data.frame(common_name = bam_spp_list("v4", "commonName"),
                     BAMv4 = "Yes") %>%
  left_join(data.frame(common_name = bam_spp_list("v5", "commonName"),
                       BAMv5 = "Yes")) %>%
  mutate(BAMv5 = ifelse(is.na(BAMv5), "No", BAMv5))

#create list of species for CGAM models
cgam <- data.frame(common_name = c("Baird's Sparrow", "Bobolink", "Brewer's Sparrow", "Chestnut-collared Longspur", "Clay-colored Sparrow", 
                                   "Eastern Kingbird", "Grasshopper Sparrow", "Horned Lark", "Lark Bunting", "Long-billed Curlew", "LeConte's Sparrow",
                                   "Loggerhead Shrike", "Marbled Godwit", "Savannah Sparrow", "Sedge Wren", "Sprague's Pipit", "Thick-billed Longspur", 
                                   "Upland Sandpiper", "Vesper Sparrow", "Western Kingbird", "Western Meadowlark", "Willet"),
                   CGAMv1 = "Yes")
duc <- data.frame(common_name = c("Green-winged Teal", "American Wigeon", "Blue-winged Teal", "Canvasback", "Gadwall", "Mallard", "Northern Pintail", "Redhead", "Lesser Scaup"),
                  DUC = "Yes")


# View full list if needed
# View(ebirdst_runs_selected)

# ================================================================
# 3) MATCH NAWCA SPECIES WITH DIFFERENT MODEL SOURCES
# ================================================================

# List of NAWCA species
nawca_species <- nawca_list %>%
  dplyr::select(
    common_name,
    global_or_american_predominantly,
    ACAD,
    PIF
  )

#join species list for each model type with nawca species list
nawca_species_models <- left_join(nawca_species, eBirdspp) %>%
  left_join(bam) %>%
  left_join(cgam) %>%
  left_join(duc) %>%
  mutate(across(everything(), ~coalesce(.x, "No")))

# Match NAWCA species to eBird data
match_nawca_species_ebird <- nawca_species %>%
  left_join(ebirdst_runs_selected)

View(match_nawca_species_ebird)

# ================================================================
# 4) IDENTIFY SPECIES WITH NO POPULATION ESTIMATE AND/OR MODEL SOURCE
# ================================================================

missing_species_models <- nawca_species_models %>%
  filter(if_all(all_of(c("eBird", "BAMv4", "BAMv5", "CGAMv1", "DUC")), ~ .x == "No")) %>%
  pull(common_name)

missing_species_models #only 3

missing_species_pop <- nawca_species_models %>%
  filter(if_all(all_of(c("ACAD", "PIF")), ~ .x == "No")) %>%
  pull(common_name)

missing_species_pop

missing_species_nawca_ebird <- match_nawca_species_ebird %>%
  filter(is.na(species_code)) %>%
  dplyr::select(common_name, scientific_name) %>%
  mutate(ebird_data = "no")

#export updated NAWCA species table with model sources


# ------------------------------------------------
# List of species without eBird Status & Trends data
# ------------------------------------------------
#
# 1  Ancient Murrelet
# 2  Black-crowned Night-Heron
# 3  Cassin's Auklet
# 4  Common Murre
# 5  Hudsonian Whimbrel
# 6  Masked Duck
# 7  Northern Goshawk
# 8  Red-necked Phalarope
# 9  White-rumped Sandpiper
# 10 Yellow Rail
#

# ================================================================
# 5) MANUAL TAXONOMY CHECK
# ================================================================

# Inspect full eBird species table
View(ebirdst::ebirdst_runs)

# Notes on taxonomy or data availability:
#
# Ancient Murrelet            -> No data available
# Black-crowned Night-Heron   -> Name mismatch (hyphen removed manually in NAWCA list)
# Cassin's Auklet             -> No data available
# Common Murre                -> No data available
# Hudsonian Whimbrel          -> Needs verification
# Masked Duck                 -> Needs verification
# Northern Goshawk            -> Needs verification
# Red-necked Phalarope        -> Needs verification
# White-rumped Sandpiper      -> Needs verification
# Yellow Rail                 -> Needs verification
#
# ------------------------------------------------
