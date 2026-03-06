###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
### CHECKING DATA AVAILABILITY FOR SPECIES OF INTEREST
###_###_####_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
##
## Objective:
## Check for which species of interest we have eBird Status & Trends data.
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
# 2) CHECK FOR WHICH SPECIES WE HAVE EBIRD DATA
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

# View full list if needed
# View(ebirdst_runs_selected)

# ================================================================
# 3) MATCH NAWCA SPECIES WITH EBIRD DATA
# ================================================================

# List of NAWCA species
nawca_species <- nawca_list %>%
  dplyr::select(
    common_name,
    global_or_american_predominantly,
    ACAD,
    PIF
  )

# Match NAWCA species to eBird data
match_nawca_species_ebird <- nawca_species %>%
  left_join(ebirdst_runs_selected, by = "common_name")

View(match_nawca_species_ebird)

# ================================================================
# 4) IDENTIFY SPECIES WITHOUT EBIRD DATA
# ================================================================

missing_species_nawca_ebird <- match_nawca_species_ebird %>%
  filter(is.na(species_code)) %>%
  dplyr::select(common_name, scientific_name) %>%
  mutate(ebird_data = "no")

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
