############################################################################
#Function title: Estimate bird species popultion size for given area(s)
#Input list of species and list of polygons (shapefiles) for which population estimates are needed. Function will output population 
#estiamtes using all available data sources for each species and polygon.
#Written by: Barry Robinson, barry.robinson@ec.gc.ca
#Date: March 25, 2026
############################################################################

#arguments:
#species = list of species for which population estimates are needed
#polys = list of shapefiles (sf objects) within which population estimates are needed

#Note that script 1_CreateLookupInfo.R is run first. This downloads and creates the appropriate look up tables and spatial layers to determine which data sets are available for given species and locations

#create test polygon and species lists
library(sf)
library(dplyr)
library(BAMexploreR)
library(ebirdst)

st_layers("data/test/ATWR.gpkg")
sand <- st_read(dsn = "data/test/SAND.gpkg", layer = "SAND_boundary")
atwr <- st_read(dsn = "data/test/ATWR.gpkg", layer = "ATWR_boundary")
polys <- list(sand, atwr)

species <- list("Savannah Sparrow", "White-troated Sparrow")

sp <- species[[1]]
pol <- polys[[1]]
# ------------------------------------------------------------
# 1. List of required packages
# ------------------------------------------------------------
required_pkgs <- c(
  "sf",      # for spatial functions
  "dplyr",   # for data manipulation
  "BAMexploreR",    # for raster/vector handling
  "ebirdst"
)

# ------------------------------------------------------------
# 2. Helper function to check dependencies
# ------------------------------------------------------------
check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing) > 0) {
    stop(
      "The following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nInstall them with install.packages(c(", 
      paste0("'", missing, "'", collapse = ", "), "))"
    )
  }
}

# ------------------------------------------------------------
# 3. Estimate populations function
# ------------------------------------------------------------
popEsts <- function(species, polys, season = "breeding") {
  #load source tables for distribution models and population estimates
  sdmSources <- read.csv("data/sdm_speciesList.csv")
  peSources <- read.csv("data/PopEsts/modified/popEst_speciesList.csv")
  
  #loop through species and estimate population size for all available data sources
  for (sp in species) {
    #identify sdms available
    sdm <- sdmSources %>%
      dplyr::filter(common_name == sp)
    
    #eBird population estimation workflow
    if(sdm$eBird == "Yes") {
      #create and set directory for eBird rasters
      dir.create("data/spatial/eBirdRasters", showWarnings = FALSE)
      wd <- getwd()
      Sys.setenv(EBIRDST_DATA_DIR = file.path(wd, "data/spatial/eBirdRasters"))
      
      #download relative abundance surface for species
      cat("/n>>> Downloading eBird relative abundance surface:", sp, "/n")
      # Download seasonal max abundance data at 3 km
      try({
        ebirdst_download_status(sp,pattern = "abundance_seasonal_mean_3km",download_occurrence = FALSE, dry_run = FALSE, force = TRUE)
      }, silent = TRUE)
      
      #Determine which large-scale population estimate to use. Prioritize regional estimates if possible (PIF or USFWS). Use ACAD global or US-Can if regional are unavailable
      peSp <- peSources %>%
        dplyr::filter(common_name == sp)
      
    ####################################################################
    #Breeding Season
      if(season = "breeding") {
        #################################################
        #If regional population estimates are available
        
        if(peSp$pif_reg == "Yes" | peSp$fws_reg == "Yes") {
          if(peSp$pif_reg == "Yes") {
            #load regional PIF estimates
            pe <- read.csv("data/PopEsts/modified/pif.csv") %>%
              dplyr::filter(common_name == sp)
            
            #load strata for regional PIF estimates
            peStrat <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "pif_reg")
          }
          
          if(peSp$fws_reg == "Yes") {
            #load regional USFWS estimates
            pe <- read.csv("data/PopEsts/modified/usfws.csv") %>%
              dplyr::filter(common_name == sp)
            
            #load strata for regional USFWS estimates
            peStrat <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "usfws")
          }
          
          #query regional strata that overlap with polygons. These strata will be used to crop eBird relative abundance surfaces
          cropPolys <- lapply(polys, FUN = function(x) {
            #reproject
            polytmp <- sf::st_transform(x, crs(peStrat))
            #determine intersecting strata and remove slivers
            suppressWarnings(ints <- sf::st_intersection(peStrat, polytmp))
            min_area <- units::set_units(3000000, "m^2") #using 3e6 m^2 as minimum area
            ints <- ints[sf::st_area(ints) > min_area, ] 
            intspoly <- peStrat %>%
              dplyr::filter(stratum %in% ints$stratum)
            return(intspoly)
          })
          
          #extract population estimates for regional strata associated with each polygon
          
          ##COLUMN NAMES FOR POP ESTS ARE DIFFERENT ACROSS TABLES. UPDATE IN SCRIPT 3 TO BE THE SAME
          ##FIGURE OUT HOW TO MAKE THE BELOW STEP UNIVERSAL FOR GLOBAL AND CAN-US POP OPTION TO MINIMZIE AMOUNT OF CODE
          pop_size <- lapply(cropPolys, function(x) {
            petmp <- pe %>%
              filter(stratum %in% x$stratum) %>%
              pull(pop_est)
            if(petmp > 1) {
              petmp <- sum(petmp)
            }
            return(petmp)
          })
        }
        ######################################################
      }
    ########################################################################
      

    }
    
  }
  
  #load population estimate tables and associated polygons
  usfws <- read.csv("data/PopEsts/modified/usfws.csv")
  acad <- read.csv("data/PopEsts/modified/acad.csv")
  pif_reg <- read.csv("data/PopEsts/modified/pif.csv")
  pifB <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "pif_reg")
  usfwsB <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "usfws")
  
  #load shapefile boundaries for density models
  bamv4 <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "bamv4")
  bamv5 <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "bamv5")
  cgam <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "cgam")
}
























