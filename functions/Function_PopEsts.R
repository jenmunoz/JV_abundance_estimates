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

species <- list("Savannah Sparrow", "White-throated Sparrow", "Mallard")

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
popEsts <- function(species, polys, breeding = TRUE, nonbreeding = TRUE) {
  #load source tables for distribution models and population estimates
  sdmSources <- read.csv("data/sdm_speciesList.csv")
  peSources <- read.csv("data/PopEsts/modified/popEst_speciesList.csv")
  
  #loop through species and estimate population size for all available data sources
  for (sp in species) {
    #identify sdms available
    sdm <- sdmSources %>%
      dplyr::filter(common_name == sp)
    
    #eBird population estimation workflow
    ################################################
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
      
      #load eBird relative abundance rasters
      abd <- load_raster(
        sp,
        product = "abundance",
        period = "seasonal",
        metric = "mean", # note that here we can use mean or max
        resolution = "3km")
      
      #Determine which large-scale population estimate to use. Prioritize regional estimates if possible (PIF or USFWS). Use ACAD global (non breeding) or US-Can (breeding) if regional are unavailable
      peSp <- peSources %>%
        dplyr::filter(common_name == sp)
      
      
      #Breeding Season Start
      ####################################################################
      if(breeding = T) {
        #extract breeding season raster
        abd_breeding <- abd[["breeding"]]
        
        #Regional Pop Ests Start
        #################################################
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
          
          #1. query regional strata that intersect with polygons. These strata will be used to crop eBird relative abundance surfaces
          #2. extract regional population estimates for those strata
          #3. use eBird relative abundace surface to estimate population size of conservation polygons
          pop_est_breeding <- lapply(polys, FUN = function(pol) {
            #1. query regional strata that intersect
            poly_prj <- sf::st_transform(pol, crs(peStrat))
            #determine intersecting strata and remove slivers
            suppressWarnings(ints <- sf::st_intersection(peStrat, poly_prj))
            min_area <- units::set_units(3000000, "m^2") #using 3e6 m^2 as minimum area
            ints <- ints[sf::st_area(ints) > min_area, ] 
            strata <- peStrat %>%
              dplyr::filter(stratum %in% ints$stratum)
            
            #2. extract population estimates
            pepoly <- pe %>%
              filter(stratum %in% strata$stratum) %>%
              pull(pop_est)
            if(pepoly > 1) {
              pepoly <- sum(pepoly)
            }
            
            #3. estimate population size
            #match projections, crop and mask eBird surface with pop est strata
            strata_v <- terra::vect(strata)
            strata_prj <- terra::project(strata_v, crs(abd_breeding))
            abd_strata <- terra::crop(abd_breeding, strata_prj) %>%
              terra::mask(mask = strata_prj)
            
            total_strata <- terra::global(abd_strata, fun = "sum", na.rm = TRUE)
            prop_strata <- abd_strata / total_strata$sum
            
            poly_v <- terra::vect(poly_prj) %>%
              project(crs(prop_strata))
            
            prop_poly <- terra::crop(prop_strata, poly_v) %>%
              terra::mask(mask = poly_v)
            prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
            
            abundance_est <- prop_poly_sum * pepoly  # Compute absolute abundance
            return(abundance_est)
          })
        }
        #Regional Pop Ests End
        ######################################################
        
        #ACAD Canada/US estimates Start
        ######################################################
        if(peSp$acad == "Yes" & peSp$pif_reg == "No" & peSp$fws_reg == "No") {
          #load ACAD population estimates
          pe <- read.csv("data/PopEsts/modified/acad.csv") %>%
            dplyr::filter(common_name == sp)
          
          #load polygon for Canada/USA. Using pif_reg and manipulating to have only 1 stratum (Can/US)
          canus <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "pif_reg") %>%
            dplyr::mutate(stratum = "canus") %>%
            dplyr::group_by(stratum) %>%
            dplyr::summarize(geom = sf::st_union(geom))
          
          #match projections, crop and mask eBird surface with Canada/USA polygon
          canus_v <- terra::vect(canus) #change to SpatVect
          canus_prj <- terra::project(canus_v, crs(abd_breeding)) #project Canada/USA to match eBird before cropping
          abd_canus <- terra::crop(abd_breeding, canus_prj) %>%
            mask(mask = canus_prj) #crop and mask eBird to Canada/USA
          
          #estimate population size for polygons
          total_canus <- global(abd_canus, fun = "sum", na.rm = TRUE)
          prop_canus <- abd_canus / total_canus$sum
          
          pop_est_breeding <- lapply(polys, FUN = function(pol) {
            poly_v <- pol %>%
              sf::st_transform(crs = sf::st_crs(prop_canus)) %>%
              terra::vect()
            prop_poly <- terra::crop(prop_canus, poly_v) %>%
              terra::mask(mask = poly_v)
            prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
            
            abundance_est <- prop_poly_sum * pe$acad_uscan  # Compute absolute abundance
            return(abundance_est)
          })
        }
        #ACAD Canada/US estimates End
        #######################################################
        
      }
      #Breeding season end
      ###########################################################
      
      #Non-breeding season Start
      ###########################################################
      if(nonbreeding = T) {
        abd_nonbreed <- abd[[c("nonbreeding", "prebreeding_migration", "postbreeding_migration")]]
        
        if(peSp$acad == "Yes") {
          #load ACAD population estimates
          pe <- read.csv("data/PopEsts/modified/acad.csv") %>%
            dplyr::filter(common_name == sp)
          
          #estimate population size for polygons
          total_global <- global(abd_nonbreed, fun = "sum", na.rm = TRUE)
          prop_global <- abd_nonbreed / total_global$sum
          
          pop_est_nonbreed <- lapply(polys, FUN = function(pol) {
            poly_v <- pol %>%
              sf::st_transform(crs = sf::st_crs(prop_global)) %>%
              terra::vect()
            prop_poly <- terra::crop(prop_global, poly_v) %>%
              terra::mask(mask = poly_v)
            prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
            
            abundance_est <- prop_poly_sum * pe$acad_global  # Compute absolute abundance
            return(abundance_est)
          })
        } else {stop("Species not found in ACAD database")}
      }
      
      
   

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
























