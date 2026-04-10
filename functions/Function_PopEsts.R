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
library(terra)
library(osfr)


st_layers("data/test/ATWR.gpkg")
sand <- st_read(dsn = "data/test/SAND.gpkg", layer = "SAND_boundary")
atwr <- st_read(dsn = "data/test/ATWR.gpkg", layer = "ATWR_boundary")
polys <- list("SAND" = sand, "ATWR" = atwr)

species <- list("Savannah Sparrow", "White-throated Sparrow", "Mallard")

sp <- species[[3]]
pol <- polys[[1]]
# ------------------------------------------------------------
# 1. List of required packages
# ------------------------------------------------------------
required_pkgs <- c(
  "sf",      # for spatial functions
  "dplyr",   # for data manipulation
  "BAMexploreR",    # for raster/vector handling
  "ebirdst",
  "terra",
  "osfr"
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
popEsts <- function(species, polys) {
  #load source tables for distribution models and population estimates
  sdmSources <- read.csv("data/sdm_speciesList.csv")
  peSources <- read.csv("data/PopEsts/modified/popEst_speciesList.csv")
  
  #loop through species and estimate population size for all available data sources
  results <- list()
  for (sp in species) {
    #identify sdms available
    cat('Estimating population size for "', sp, '" using the following data sources: \n', sep = "")
    sdm <- sdmSources %>%
      dplyr::filter(common_name == sp)
    
    #eBird workflow
    ################################################
    if(sdm$eBird == "Yes") {
      cat("\t eBird relative abundance surfaces \n")
      #create and set directory for eBird rasters
      ebirdFolder <- "data/spatial/eBirdRasters"
      dir.create(ebirdFolder, showWarnings = FALSE)
      wd <- getwd()
      Sys.setenv(EBIRDST_DATA_DIR = file.path(wd, ebirdFolder))
      
      #download eBird relative abundance surface if no already present
      spCode <- ebirdst::get_species(sp)
      if(!any(basename(list.dirs(ebirdFolder, recursive = TRUE)) == spCode)) {
        cat("Downloading eBird relative abundance surface for:", sp, "\n")
        # Download seasonal max abundance data at 3 km
        suppressMessages(
          try({ebirdst_download_status(sp, pattern = "abundance_seasonal_mean_3km", download_occurrence = FALSE, dry_run = FALSE, force = TRUE)}, silent = TRUE)
        )
      }

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
          invisible(capture.output({
            peStrat <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "pif_reg")
          }))
          
          #set population estimate source for export later
          psource <- "PIF regional"
        }
        
        if(peSp$fws_reg == "Yes") {
          #load regional USFWS estimates
          pe <- read.csv("data/PopEsts/modified/usfws.csv") %>%
            dplyr::filter(common_name == sp)
          
          #load strata for regional USFWS estimates
          invisible(capture.output({
            peStrat <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "usfws")
          }))
          
          #set population estimate source for export later
          psource <- "USFWS"
        }
        
        #determine if conservation polygons overlap with strata polygons
        polyTest <- lapply(polys, function (pol) {
          poly_prj <- sf::st_transform(pol, crs(peStrat))
          test <- st_within(poly_prj, st_union(peStrat), sparse = FALSE)
          return(as.logical(test))
        })
        
        #split polygons that overlap with regional strata and those that don't
        polys_reg <- polys[unlist(polyTest)]
        polys_canus <- polys[!unlist(polyTest)]
        
        #1. query regional strata that intersect with polygons. These strata will be used to crop eBird relative abundance surfaces
        #2. extract regional population estimates for those strata
        #3. use eBird relative abundace surface to estimate population size of conservation polygons
        if(length(polys_reg) > 0) {
          pop_est_breeding_reg <- lapply(polys_reg, FUN = function(pol) {
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
            
            prop_poly <- terra::crop(prop_strata, poly_v, snap = "near", mask = T)
            prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
            
            abundance_est <- round(prop_poly_sum * pepoly, -2)  # Compute absolute abundance
            abundance_est$season <- rownames(abundance_est)
            abundance_est$popEstSource <- psource
            return(abundance_est)
          })
          names(pop_est_breeding_reg) <- names(polys_reg)
        }
        
      }#Regional Pop Ests End
      ######################################################
      
      #ACAD Canada/US estimates Start
      ######################################################
      if((peSp$acad == "Yes" & peSp$pif_reg == "No" & peSp$fws_reg == "No") | 
         (peSp$acad == "Yes" & length(polys_canus) > 0)) {
        
        if(peSp$acad == "Yes" & length(polys_canus) > 0) {
          polys_tmp <- polys_canus
          }  else {polys_tmp <- polys}
        
        #load ACAD population estimates
        pe <- read.csv("data/PopEsts/modified/acad.csv") %>%
          dplyr::filter(common_name == sp)
        
        #load polygon for Canada/USA. Using pif_reg and manipulating to have only 1 stratum (Can/US)
        invisible(capture.output({
        canus <- sf::st_read(dsn = "data/spatial/modified/modelExtents.gpkg", layer = "pif_reg") %>%
          dplyr::mutate(stratum = "canus") %>%
          dplyr::group_by(stratum) %>%
          dplyr::summarize(geom = sf::st_union(geom))
        }))
        
        #match projections, crop and mask eBird surface with Canada/USA polygon
        canus_v <- terra::vect(canus) #change to SpatVect
        canus_prj <- terra::project(canus_v, crs(abd_breeding)) #project Canada/USA to match eBird before cropping
        abd_canus <- terra::crop(abd_breeding, canus_prj, snap = "near", mask = T) #crop and mask eBird to Canada/USA
        
        #estimate population size for polygons
        total_canus <- global(abd_canus, fun = "sum", na.rm = TRUE)
        prop_canus <- abd_canus / total_canus$sum
        
        pop_est_breeding_canus <- lapply(polys_tmp, FUN = function(pol) {
          poly_v <- pol %>%
            sf::st_transform(crs = sf::st_crs(prop_canus)) %>%
            terra::vect()
          prop_poly <- terra::crop(prop_canus, poly_v, snap = "near", mask = T)
          prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
          
          abundance_est <- round(prop_poly_sum * pe$acad_uscan, -2)  # Compute absolute abundance
          abundance_est$season <- rownames(abundance_est)
          abundance_est$popEstSource <- "ACAD Can/USA"
          return(abundance_est)
        })
        names(pop_est_breeding_canus) <- names(polys_tmp)
        
        
        if(exists("pop_est_breeding_reg")) {
          pop_est_breeding <- c(pop_est_breeding_reg,pop_est_breeding_canus)
        } else {pop_est_breeding <- pop_est_breeding_canus}
      }#ACAD Canada/US estimates End
      
      #combine results with regional and canus estimates if both exists for a species (likely only waterfowl when one polygon is within USFWS strata and one is without)
      if (exists("pop_est_breeding_reg") && exists("pop_est_breeding_canus")) {
        pop_est_breeding <- c(pop_est_breeding_reg, pop_est_breeding_canus)
      } else if (exists("pop_est_breeding_reg")) {
        pop_est_breeding <- pop_est_breeding_reg
      } else if (exists("pop_est_breeding_canus")) {
        pop_est_breeding <- pop_est_breeding_canus
      }

      #remove unneeded polygon lists to save memory
      rm(polys_reg, polys_canus, polys_tmp)
      
      #Non-breeding season Start
      ###########################################################
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
          
          abundance_est <- round(prop_poly_sum * pe$acad_global, -2)  # Compute absolute abundance
          abundance_est$season <- rownames(abundance_est)
          abundance_est$popEstSource <- "ACAD global"
          return(abundance_est)
        })
        names(pop_est_nonbreed) <- names(polys)
      }
      
      #combine results from eBird for breeding and non-breeding season
      ebirdResults <- dplyr::bind_rows(
        c(pop_est_breeding,pop_est_nonbreed),
        .id = "polyID") %>%
        rename(pop_est = sum) %>%
        mutate(sdmSource = "eBird",
               species = sp) %>%
        select(species, polyID,sdmSource, popEstSource, season, pop_est) %>%
        arrange(polyID)
      rownames(ebirdResults) <- NULL
      
    } else {
      ebirdResults <- data.frame(species = NA,
                                 polyID = NA,
                                 sdmSource = NA,
                                 popEstSource = NA,
                                 season = NA,
                                 pop_est = NA)
      }#END OF eBIRD WORKFLOW
    
    #BAM workflow
    if(sdm$BAMv4 == "Yes" | sdm$BAMv5 == "Yes") {
      
      #confirm that conservation polygons overlap with BAM AOI
      
      polyTest <- lapply(polys, function (pol) {
        poly_prj <- sf::st_transform(pol, crs(peStrat))
        test <- st_within(poly_prj, st_union(peStrat), sparse = FALSE)
        return(as.logical(test))
      })
      
      cat("\t BAM density model \n")
      #Determine which version is available for given species. Use V5 if available
      if(sdm$BAMv5 == "Yes") {
        ver <- "v5"
      } else {ver <- "v4"}
      
      #download BAM raster for given species
      dir.create("data/spatial/bamRasters", showWarnings = FALSE)
      #extract 4-letter code from species table and download BAM raster layer
      spCode <- spp_tbl %>%
        filter(commonName == sp) %>%
        pull(speciesCode)
      
      #check if species raster is already downloaded, then download if needed
      if(!file.exists(file.path("data/spatial/bamRasters",paste0("pred-", spCode, "-CAN-Mean.tif")))) {
        abd <- bam_get_layer(spCode, ver, "data/spatial/bamRasters")
      } else {
        abd <- list()
        abd[[spCode]] <- rast(file.path("data/spatial/bamRasters",paste0("pred-", spCode, "-CAN-Mean.tif")))
      }
      
      #estimate population size for each polygon using existing function
      pop_est_bam <- lapply(polys, function(pol) {
        poly_v <- terra::vect(pol)
        invisible(capture.output({
          abundance_est <- bam_pop_size(abd, crop_ext = poly_v) %>%
            rename(pop_est = total_pop) %>%
            mutate(species = sp,
                   season = "breeding",
                   sdmSource = paste0("BAM", ver),
                   popEstSource = paste0("BAM", ver))
        }))
        
        return(abundance_est)
      })
      names(pop_est_bam) <- names(polys)
      
      #combine BAM results
      bamResults <- dplyr::bind_rows(pop_est_bam, .id = "polyID") %>%
        select(species, polyID,sdmSource, popEstSource, season, pop_est) %>%
        arrange(polyID)
    } else {
      bamResults <- data.frame(species = NA,
                            polyID = NA,
                            sdmSource = NA,
                            popEstSource = NA,
                            season = NA,
                            pop_est = NA)
      }#END OF BAM WORKFLOW 
    
    #CGAM Workflow
    if(sdm$CGAMv1 == "Yes") {
      cat("\t CGAM density model \n")
      
      #download CGAM raster for given species
      cgamDir <- "data/spatial/cgamRasters"
      dir.create(cgamDir, showWarnings = FALSE)
      spCode <- read.csv("data/IBPSpeciesCodes.csv") %>%
        dplyr::filter(COMMONNAME == sp) %>%
        pull(SPEC)
      
      file <- osfr::osf_retrieve_node("csugd") %>%
        osfr::osf_ls_files(pattern = spCode)
      file_path <- file.path(cgamDir, file$name)
      
      if(!file.exists(file_path)) {
        osfr::osf_download(file, path = cgamDir)
      }
      abd <- terra::rast(file_path)
      
      #estimate population size
      pop_est_cgam <- lapply(polys, function(pol) {
        poly_v <- terra::vect(pol) %>%
          terra::project(crs(abd))
        abd <- crop(abd, poly_v, snap = "near", mask = T)
        abundance_est <- round(global(abd, fun = "sum", na.rm = TRUE), -2) %>%
          rename(pop_est = sum) %>%
          mutate(species = sp,
                 season = "breeding",
                 sdmSource = "CGAMv1",
                 popEstSource = "CGAMv1")
        return(abundance_est)
      })
      names(pop_est_cgam) <- names(polys)
      
      #combine CGAM results
      cgamResults <- dplyr::bind_rows(pop_est_cgam, .id = "polyID") %>%
        select(species, polyID,sdmSource, popEstSource, season, pop_est) %>%
        arrange(polyID)
    } else {
      cgamResults <- data.frame(species = NA,
                               polyID = NA,
                               sdmSource = NA,
                               popEstSource = NA,
                               season = NA,
                               pop_est = NA)
    } #END of CGAM Workflow
    
    #combine results across data sources
    results[[sp]] <- rbind(ebirdResults,
                           bamResults,
                           cgamResults) %>%
      filter(rowSums(!is.na(.)) > 0) %>% #remove rows with all NAs
      dplyr::arrange(polyID, season, sdmSource)
    #clear environment and RAM before running next species
    keep <- c("results", "sdmSources", "peSources")
    rm(list = setdiff(ls(), keep),
       envir = environment())
    invisible(capture.output({gc()}))
    
  } #END OF SPECIES LOOP
  #combine results across species
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  return(results)
} #END POP EST FUNCTION

test <- popEsts(species, polys)

final <- do.call(rbind, test)

rownames(final) <- NULL


sdmSources <- read.csv("data/sdm_speciesList.csv")
peSources <- read.csv("data/PopEsts/modified/popEst_speciesList.csv")
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















