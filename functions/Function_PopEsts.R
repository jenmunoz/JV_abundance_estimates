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
#All lookup tables and spatial layers will be uploaded to the GitHub repo, so this step can be skipped by general users.

# ------------------------------------------------------------
# 1. Helper function to check dependencies
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
# 2. Estimate populations function
# ------------------------------------------------------------
popEsts <- function(species, polys) {
  
  #check for required packages
  required_pkgs <- c("sf", "dplyr", "BAMexploreR", "ebirdst", "terra", "osfr", "purrr")
  check_packages(required_pkgs)
  
  #load source tables for distribution models and population estimates
  sdmSources <- read.csv("LookupData/sdm_speciesList.csv")
  peSources <- read.csv("LookupData/popEst_speciesList.csv")
  
  #create function to determine if conservation polygons overlap with strata polygons
  polyOverlap <- function(pol) {
    poly_prj <- sf::st_transform(pol, terra::crs(peStrat))
    test <- sf::st_intersects(poly_prj, sf::st_union(peStrat), sparse = FALSE)
    return(as.logical(test))
  }
  
  #function to estimate population size from density model
  #pol = conservation polygon, sdm = species raster, fact = factor to transform units into individuals/pixel
  popEst_DensityModel <- function(pol, sdm, fact = 1, dataSource) {
    poly_v <- terra::vect(pol) %>%
      terra::project(terra::crs(sdm))
    abd <- terra::crop(sdm, poly_v, snap = "near", mask = T)
    abundance_est <- round(terra::global(abd, fun = "sum", na.rm = TRUE) * fact, -1) %>%
      dplyr::rename(pop_est = sum) %>%
      dplyr::mutate(species = sp,
                    season = "breeding",
                    sdmSource = dataSource,
                    popEstSource = dataSource)
    return(abundance_est)
  }
  
  #loop through species and estimate population size for all available data sources
  results <- list()
  for (sp in species) {
    #check if any data sources for species
    sdm <- sdmSources %>%
      dplyr::filter(common_name == sp)
    
    if(nrow(sdm) == 0) {cat("No data available for", sp, "\n")} else {
      #identify sdms available
      cat('Estimating population size for "', sp, '" using the following data sources: \n', sep = "")
      
      #--------------------------------------
      #eBird workflow
      #--------------------------------------
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
          cat("\t Downloading eBird relative abundance surface for:", sp, "\n")
          # Download seasonal max abundance data at 3 km
          suppressMessages(
            try({ebirdst::ebirdst_download_status(sp, pattern = "abundance_seasonal_mean_3km", download_occurrence = FALSE, dry_run = FALSE, force = TRUE)}, silent = TRUE)
          )
        }
        
        #load eBird relative abundance rasters
        abd <- ebirdst::load_raster(
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
            pe <- read.csv("LookupData/pif.csv") %>%
              dplyr::filter(common_name == sp)
            
            #load strata for regional PIF estimates
            invisible(capture.output({
              peStrat <- sf::st_read(dsn = "LookupData/modelExtents.gpkg", layer = "pif_reg")
            }))
            
            #set population estimate source for export later
            psource <- "PIF regional"
          }
          
          if(peSp$fws_reg == "Yes") {
            #load regional USFWS estimates
            pe <- read.csv("LookupData/usfws.csv") %>%
              dplyr::filter(common_name == sp)
            
            #load strata for regional USFWS estimates
            invisible(capture.output({
              peStrat <- sf::st_read(dsn = "LookupData/modelExtents.gpkg", layer = "usfws")
            }))
            
            #set population estimate source for export later
            psource <- "USFWS"
          }
          
          polyTest <- purrr::map(polys, polyOverlap)
          
          #split polygons that overlap with regional strata and those that don't
          polys_reg <- polys[unlist(polyTest)]
          polys_canus <- polys[!unlist(polyTest)]
          
          #1. query regional strata that intersect with polygons. These strata will be used to crop eBird relative abundance surfaces
          #2. extract regional population estimates for those strata
          #3. use eBird relative abundace surface to estimate population size of conservation polygons
          if(length(polys_reg) > 0) {
            pop_est_breeding_reg <- lapply(polys_reg, FUN = function(pol) {
              #1. query regional strata that intersect
              poly_prj <- sf::st_transform(pol, terra::crs(peStrat))
              #determine intersecting strata and remove slivers
              suppressWarnings(ints <- sf::st_intersection(peStrat, poly_prj))
              min_area <- units::set_units(3000000, "m^2") #using 3e6 m^2 as minimum area
              ints <- ints[sf::st_area(ints) > min_area, ] 
              strata <- peStrat %>%
                dplyr::filter(stratum %in% ints$stratum)
              
              #2. extract population estimates
              pepoly <- pe %>%
                dplyr::filter(stratum %in% strata$stratum) %>%
                dplyr::pull(pop_est)
              if(length(pepoly) > 0) {
                pepoly <- sum(pepoly)
              } else {
                pepoly <- 0
              }
              
              #3. estimate population size
              #match projections, crop and mask eBird surface with pop est strata
              strata_v <- terra::vect(strata)
              strata_prj <- terra::project(strata_v, terra::crs(abd_breeding))
              abd_strata <- terra::crop(abd_breeding, strata_prj) %>%
                terra::mask(mask = strata_prj)
              
              total_strata <- terra::global(abd_strata, fun = "sum", na.rm = TRUE)
              prop_strata <- abd_strata / total_strata$sum
              
              poly_v <- terra::vect(poly_prj) %>%
                terra::project(terra::crs(prop_strata))
              
              prop_poly <- terra::crop(prop_strata, poly_v, snap = "near", mask = T)
              prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
              
              abundance_est <- round(prop_poly_sum * pepoly, -1)  # Compute absolute abundance
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
        use_canus <- exists("polys_canus") && peSp$acad == "Yes" && length(polys_canus) > 0
        condition <- (peSp$acad == "Yes" && peSp$pif_reg == "No" && peSp$fws_reg == "No") ||
          use_canus
        
        if(condition) {
          polys_tmp <- if(use_canus) {polys_canus} else {polys}
          
          #load ACAD population estimates
          pe <- read.csv("LookupData/acad.csv") %>%
            dplyr::filter(common_name == sp)
          
          #load polygon for Canada/USA. Using pif_reg and manipulating to have only 1 stratum (Can/US)
          invisible(capture.output({
            canus <- sf::st_read(dsn = "LookupData/modelExtents.gpkg", layer = "pif_reg") %>%
              dplyr::mutate(stratum = "canus") %>%
              dplyr::group_by(stratum) %>%
              dplyr::summarize(geom = sf::st_union(geom))
          }))
          
          #match projections, crop and mask eBird surface with Canada/USA polygon
          canus_v <- terra::vect(canus) #change to SpatVect
          canus_prj <- terra::project(canus_v, terra::crs(abd_breeding)) #project Canada/USA to match eBird before cropping
          abd_canus <- terra::crop(abd_breeding, canus_prj, snap = "near", mask = T) #crop and mask eBird to Canada/USA
          
          #estimate population size for polygons
          total_canus <- terra::global(abd_canus, fun = "sum", na.rm = TRUE)
          prop_canus <- abd_canus / total_canus$sum
          
          pop_est_breeding_canus <- lapply(polys_tmp, FUN = function(pol) {
            poly_v <- pol %>%
              sf::st_transform(crs = sf::st_crs(prop_canus)) %>%
              terra::vect()
            prop_poly <- terra::crop(prop_canus, poly_v, snap = "near", mask = T)
            prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
            
            abundance_est <- round(prop_poly_sum * pe$acad_uscan, -1)  # Compute absolute abundance
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
        rm(list = intersect(c("polys_reg", "polys_canus", "polys_tmp"), ls()))
        
        #Non-breeding season Start
        ###########################################################
        abd_nonbreed <- abd[[names(abd) != "breeding"]] 
        
        if(peSp$acad == "Yes") {
          #load ACAD population estimates
          pe <- read.csv("LookupData/acad.csv") %>%
            dplyr::filter(common_name == sp)
          
          #estimate population size for polygons
          total_global <- terra::global(abd_nonbreed, fun = "sum", na.rm = TRUE)
          prop_global <- abd_nonbreed / total_global$sum
          
          pop_est_nonbreed <- lapply(polys, FUN = function(pol) {
            poly_v <- pol %>%
              sf::st_transform(crs = sf::st_crs(prop_global)) %>%
              terra::vect()
            prop_poly <- terra::crop(prop_global, poly_v) %>%
              terra::mask(mask = poly_v)
            prop_poly_sum <- terra::global(prop_poly, fun = "sum", na.rm = TRUE)
            
            abundance_est <- round(prop_poly_sum * pe$acad_global, -1)  # Compute absolute abundance
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
          dplyr::rename(pop_est = sum) %>%
          dplyr::mutate(sdmSource = "eBird",
                 species = sp) %>%
          dplyr::select(species, polyID, sdmSource, popEstSource, season, pop_est) %>%
          dplyr::arrange(polyID)
        rownames(ebirdResults) <- NULL
        
      } else {
        ebirdResults <- data.frame(species = NA,
                                   polyID = NA,
                                   sdmSource = NA,
                                   popEstSource = NA,
                                   season = NA,
                                   pop_est = NA)
      }#END OF eBIRD WORKFLOW
      
      #--------------------------------------
      #BAM workflow
      #--------------------------------------
      #create blank database
      bamResults <- data.frame(species = NA,
                               polyID = NA,
                               sdmSource = NA,
                               popEstSource = NA,
                               season = NA,
                               pop_est = NA)
      
      if(sdm$BAMv4 == "Yes") {# | sdm$BAMv5 == "Yes") { #### Excluding V5 for now while there's a glitch in BAMexploreR
        #Determine which version is available for given species. Use V5 if available
        #ver <- if(sdm$BAMv5 == "Yes") {"v5"} else {"v4"}
        ver = "v4"
        
        #confirm that conservation polygons overlap with BAM AOI
        invisible(capture.output({
          peStrat <- sf::st_read(dsn = "LookupData/modelExtents.gpkg", layer = paste0("bam", ver))
        }))
        polyTest <- purrr::map(polys, polyOverlap)
        polys_tmp <- polys[unlist(polyTest)]
        
        if(length(polys_tmp) > 0) {
          cat("\t BAM density model \n")
          #download BAM raster for given species
          dir.create("data/spatial/bamRasters", showWarnings = FALSE)
          #extract 4-letter code from species table and download BAM raster layer
          spCode <- BAMexploreR::spp_tbl %>%
            dplyr::filter(commonName == sp) %>%
            dplyr::pull(speciesCode)
          
          #check if species raster is already downloaded, then download if needed
          if(!file.exists(file.path("data/spatial/bamRasters",paste0("pred-", spCode, "-CAN-Mean.tif")))) {
            cat("\t Downloading BAM density model for:", sp, "\n")
            suppressMessages({
              abd <- BAMexploreR::bam_get_layer(spCode, ver, "data/spatial/bamRasters")
            })
          } else {
            abd <- list()
            abd[[spCode]] <- terra::rast(file.path("data/spatial/bamRasters",paste0("pred-", spCode, "-CAN-Mean.tif")))
          }
          
          #estimate population size for each polygon using existing function
          pop_est_bam <- purrr::map(polys_tmp, popEst_DensityModel,
                                    sdm = abd[[spCode]],
                                    fact = 100,
                                    dataSource = paste0("BAM", ver))
          # pop_est_bam <- lapply(polys_tmp, function(pol) {
          #   poly_v <- terra::vect(pol)
          #   invisible(capture.output({
          #     abundance_est <- BAMexploreR::bam_pop_size(abd, crop_ext = poly_v) %>%
          #       dplyr::rename(pop_est = total_pop) %>%
          #       dplyr::mutate(species = sp,
          #              season = "breeding",
          #              sdmSource = paste0("BAM", ver),
          #              popEstSource = paste0("BAM", ver))
          #   }))
          # 
          #   return(abundance_est)
          # })
          names(pop_est_bam) <- names(polys_tmp)
          
          #combine BAM results
          bamResults <- dplyr::bind_rows(pop_est_bam, .id = "polyID") %>%
            dplyr::select(species, polyID,sdmSource, popEstSource, season, pop_est) %>%
            dplyr::arrange(polyID)
        }
      }#END OF BAM WORKFLOW 
      
      #--------------------------------------
      #CGAM Workflow
      #--------------------------------------
      #create empty data
      cgamResults <- data.frame(species = NA,
                                polyID = NA,
                                sdmSource = NA,
                                popEstSource = NA,
                                season = NA,
                                pop_est = NA)
      
      if(sdm$CGAMv1 == "Yes") {
        #confirm that conservation polygons overlap with BAM AOI
        invisible(capture.output({
          peStrat <- sf::st_read(dsn = "LookupData/modelExtents.gpkg", layer = "cgam")
        }))
        
        polyTest <- purrr::map(polys, polyOverlap)
        polys_tmp <- polys[unlist(polyTest)]
        
        if(length(polys_tmp) > 0) {
          cat("\t CGAM density model \n")
          
          #download CGAM raster for given species
          cgamDir <- "data/spatial/cgamRasters"
          dir.create(cgamDir, showWarnings = FALSE)
          spCode <- read.csv("LookupData/IBPSpeciesCodes.csv") %>%
            dplyr::filter(COMMONNAME == sp) %>%
            dplyr::pull(SPEC)
          
          file <- osfr::osf_retrieve_node("csugd") %>%
            osfr::osf_ls_files(pattern = spCode)
          file_path <- file.path(cgamDir, file$name)
          
          if(!file.exists(file_path)) {
            cat("\t Downloading CGAM density model for:", sp, "\n")
            osfr::osf_download(file, path = cgamDir)
          }
          abd <- terra::rast(file_path)
          
          #estimate population size
          pop_est_cgam <- purrr::map(polys_tmp, popEst_DensityModel,
                                     sdm = abd,
                                     fact = 1,
                                     dataSource = "CGAMv1")
          
          names(pop_est_cgam) <- names(polys_tmp)
          
          #combine CGAM results
          cgamResults <- dplyr::bind_rows(pop_est_cgam, .id = "polyID") %>%
            dplyr::select(species, polyID,sdmSource, popEstSource, season, pop_est) %>%
            dplyr::arrange(polyID)
        }
        
      } #END of CGAM Workflow
      
      
      #--------------------------------------
      #DUC Workflow
      #--------------------------------------
      #create empty data
      ducResults <- data.frame(species = NA,
                                polyID = NA,
                                sdmSource = NA,
                                popEstSource = NA,
                                season = NA,
                                pop_est = NA)
      if(sdm$DUC == "Yes") {
        invisible(capture.output({
          peStrat <- sf::st_read(dsn = "LookupData/modelExtents.gpkg", layer = "duc")
        }))
        
        polyTest <- purrr::map(polys, polyOverlap)
        polys_tmp <- polys[unlist(polyTest)]
        
        if(length(polys_tmp) > 0) {
          cat("\t DUC density model \n")
          
          #load DUC raster for given species
          spCode <- read.csv("LookupData/IBPSpeciesCodes.csv") %>%
            dplyr::filter(COMMONNAME == sp) %>%
            dplyr::pull(SPEC)
          
          ducDir <- "data/spatial/ducRasters"
          file <- list.files(ducDir, pattern = spCode)
          abd <- terra::rast(file.path(ducDir, file))
          
          #estimate population size
          pop_est_duc <- purrr::map(polys_tmp, popEst_DensityModel,
                                     sdm = abd,
                                     fact = 0.16, #DUC models have 0.4 x 0.4 km pixels = 0.16km^2
                                     dataSource = "DUC")
          
          names(pop_est_duc) <- names(polys_tmp)
          
          #combine CGAM results
          ducResults <- dplyr::bind_rows(pop_est_duc, .id = "polyID") %>%
            dplyr::select(species, polyID, sdmSource, popEstSource, season, pop_est) %>%
            dplyr::arrange(polyID)
        }
        
      } #END DUC Workflow
      
      #combine results across data sources
      results[[sp]] <- rbind(ebirdResults,
                             bamResults,
                             cgamResults,
                             ducResults) %>%
        dplyr::filter(rowSums(!is.na(.)) > 0) %>% #remove rows with all NAs
        dplyr::arrange(polyID, season, sdmSource)
      
      #clear environment and RAM before running next species
      keep <- c("results", "sdmSources", "peSources", "polyOverlap", "popEst_DensityModel")
      rm(list = setdiff(ls(), keep),
         envir = environment())
      invisible(capture.output({gc()}))
    }
    
    
  } #END OF SPECIES LOOP
  
  #combine results across species
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  return(results)
} #END POP EST FUNCTION














