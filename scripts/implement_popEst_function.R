#create test polygon and species lists
library(sf)
# library(dplyr)
# library(BAMexploreR)
# library(ebirdst)
# library(terra)
# library(osfr)
# library(purrr)


#Import layers for testing
st_layers("data/test/ATWR.gpkg")
sand <- st_read(dsn = "data/test/SAND.gpkg", layer = "SAND_boundary")
atwr <- st_read(dsn = "data/test/ATWR.gpkg", layer = "ATWR_boundary")
polys <- list("SAND" = sand, "ATWR" = atwr)

polysf <- st_read(dsn = "data/spatial/examplePolys2.shp")
#break polys into a list
polys <- split(polysf, seq_len(nrow(polysf)))
names(polys) <- polysf$NAME_E

#Jenny and Andrew's test polygon
polys <- list("Columbia Wetland" = st_read("data/spatial/TestPolys/Columbia Wetland Corridor area FINAL.shp"))
species <- list("Sora", "American Coot", "Pied-billed Grebe")
#test random list of species from the NAWCA priority species list
nawca_list <- read.csv("data/nawca_acad_species_match.csv", stringsAsFactors = FALSE)
species <- nawca_list$NAWCA_species[sample(1:nrow(nawca_list), 3)]

source("functions/Function_PopEsts.R")
test <- popEsts(species, polys)






#testing Jenny's function
JennyGrebe <- estimate_pop_conservation_area_breeding("Pied-billed Grebe", polys[[1]], ACAD_clean)
JennySora <- estimate_pop_conservation_area_breeding("Sora", polys[[1]], ACAD_clean)
JennyCoot <- estimate_pop_conservation_area_breeding("American Coot", polys[[1]], ACAD_clean)
jenny <- rbind(JennySora, JennyCoot, JennyGrebe)
testJennynonBreed <- estimate_pop_conservation_area_global("Pied-billed Grebe", polys[[1]], ACAD_clean)
