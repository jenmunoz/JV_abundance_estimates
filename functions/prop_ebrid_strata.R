prop_ebird_strata <- function(pol, sdm) {
  
  #1. query regional strata that intersect
  poly_prj <- sf::st_transform(pol, terra::crs(peStrat))
  #determine intersecting strata and remove slivers
  suppressWarnings(ints <- sf::st_intersection(peStrat, poly_prj))
  min_area <- units::set_units(3000000, "m^2") #using 3e6 m^2 as minimum area
  ints <- ints[sf::st_area(ints) > min_area, ] 
  strata <- peStrat %>%
    dplyr::filter(stratum %in% ints$stratum)
  
  #2. extract population estimates for intersecting strata
  pepoly <- pe %>%
    dplyr::filter(stratum %in% strata$stratum) %>%
    dplyr::pull(pop_est)
  if(length(pepoly) > 0) {
    pepoly <- sum(pepoly)
  } else {
    pepoly <- 0
  }
  
  
  #3. Create proportional relative abudnace raster for each polygon, matched to the appropriate strata
  strata_v <- terra::vect(strata) %>%
    terra::project(terra::crs(sdm))
  abd_strata <- terra::crop(sdm, strata_v, mask = T)
  
  total_strata <- terra::global(abd_strata, fun = "sum", na.rm = TRUE)
  prop_strata <- abd_strata / total_strata$sum
  
  return(list("pe" = pepoly, "prop_strata" = prop_strata))
}


test<- purrr::map(polys_reg, prop_ebird_strata, sdm = abd_breeding)
