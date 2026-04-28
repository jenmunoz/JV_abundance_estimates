pop_est_ebird <- function(pol, pe_sdm, popEstSource) {
  #extract population estimate and sdm from list
  pe <- pe_sdm[["pe"]]
  prop_strata <- pe_sdm[["prop_strata"]]
  
  poly_v <- terra::vect(pol) %>%
    terra::project(terra::crs(prop_strata))
  
  poly_area <- terra::expanse(poly_v, unit = "km") #calculate area of polygon
  results <- terra::extract(prop_global, poly_v, weights = T) %>% #extract values of pixels that intersect with polygon
    dplyr::summarise(dplyr::across(-c(ID, weight),
                                   ~ sum(.x * weight, na.rm = T))) %>% #sum pixel values weighted by proportion of pixel that occurs within the polygon
    tidyr::pivot_longer(cols = everything(),
                        names_to = "season",
                        values_to = "propPolySum") %>%
    dplyr::mutate(pop_est = round(propPolySum * pe, -1)) %>% #multiply weighted, summed proportions by total population size to get population estimate
    dplyr::mutate(density_sqkm = round(pop_est/poly_area, 3), #divide by area to get mean density
                  popEstSource = popEstSource) %>%
    select(season, popEstSource, pop_est, density_sqkm)
  return(results)
}


test2 <- purrr::map2(polys_reg, test, pop_est_ebird, popEstSource = "USFWS")
pol = polys_reg[[1]]
prop_strata<- test[[1]]$prop_strata
pe = test[[1]]$pe
