#' Create Critter
#' 
#' Creates a critter object. If only a scientific name is provided
#' create_critter will try and look up relevant life history from
#' FishLife
#' 
#' Critical inputs are adult_movement, adult_movement_sigma, and resolution
#' 
#'
#' @param common_name 
#' @param scientific_name 
#' @param critter_type 
#' @param seasonal_habitat 
#' @param habitat_seasons 
#' @param rec_habitat 
#' @param seasons 
#' @param rec_form 
#' @param adult_movement 
#' @param adult_movement_sigma 
#' @param fished_depletion 
#' @param init_explt initial annual exploitation rate
#' @param explt_type f or fmsy
#' @param ... 

create_critter <- function(common_name = 'white seabass',
                           scientific_name = NA,
                           critter_type = "fish",
                           seasonal_habitat = list(),
                           habitat_seasons = list(),
                           rec_habitat = NA,
                           seasons = 1,
                           rec_form = 1,
                           adult_movement = 0, 
                           adult_movement_sigma = 2,
                           dd_movement = 0,
                           fished_depletion = 0.4,
                           init_explt = .1,
                           explt_type = "f",
                           burn_years = 50,
                           ...) {

  
  init_explt <- init_explt / seasons # convert to seasonal exploitation rate
  
  if (critter_type == "fish"){
    critter <-
      marlin::Fish$new(
        common_name = common_name,
        scientific_name = scientific_name,
        seasonal_habitat = seasonal_habitat,
        habitat_seasons = habitat_seasons,
        rec_habitat = rec_habitat,
        seasons = seasons,
        rec_form = rec_form,
        fished_depletion = fished_depletion,
        adult_movement = adult_movement,
        adult_movement_sigma = adult_movement_sigma,
        dd_movement = dd_movement,
        init_explt = init_explt,
        explt_type = explt_type,
        burn_years = burn_years
      )
    
  }

  return(critter)
}
