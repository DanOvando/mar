#' R6 Class Representing a Fish
#'
#' @description
#' A fish object has all the life history data needed to model a finfish population
#'
#' @details
#' Fish like to swim

Fish <- R6::R6Class(
  "Fish",
  lock_objects = FALSE,
  public = list(
    #' @description
    #' fill in fish object
    #' @param common_name common name
    #' @param scientific_name
    #' @param linf
    #' @param vbk
    #' @param t0
    #' @param cv_len
    #' @param length_units
    #' @param min_age
    #' @param max_age
    #' @param weight_a
    #' @param weight_b
    #' @param weight_units
    #' @param length_50_mature
    #' @param length_95_mature
    #' @param delta_mature
    #' @param age_50_mature
    #' @param age_95_mature
    #' @param age_mature
    #' @param length_mature
    #' @param m
    #' @param steepness
    #' @param r0
    #' @param ssb0
    #' @param density_dependence_form
    #' @param adult_movement
    #' @param adult_movement_sigma
    #' @param larval_movement
    #' @param query_fishlife
    #' @param sigma_r
    #' @param rec_ac
    #' @param cores
    #' @param mat_mode
    #' @param default_wb
    #' @param tune_weight
    #' @param density_movement_modifier
    #' @param linf_buffer
    #' @param resolution
    #' @param seasonal_habitat
    #' @param habitat_seasons
    #' @param rec_habitat
    #' @param fished_depletion
    #' @param rec_form
    #' @param burn_years
    #' @param seasonal_hab
    #' @param seasons
    #' @param init_explt
    #' @param explt_type
    initialize = function(common_name = 'white seabass',
                          scientific_name = NA,
                          linf = NA,
                          vbk = NA,
                          t0 = -0.1,
                          cv_len = 0.1,
                          length_units = 'cm',
                          min_age = 0,
                          max_age = NA,
                          weight_a = NA,
                          weight_b = NA,
                          weight_units = 'kg',
                          length_50_mature = NA,
                          length_95_mature = NA,
                          delta_mature = .1,
                          age_50_mature = NA,
                          age_95_mature = NA,
                          age_mature = NA,
                          length_mature = NA,
                          m = NA,
                          steepness = 0.8,
                          r0 = 10000,
                          ssb0 = NA,
                          density_dependence_form = 1,
                          adult_movement = 2,
                          adult_movement_sigma = 2,
                          larval_movement = 2,
                          query_fishlife = T,
                          sigma_r = 0,
                          rec_ac = 0,
                          cores = 4,
                          mat_mode = "age",
                          default_wb = 2.8,
                          tune_weight = FALSE,
                          dd_movement = 0,
                          linf_buffer = 1.2,
                          resolution = 25,
                          seasonal_habitat = list(),
                          habitat_seasons = list(),
                          rec_habitat = NA,
                          fished_depletion = 1,
                          rec_form = 1,
                          burn_years = 50,
                          seasonal_hab = NA,
                          seasons = 1,
                          explt_type = "f",
                          init_explt = .1) {
      seasons <- as.integer(seasons)
      
      if (seasons < 1) {
        stop("seasons must be greater than or equal to 1")
      }
      #
      
      if (length(seasonal_habitat) == 0) {
        seasonal_habitat <-
          purrr::map(1:seq_along(seasons), function(x, res)
            matrix(1, nrow = res, ncol = res), res = resolution)
        
      }
      
      if (length(habitat_seasons) == 0) {
        
        seasons_per_habitat <- seasons / length(seasonal_habitat) # determine how many seasons to assign to each habitat block
       
        
        tmp <- data.frame(block =  rep(seq_along(seasonal_habitat), each = seasons_per_habitat), season =  1:seasons)
        
        habitat_seasons <- vector(mode = "list", length = length(seasonal_habitat))
        
        for (i in seq_along(seasonal_habitat)){
          
          habitat_seasons[[i]] <- tmp$season[tmp$block == i]
          
        }
  
      }
      
      
      if (length(seasonal_habitat) != length(habitat_seasons)) {
        stop("length of seasonal_habitat must equal length of habitat_seasons")
      }
      
      
      
      
      time_step = 1 / seasons
      
      burn_steps <-
        burn_years * seasons # convert years to time steps
      
      season_foo <-
        function(x, seasons, time_step) {
          if (any(x > seasons)) {
            stop("habitat seasons must be less than or equal to the number of seasons")
            
          }
          
          y = (x / seasons) - time_step
          
        }
      

      all_habitat_seaons <- purrr::map_dfr(habitat_seasons, ~ data.frame(season = .x))
      
      if (!all(1:seasons %in% all_habitat_seaons$season)){
        stop("all seasons must be represented in habitat_seasons")
      }
      
      habitat_seasons <-
        purrr::map(habitat_seasons,
                   season_foo,
                   time_step = time_step,
                   seasons = seasons)
      
      if (!is.null(dim(seasonal_habitat[[1]]))) {
        resolution <- nrow(seasonal_habitat[[1]])
        
      }
      
      patches <- resolution ^ 2
      
      if (!is.na(scientific_name)) {
        common_name <-
          taxize::sci2comm(scientific_name, db = "ncbi")[[1]][1]
        
      }
      
      if (is.na(scientific_name) &
          !is.na(common_name)) {
        scientific_name <-
          taxize::comm2sci(common_name, db = "ncbi")[[1]][1]
        
      }
      
      # check fishbase -------------
      if (is.na(scientific_name) == F &
          query_fishlife == T) {
        sq <- purrr::safely(purrr::quietly(get_traits))
        
        
        genus_species <-
          stringr::str_split(scientific_name, " ", simplify = T) %>%
          as.data.frame() %>%
          rlang::set_names(c("genus", "species"))
        
        fish_life <- genus_species %>%
          dplyr::mutate(life_traits = purrr::pmap(list(
            Genus = genus, Species = species
          ),
          sq))
        
        
        if (!is.null(fish_life$life_traits[[1]]$error)) {
          stop("No match in FishLife: check spelling or supply your own life history values")
        }
        
        fish_life <- fish_life %>%
          dplyr::mutate(fish_life_worked = purrr::map(life_traits, 'error') %>% purrr::map_lgl(is.null)) %>%
          dplyr::filter(fish_life_worked) %>%
          dplyr::mutate(life_traits = purrr::map(life_traits, c('result', "result"))) %>%
          tidyr::unnest(cols = life_traits) %>%
          dplyr::mutate(taxa = glue::glue('{genus} {species}')) %>%
          rlang::set_names(tolower)
        
        
        if (weight_units == "kg") {
          fish_life$winfinity <- fish_life$winfinity / 1000
        }
        
        if (tune_weight == T) {
          weight_stan <- "
   data {
    real winf;
    real linf;
  }
    parameters {
    real<lower = 0> wa;
    real<lower = 2.7, upper = 3.2> wb;
    real<lower = 0, upper = 1> sigma;
    }
    transformed parameters{
    real w_hat;
    w_hat = wa*linf^wb;

    }
    model {
    winf ~ normal(w_hat, sigma);
    wb ~ normal(3,.1);
    }
    "
          
          weight_fit <-
            rstan::stan(
              model_code = weight_stan,
              data = list(winf = fish_life$winfinity * 2, linf = fish_life$loo),
              verbose = F,
              cores = cores
            )
          
          weight_fit <-
            broom::tidy(weight_fit) %>%
            dplyr::select(term, estimate) %>%
            tidyr::spread(term, estimate)
        } else{
          weight_fit <-
            dplyr::tibble(wa = fish_life$winfinity / (fish_life$loo ^ default_wb),
                          wb = default_wb)
        }
        # process lengths ---------------------------------------------------------
        
        if (is.na(linf)) {
          linf <- fish_life$loo
        }
        
        if (is.na(vbk)) {
          vbk <- fish_life$k
          
        }
        
        if (is.na(weight_a)) {
          weight_a <- weight_fit$wa
          
        }
        if (is.na(weight_b)) {
          weight_b <- weight_fit$wb
          
        }
        
        if (is.na(max_age)) {
          max_age <- ceiling(fish_life$tmax)
          
        }
        
        if (is.na(age_mature)) {
          age_mature <- fish_life$tm
          
        }
        
        if (is.na(length_mature)) {
          length_mature <- fish_life$lm
          
        }
        
        if (is.na(m)) {
          m <- fish_life$m
          
        }
        
      } #close fishlife query
      
      length_at_age <-
        linf * (1 - exp(-vbk * (seq(
          min_age, max_age, by = time_step
        ) - t0)))
      
      # process weight
      
      weight_at_age <-
        weight_a * length_at_age ^ weight_b
      
      lmat_to_linf_ratio <- length_mature / linf
      
      m_at_age <-
        rep(m * time_step, length(weight_at_age)) # place holder to allow for different m at age
      
      #
      
      length_at_age_key <-
        generate_length_at_age_key(
          min_age = min_age,
          max_age = max_age,
          cv = cv_len,
          linf = linf,
          k = vbk,
          t0 = t0,
          time_step = time_step,
          linf_buffer = linf_buffer
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(age, length_bin, p_bin) %>%
        tidyr::spread(length_bin, p_bin) %>%
        dplyr::select(-age)
      
      # process maturity
      if ((is.na(age_50_mature) |
           is.na(age_95_mature)) &
          is.na(age_mature) == F) {
        age_50_mature <- age_mature
        
        age_95_mature <-
          age_50_mature + delta_mature
        
        maturity_at_age <-
          ((1 / (1 + exp(-log(
            19
          ) * ((seq(min_age, max_age, by = time_step) - age_50_mature) / (age_95_mature - age_50_mature)
          )))))
        
      } else if (is.na(age_mature) |
                 mat_mode == "length") {
        if (is.na(length_mature)) {
          length_mature <-  linf * lmat_to_linf_ratio
        }
        
        length_bins <-
          as.numeric(colnames(length_at_age_key))
        
        mat_at_bin <- ((1 / (1 + exp(-log(
          19
        ) * ((length_bins - length_mature) / (delta_mature)
        )))))
        
        p_mat_at_age <-
          (as.matrix(length_at_age_key) %*% mat_at_bin)
        
        mat_at_age <-
          dplyr::tibble(age = seq(min_age, max_age, by = time_step),
                        mean_mat_at_age = p_mat_at_age)
        
        age_mature <-
          mat_at_age$age[mat_at_age$mean_mat_at_age >= 0.5][1]
        
        age_50_mature <- age_mature
        
        age_95_mature <-
          mat_at_age$age[mat_at_age$mean_mat_at_age >= 0.95][1]
        
        maturity_at_age <-
          mat_at_age$mean_mat_at_age
      }
      
      
      if (is.na(length_50_mature)) {
        length_50_mature <- length_mature
        
        length_95_mature <-
          length_50_mature + delta_mature
        
      }
      
      self$maturity_at_age <- maturity_at_age
      
      self$weight_at_age <- weight_at_age
      
      self$ssb_at_age <-
        maturity_at_age * weight_at_age
      
      # create habitat and movement matrices
      
      calc_move_mat <-
        function(habitat,
                 movement,
                 movement_sigma,
                 time_step,
                 resolution) {
          distance <-
            tidyr::expand_grid(x = 1:resolution, y = 1:resolution) %>%
            dist() %>%
            as.matrix()
          
          p_move <-
            dnorm(distance, movement * time_step, movement_sigma * time_step) # movement represents amount moved in a year, so scale down by time step
          
          p_move <- p_move / rowSums(p_move)
          
          if ((!is.null(dim(habitat)))) {
            habitat <- habitat %>%
              as.data.frame() %>%
              dplyr::mutate(x = 1:nrow(.)) %>%
              tidyr::pivot_longer(
                -x,
                names_to = "y",
                values_to = "habitat",
                names_prefix = "V",
                names_ptypes = list(habitat = integer())
              )
            
            move_habitat <- matrix(
              rep(habitat$habitat, resolution),
              nrow = resolution ^ 2,
              ncol = resolution ^ 2,
              byrow = TRUE
            )
            
            move_habitat <-
              move_habitat / rowSums(move_habitat)
            
            move_mat <-  p_move * move_habitat
            
            move_mat <-
              move_mat / rowSums(move_mat)
            
          } else {
            move_mat <- p_move
          }
          
          move_mat <-
            t(move_mat) # needs to be transposed for use in population function
          
        } # close calc_move_mat
      
      self$seasonal_movement <-
        purrr::map(
          seasonal_habitat,
          calc_move_mat,
          movement = adult_movement,
          movement_sigma = adult_movement_sigma,
          resolution = resolution,
          time_step = time_step
        )
      
      self$movement_seasons <- habitat_seasons
      
      self$dd_movement <- dd_movement
      
      # set up unfished recruitment by patch
      if (is.null(dim(rec_habitat))) {
        rec_habitat <- matrix(1, nrow = resolution, ncol = resolution)
        
      }
      
      r0s <- rec_habitat %>%
        as.data.frame() %>%
        dplyr::mutate(x = 1:nrow(.)) %>%
        tidyr::pivot_longer(
          -x,
          names_to = "y",
          values_to = "rec_habitat",
          names_prefix = "V",
          names_ptypes = list(rec_habitat = integer())
        ) %>%
        dplyr::mutate(rec_habitat = rec_habitat / sum(rec_habitat))
      
      self$r0s <- r0 * r0s$rec_habitat
      
      # tune SSB0 and unfished equilibrium
      init_pop <-
        self$r0s * matrix(
          rep(exp(-m * seq(0, max_age, by = time_step)), patches),
          nrow = patches,
          ncol = length(length_at_age),
          byrow = TRUE
        )
      f_p_a <-
        matrix(0, nrow = patches, ncol = length(length_at_age))
      
      self$max_age <- max_age # in years
      
      self$min_age <- min_age # starting age for the model
      
      self$linf <- linf
      
      self$m <- m
      
      self$density_dependence_form <-
        density_dependence_form
      
      self$time_step <- time_step
      
      self$seasons <- seasons
      
      self$rec_form <- rec_form
      
      self$common_name <- common_name
      
      self$scientific_name <- scientific_name
      
      self$length_at_age <- length_at_age
      
      self$length_at_age_key <- length_at_age_key
      
      self$length_50_mature <- length_50_mature
      
      self$patches <- patches
      
      self$steepness <- steepness
      
      self$m_at_age <- m_at_age
      
      self$fished_depletion <- fished_depletion
      
      self$explt_type <- explt_type
      
      self$init_explt <- init_explt
      
      unfished <- marlin::sim_fish(
        length_at_age = length_at_age,
        weight_at_age = weight_at_age,
        maturity_at_age = maturity_at_age,
        steepness = steepness,
        m_at_age = m_at_age,
        patches = resolution ^ 2,
        burn_steps = burn_steps,
        time_step = time_step,
        season = 0,
        r0s = self$r0s,
        ssb0 = NA,
        ssb0_p = rep(-999, patches),
        f_p_a = f_p_a,
        seasonal_movement = self$seasonal_movement,
        movement_seasons = self$movement_seasons,
        dd_movement = self$dd_movement,
        last_n_p_a = init_pop,
        tune_unfished = 1,
        rec_form = rec_form
      )
      self$ssb0 <- unfished$ssb0
      
      self$ssb0_p <- unfished$ssb0_p
      
      self$n_p_a_0 <- unfished$tmppop$n_p_a
      
      self$unfished <- unfished$tmppop
      
     
      
    }, # close initialize
    plot = function(type = 2) {
     
      tmp <- as.list(self)
      
      ogives <- tmp[which(str_detect(names(tmp), "at_age$"))]
      
      tidy_ogives <-
        purrr::map_df(ogives, ~ data.frame(
          val_at_age = .x,
          age = seq(self$min_age, self$max_age, by = self$time_step)
        ), .id = "trait")
      
      tidy_ogives %>%
        ggplot2::ggplot(aes(age, val_at_age, color = trait)) +
        ggplot2::geom_line(show.legend = FALSE) +
        ggplot2::facet_wrap( ~ trait, scales = "free_y") + 
        marlin::theme_marlin() + 
        ggplot2::scale_fill_manual(values = marlin::marlin_pal("diverging_fish")(length(unique(tidy_ogives$trait))))    
      
    },
    #' Swim
    #'
    #' Swim advances the population one time step
    #'
    #' @param burn_steps number of steps for burn in period if applicable
    #' @param season the current season
    #' @param f_p_a matrix of fishing mortality by patch and age
    #' @param last_n_p_a matrix of initial numbers by patch and age
    #' @param tune_unfished boolean indicating whether to tune unfished
    #'
    #' @return the population in the next time step
    swim = function(burn_steps = 0,
               season = 1,
               f_p_a = NULL,
               last_n_p_a = NULL,
               tune_unfished = 0) {
      
      season <- (season / self$seasons) - self$time_step
        
        if (is.null(f_p_a)) {
          f_p_a <-
            matrix(0,
                   nrow = self$patches,
                   ncol = length(self$length_at_age))
          
        }
        
        if (is.null(last_n_p_a)) {
          last_n_p_a <- self$n_p_a_0
        }
        pop <- marlin::sim_fish(
          length_at_age = self$length_at_age,
          weight_at_age = self$weight_at_age,
          maturity_at_age = self$maturity_at_age,
          steepness = self$steepness,
          m_at_age = self$m_at_age,
          patches = self$patches,
          burn_steps = burn_steps,
          time_step = self$time_step,
          season = season,
          r0s = self$r0s,
          ssb0 = self$ssb0,
          ssb0_p = self$ssb0_p,
          seasonal_movement = self$seasonal_movement,
          movement_seasons = self$movement_seasons,
          dd_movement = self$dd_movement,
          f_p_a = f_p_a,
          last_n_p_a = last_n_p_a,
          tune_unfished = tune_unfished,
          rec_form = self$rec_form
        )
        
      } # close swim
    
  ) # close public
) # close object