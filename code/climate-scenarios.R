# helper functions to define climate scenarios

# Scenario options:
#    - wet
#    - dry
#    - ave
#    - dry punctuated by wet
#    - alternating wet and dry
#
# Based on data for wet, average, and dry years
# transition <- matrix(c(dry_to_dry, dry_to_ave, dry_to_wet,
#                        ave_to_dry, ave_to_ave, ave_to_wet,
#                        wet_to_dry, wet_to_ave, wet_to_wet),
#                      ncol = 3,
#                      byrow = FALSE)
#

## CALCULATE LONG-TERM STATE PROBABILITIES FOR EACH

# wet transition matrix based on RCP4.5 CESMI-BGC
wet <- matrix(c(0.4, 0.4, 0.2,
                0.19, 0.37, 0.44,
                0.1, 0.37, 0.53),
              ncol = 3,
              byrow = FALSE)

# dry transition matrix based on RCP4.5 CSIRO-MK3-6-0
dry <- matrix(c(0.84, 0.14, 0.02,
                0.32, 0.5, 0.18,
                0.17, 0.33, 0.5),
              ncol = 3,
              byrow = FALSE)

# ave transition matrix
ave <- matrix(c(0.1, 0.8, 0.1,
                0.2, 0.6, 0.2,
                0.1, 0.8, 0.1),
              ncol = 3,
              byrow = FALSE)

# dry punctuated by wet transition matrix based on Bendigo 1992 to 2016
dry2wet <- matrix(c(0.6, 0.27, 0.13,
                    0.6, 0.2, 0.2,
                    0.79, 0.01, 0.1),
                  ncol = 3,
                  byrow = FALSE)

# alternating wet and dry transition matrix based on RCP4.5 ACCESS1.3
wet2dry <- matrix(c(0.58, 0.23, 0.19,
                    0.31, 0.23, 0.46,
                    0.5, 0.33, 0.17),
                  ncol = 3,
                  byrow = FALSE)

# function to update climate state
update_climate <- function(state, transmat, nsim) {
  
  for (st in seq_along(state)) {
    state[st] <- sample(c(1:3),
                        size = 1,
                        prob = transmat[, state[st]])
  }
  
  state
  
}

prepare_climate_data <- function(climate_in) {
  
  climate_data <- vector("list", length = 3)
  
  climate_data[[1]] <- list(rain_summer = climate_in$rainfall[which((climate_in$level == "Below") &
                                                                      (climate_in$season == "Summer"))],
                            rain_autumn = climate_in$rainfall[which((climate_in$level == "Below") &
                                                                      (climate_in$season == "Autumn"))],
                            rain_winter = climate_in$rainfall[which((climate_in$level == "Below") &
                                                                      (climate_in$season == "Winter"))],
                            rain_spring = climate_in$rainfall[which((climate_in$level == "Below") &
                                                                      (climate_in$season == "Spring"))],
                            maxt_summer = climate_in$maxt[which((climate_in$level == "Above") &
                                                                  (climate_in$season == "Summer"))],
                            maxt_autumn = climate_in$maxt[which((climate_in$level == "Above") &
                                                                  (climate_in$season == "Autumn"))],
                            maxt_winter = climate_in$maxt[which((climate_in$level == "Above") &
                                                                  (climate_in$season == "Winter"))],
                            maxt_spring = climate_in$maxt[which((climate_in$level == "Above") &
                                                                  (climate_in$season == "Spring"))])
  climate_data[[2]] <- list(rain_summer = climate_in$rainfall[which((climate_in$level == "Average") &
                                                                      (climate_in$season == "Summer"))],
                            rain_autumn = climate_in$rainfall[which((climate_in$level == "Average") &
                                                                      (climate_in$season == "Autumn"))],
                            rain_winter = climate_in$rainfall[which((climate_in$level == "Average") &
                                                                      (climate_in$season == "Winter"))],
                            rain_spring = climate_in$rainfall[which((climate_in$level == "Average") &
                                                                      (climate_in$season == "Spring"))],
                            maxt_summer = climate_in$maxt[which((climate_in$level == "Average") &
                                                                  (climate_in$season == "Summer"))],
                            maxt_autumn = climate_in$maxt[which((climate_in$level == "Average") &
                                                                  (climate_in$season == "Autumn"))],
                            maxt_winter = climate_in$maxt[which((climate_in$level == "Average") &
                                                                  (climate_in$season == "Winter"))],
                            maxt_spring = climate_in$maxt[which((climate_in$level == "Average") &
                                                                  (climate_in$season == "Spring"))])
  climate_data[[3]] <- list(rain_summer = climate_in$rainfall[which((climate_in$level == "Above") &
                                                                      (climate_in$season == "Summer"))],
                            rain_autumn = climate_in$rainfall[which((climate_in$level == "Above") &
                                                                      (climate_in$season == "Autumn"))],
                            rain_winter = climate_in$rainfall[which((climate_in$level == "Above") &
                                                                      (climate_in$season == "Winter"))],
                            rain_spring = climate_in$rainfall[which((climate_in$level == "Above") &
                                                                      (climate_in$season == "Spring"))],
                            maxt_summer = climate_in$maxt[which((climate_in$level == "Below") &
                                                                  (climate_in$season == "Summer"))],
                            maxt_autumn = climate_in$maxt[which((climate_in$level == "Below") &
                                                                  (climate_in$season == "Autumn"))],
                            maxt_winter = climate_in$maxt[which((climate_in$level == "Below") &
                                                                  (climate_in$season == "Winter"))],
                            maxt_spring = climate_in$maxt[which((climate_in$level == "Below") &
                                                                  (climate_in$season == "Spring"))])
  
  # return outputs
  climate_data
  
}

sample_seasons <- function(x, climate_data, standards) {
  
  out <- list()
  out$rain <- c(sample(climate_data[[x]]$rain_summer, size = 3, replace = TRUE),
                sample(climate_data[[x]]$rain_autumn, size = 3, replace = TRUE),
                sample(climate_data[[x]]$rain_winter, size = 3, replace = TRUE),
                sample(climate_data[[x]]$rain_spring, size = 3, replace = TRUE))
  out$maxt <- c(sample(climate_data[[x]]$maxt_summer, size = 3, replace = TRUE),
                sample(climate_data[[x]]$maxt_autumn, size = 3, replace = TRUE),
                sample(climate_data[[x]]$maxt_winter, size = 3, replace = TRUE),
                sample(climate_data[[x]]$maxt_spring, size = 3, replace = TRUE))
  out$rain <- (out$rain - standards$rain_mean) / standards$rain_sd
  out$maxt <- (out$maxt - standards$maxt_mean) / standards$maxt_sd
  out
  
}

sample_climates <- function(climate_set,
                            climate_data,
                            standards) {
  
  # initialise climate
  climate <- list()
  
  # sample rainfall and maxt data
  climate_tmp <- lapply(climate_set, sample_seasons,
                        climate_data,
                        standards)
  climate$rain <- sapply(climate_tmp, function(x) x$rain)
  climate$maxt <- sapply(climate_tmp, function(x) x$maxt)
  
  # return outputs
  climate
  
}
