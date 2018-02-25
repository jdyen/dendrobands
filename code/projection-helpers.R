# extract outputs from fitted stan model
extract_fitted <- function(x, sp_names = NULL) {
  
  # extract means and SDs from fitted model
  alpha <- x[grep("alpha$", rownames(x)), c("mean", "sd")]
  alpha_sp <- x[grep("alpha_sp\\[", rownames(x)), c("mean", "sd")]
  beta_rain <- x[grep("beta_rain\\[", rownames(x)), c("mean", "sd")]
  beta_maxt <- x[grep("beta_maxt\\[", rownames(x)), c("mean", "sd")]
  beta_uthor <- x[grep("beta_uthor\\[", rownames(x)), c("mean", "sd")]
  beta_soilk <- x[grep("beta_soilk\\[", rownames(x)), c("mean", "sd")]
  beta_rich <- x[grep("beta_rich\\[", rownames(x)), c("mean", "sd")]
  beta_age <- x[grep("beta_age*", rownames(x)), c("mean", "sd")]
  beta_size <- x[grep("beta_size", rownames(x)), c("mean", "sd")]
  beta_nfix <- x[grep("beta_nfix", rownames(x)), c("mean", "sd")]
  beta_dens <- x[grep("beta_dens", rownames(x)), c("mean", "sd")]
  beta_basal <- x[grep("beta_basal", rownames(x)), c("mean", "sd")]
  beta_soil <- x[grep("beta_soil\\[", rownames(x)), c("mean", "sd")]
  
  # add species names if provided
  if (!is.null(sp_names)) {
    rownames(alpha_sp) <- sp_names
    rownames(beta_rain) <- rownames(beta_maxt) <- rownames(beta_uthor) <- sp_names
    rownames(beta_soilk) <- rownames(beta_rich) <- rownames(beta_age) <- sp_names
    rownames(beta_size) <- rownames(beta_nfix) <- rownames(beta_dens) <- sp_names
    rownames(beta_basal) <- sp_names
  }
  
  # collate outputs
  params <- list(alpha_mean = alpha[1], alpha_sd = alpha[2],
                 alpha_sp_mean = alpha_sp[, 1], alpha_sp_sd = alpha_sp[, 2],
                 beta_rain_mean = beta_rain[, 1], beta_rain_sd = beta_rain[, 2],
                 beta_maxt_mean = beta_maxt[, 1], beta_maxt_sd = beta_maxt[, 2],
                 beta_uthor_mean = beta_uthor[, 1], beta_uthor_sd = beta_uthor[, 2],
                 beta_soilk_mean = beta_soilk[, 1], beta_soilk_sd = beta_soilk[, 2],
                 beta_rich_mean = beta_rich[, 1], beta_rich_sd = beta_rich[, 2],
                 beta_age_mean = beta_age[, 1], beta_age_sd = beta_age[, 2],
                 beta_size_mean = beta_size[, 1], beta_size_sd = beta_size[, 2],
                 beta_nfix_mean = beta_nfix[, 1], beta_nfix_sd = beta_nfix[, 2],
                 beta_dens_mean = beta_dens[, 1], beta_dens_sd = beta_dens[, 2],
                 beta_basal_mean = beta_basal[, 1], beta_basal_sd = beta_basal[, 2],
                 beta_soil_mean = matrix(beta_soil[, 1],
                                         nrow = 3, byrow = TRUE)[1, ],
                 beta_soil_sd = matrix(beta_soil[, 2],
                                       nrow = 3, byrow = TRUE)[1, ])
  
  # return outputs
  params
  
}

# predict one time step update from fitted params
predict_growth <- function(params,
                           scenario,
                           climate,
                           tree_params,
                           standards) {
  
  # unpack parameters
  alpha <- rnorm(scenario$nsim,
                 mean = params$alpha_mean,
                 sd = params$alpha_sd)
  alpha_sp <- rnorm(scenario$nsim,
                    mean = params$alpha_sp_mean[scenario$sp_id],
                    sd = params$alpha_sp_sd[scenario$sp_id])
  beta_rain <- rnorm(scenario$nsim,
                     mean = params$beta_rain_mean[scenario$sp_id],
                     sd = params$beta_rain_sd[scenario$sp_id])
  beta_maxt <- rnorm(scenario$nsim,
                     mean = params$beta_maxt_mean[scenario$sp_id],
                     sd = params$beta_maxt_sd[scenario$sp_id])
  beta_rich <- rnorm(scenario$nsim,
                     mean = params$beta_rich_mean[scenario$sp_id],
                     sd = params$beta_rich_sd[scenario$sp_id])
  beta_size <- rnorm(scenario$nsim,
                     mean = params$beta_size_mean[scenario$sp_id],
                     sd = params$beta_size_sd[scenario$sp_id])
  beta_age <- rnorm(scenario$nsim,
                     mean = params$beta_age_mean[scenario$sp_id],
                     sd = params$beta_age_sd[scenario$sp_id])
  beta_nfix <- rnorm(scenario$nsim,
                     mean = params$beta_nfix_mean[scenario$sp_id],
                     sd = params$beta_nfix_sd[scenario$sp_id])
  
  # calculate growth rates (per day)
  growth_mean <- alpha + alpha_sp +
    beta_rain * climate$rain + # rnorm(scenario$nsim, mean = climate$rain, sd = (diff(range(climate$rain)) / 2)) +
    beta_maxt * climate$maxt + # rnorm(scenario$nsim, mean = climate$maxt, sd = (diff(range(climate$maxt)) / 2)) +
    beta_rich * ((scenario$sp_rich - standards$sp_rich_mean) / standards$sp_rich_sd) +
    beta_age * ((scenario$site_age - standards$site_age_mean) / standards$site_age_sd) +
    beta_size * ((tree_params$size - standards$size_mean[scenario$sp_id]) / standards$size_sd[scenario$sp_id]) +
    beta_nfix * ((scenario$nfix))

  # update size
  tree_params$size <- tree_params$size +
    scenario$ndays * growth_mean # rnorm(scenario$nsim, mean = growth_mean, sd = (diff(range(growth_mean)) / 2))
  
  # set any negative values to zero
  tree_params$size <- ifelse(tree_params$size < 0, 0, tree_params$size)
  
  # return outputs
  tree_params
  
}

# convert circumferences to carbon with species-specific corrections
convert_to_carbon <- function(x, correction_cm) {
  
  diam_cm <- x / (10 * pi)
  diam_cm <- sweep(diam_cm, 1, correction_cm, "-")
  out <- 0.5 * exp(-2.3243) * (diam_cm ** 2.4891)
  
  out
  
}

# convert circumferences to carbon with species-specific corrections and conversions
convert_to_carbon_spp <- function(x, correction_cm, sp_id) {
  
  diam_cm <- x / (10 * pi)
  diam_cm <- diam_cm - correction_cm

  out <- get(paste0(sp_id, "_carbon_conversion"))(diam_cm = diam_cm)

  out
  
}

# general carbon eqn
carbon_convert_internal <- function(diam_cm, a, b, cf) {
  
  0.5 * cf * exp((b * log(diam_cm)) + a)
  
}

# species-specific conversions
ALLO_carbon_conversion <- function(diam_cm) {

  # set parameters
  a <- -3.02
  b <- 2.51
  cf <- 1.00494
  size_max <- Inf
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out

}

# species-specific conversions
BM_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.62
  b <- 2.26
  cf <- 0.97384
  size_max <- 94.25
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
BRG_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.83
  b <- 2.15
  cf <- 1.08072
  size_max <- 62.83
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
GB_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.92
  b <- 2.36
  cf <- 1.16698
  size_max <- 345.58
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
IB_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -2.39
  b <- 2.40
  cf <- 1.09881
  size_max <- 188.5
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
LW_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.53
  b <- 2.13
  cf <- 1.02590
  size_max <- 47.12
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
RB_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.46
  b <- 2.05
  cf <- 1.07403
  size_max <- 78.54
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
RG_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.66
  b <- 2.20
  cf <- 1.42022
  size_max <- 219.91
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
RSB_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.71
  b <- 2.21
  cf <- 1.28937
  size_max <- Inf
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
SG_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.36
  b <- 2.30
  cf <- 1.12948
  size_max <- Inf
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
SW_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.21
  b <- 2.11
  cf <- 1.05513
  size_max <- 97.39
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
WB_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.71
  b <- 2.21
  cf <- 1.28937
  size_max <- Inf
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}

# species-specific conversions
YG_carbon_conversion <- function(diam_cm) {
  
  # set parameters
  a <- -1.37
  b <- 2.07
  cf <- 1.03770
  size_max <- 78.54
  
  # calculate carbon
  out <- ifelse((pi * diam_cm) < size_max,
                carbon_convert_internal(diam_cm, a = a, b = b, cf = cf),
                carbon_convert_internal(diam_cm, a = -1.71, b = 2.21, cf = 1.28937))
  
  # return outputs
  out
  
}
