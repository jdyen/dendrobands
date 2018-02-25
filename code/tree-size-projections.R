# project dendrobands based on fitted model and climate scenarios

# optional: clear workspace
# rm(list = ls())

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
if (!require(viridis)) {
  install.packages("viridis")
  library(viridis)
}
  
# load required functions
source("./code/projection-helpers.R")
source("./code/climate-scenarios.R")
source("./code/plot-trajectories.R")

# load data
tree_data <- get(load("./data/tree_data.R"))
data_set <- get(load("./data/dendroband_data.R"))
climate_scenarios <- read.csv("./data/climate_scenarios.csv",
                              stringsAsFactors = FALSE)
climate_data <- prepare_climate_data(climate_scenarios)
standards_list <- get(load("./data/data_standardizations.R"))

# load fitted models
fitted_model <- get(load("./outputs/fitted_model_summary.R"))

# setup scenario options
nsim <- 10              # number of replicate runs
nyear <- 10             # number of years to project over
nmonth <- nyear * 12

# set input parameters
scenario <- list(nsim = nsim,
                 nyear = nyear,
                 ndays = 30)
climate_options <- c("dry", "wet2dry", "wet", "ave")
sp_rich_options <- c(1, 5, 10)
nfixers_present <- TRUE
sp_list <- levels(factor(tree_data$species))
params <- extract_fitted(fitted_model,
                         sp_names = sp_list)

## parallelise species/climates
final_sizes <- array(NA, dim = c(length(sp_list),
                                 length(climate_options),
                                 length(sp_rich_options),
                                 (12 * nyear),
                                 scenario$nsim))
initial_sizes <- array(NA, dim = c(length(sp_list),
                                   length(climate_options),
                                   length(sp_rich_options),
                                   scenario$nsim))
tree_size_range <- tapply(tree_data$initial_circum[which(tree_data$site_age < 10)],
                          tree_data$species[which(tree_data$site_age < 10)],
                          range, na.rm = TRUE)
tree_size_range$BM <- tree_size_range$RG
tree_size_range$LW <- tree_size_range$SW
tree_size_range$RSB <- tree_size_range$GB
tree_size_range$SG <- tree_size_range$RG
tree_size_range$YG <- tree_size_range$RG
for (i in seq_along(sp_list)) {
  sp_set <- sp_list[i]
  
  # set nitrogen fixers
  if (nfixers_present) {
    scenario$nfix <- 1
  } else {
    scenario$nfix <- 0
  }
  # set nitrogen fixers to present if LW, SW or ALLO are being modelled
  if (!is.na(match(sp_set, c("LW", "SW", "ALLO"))))
    scenario$nfix <- 1
  
  # set initial tree sizes for each species (from sites <10 years old)
  scenario$tree_size <- runif(scenario$nsim,
                              min = tree_size_range[[sp_set]][1],
                              max = tree_size_range[[sp_set]][2])
  
  # add scenario details
  scenario$sp_id <- i
  tree_params <- list()
  
  for (rich in seq_along(sp_rich_options)) {
    
    # set species richness scenario
    scenario$sp_rich <- sp_rich_options[rich]
    
    for (j in seq_along(climate_options)) {
      
      # reset tree size after each scenario
      tree_params$size <- scenario$tree_size
      
      # reset site age after each scenario
      scenario$site_age <- 5
      
      # randomised burn in to generate initial climate state for each scenario
      nburn <- 50
      clim_scen <- climate_options[sample(seq_along(climate_options), size = 1)]
      climate_set <- sample(1:3, size = scenario$nsim, replace = TRUE)
      for (ii in seq_len(nburn)) {
        climate_set <- update_climate(state = climate_set,
                                      transmat = get(clim_scen),
                                      nsim = scenario$nsim)
      }
      
      initial_sizes[scenario$sp_id, j, rich, ] <- tree_params$size
      clim_scen <- climate_options[j]
      
      for (k in 1:nyear) {
        
        climate_set <- update_climate(state = climate_set,
                                      transmat = get(clim_scen))
        climate <- sample_climates(climate_set,
                                   climate_data,
                                   standards_list)
        
        for (m in 1:12) {
          
          # predict one time step update from fitted params
          tree_params <- predict_growth(params = params,
                                        scenario = scenario,
                                        climate = lapply(climate, function(x) x[m, ]),
                                        tree_params = tree_params,
                                        standards = standards_list)
          
          final_sizes[scenario$sp_id, j, rich, (k - 1) * 12 + m, ] <- tree_params$size
          
        }
        
        # update site age after each year
        scenario$site_age <- scenario$site_age + 1
        
      }
      
    }
    
  }
  
}

# plot projected tree sizes over full time frame
plot_sims(final_sizes, 
          initial_sizes,
          output = "circumference",  # option: "circumference" for circumference,
                                     #         "carbon_simple" for basic carbon calc,
                                     #         "carbon_complex" for species-specific carbon calculation
          nyear_plot = nyear,
          nsim_plot = 5,             # number of simulated trajectories to plot (more = slower)
          file = NULL,               # give a filename for plots, e.g. "circumference_trajectories.pdf"
          table_file = NULL,         # give a filename for output values, e.g. "circumference_trajectories.csv"
          correction_cm = NULL,      # set the bark correction (in cm) for all species
          sp_list = sp_list,
          climate_scenarios = climate_options,
          sp_rich_set = 3,           # set the index of sp_rich_options that you want to plot
          plot_dim = c(1, 1))        # can change to have more than one plot per page: plot_dim = c(no_rows, no_columns)
