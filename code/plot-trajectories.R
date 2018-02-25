plot_sims <- function(final_sizes,
                      initial_sizes,
                      output = c("circumference", "carbon_simple", "carbon_complex"),
                      nyear_plot = NULL,
                      nsim_plot = 25,
                      file = NULL,
                      table_file = NULL,
                      correction_cm = NULL,
                      sp_list,
                      climate_scenarios,
                      sp_rich_set = 1,
                      plot_dim) {
  
  # subset final and initial sizes based on species richness settings
  final_sizes <- final_sizes[, , sp_rich_set, , ]
  initial_sizes <- initial_sizes[, , sp_rich_set, ]
  
  # open pdf file if requested
  if (!is.null(file)) 
    pdf(file = file, height = 7, width = 7)
  
  # plot settings
  old_mfrow <- par()$mfrow
  if (is.null(plot_dim))
    plot_dim <- c(2, 3)
  par(mfrow = plot_dim)
  
  # colour settings
  col_set <- viridis(256, alpha = 1)[floor(seq(250, 50, length = length(climate_scenarios)))]
  col_set_fine <- viridis(256, alpha = 0.15)[floor(seq(250, 50, length = length(climate_scenarios)))]
  
  # setup x axis
  if (is.null(nyear_plot))
    nyear_plot <- dim(final_sizes)[3] / 12
  time_seg <- c(1:(12 * nyear_plot))
  x_ids <- 2017 + seq(1, nyear_plot, length = length(time_seg))
  
  # calculate carbon if needed
  if (output != "circumference") {
    ylab_set <- "Carbon (kg)"
    
    # set correction to zero if not supplied
    if (is.null(correction_cm))
      correction_cm <- rep(0, length(sp_list))
    
    # calculate carbon values with correction factors for bark width
    if (output == "carbon_simple") {
      final_vals <- convert_to_carbon(final_sizes, correction_cm = correction_cm)
      initial_vals <- convert_to_carbon(initial_sizes, correction_cm = correction_cm)
    } else {
      final_vals <- array(NA, dim = dim(final_sizes))
      initial_vals <- array(NA, dim = dim(initial_sizes))
      for (i in seq_len(dim(final_sizes)[1])) {
        final_vals[i, , , ] <- convert_to_carbon_spp(final_sizes[i, , , ],
                                                     correction_cm = correction_cm[i],
                                                     sp_id = sp_list[i])
        initial_vals[i, , ] <- convert_to_carbon_spp(initial_sizes[i, , ],
                                                     correction_cm = correction_cm[i],
                                                     sp_id = sp_list[i])
      }
    }
    
  } else {
    ylab_set <- "Circumference (mm)"
    final_vals <- final_sizes
    initial_vals <- initial_sizes
  }
  
  # plot trajectories
  for (i in seq_along(sp_list)) {
    plot_sub <- sample(seq_len(dim(final_vals)[4]), size = nsim_plot, replace = FALSE)
    ylim_set <- range(final_vals[i, , time_seg, plot_sub], na.rm = TRUE)
    plot(final_vals[i, 1, time_seg, 1] ~ x_ids,
         bty = "l'", type = "n", las = 1,
         xlab = "Year", ylab = ylab_set,
         xlim = c(x_ids[1] - diff(x_ids)[1], max(x_ids)),
         ylim = ylim_set)
    for (j in seq_along(climate_scenarios)) {
      for (k in plot_sub) {
        vals_to_plot <- c(initial_vals[i, j, k],
                           final_vals[i, j, time_seg, k])
        lines(vals_to_plot ~ c(x_ids[i] - diff(x_ids)[1], x_ids),
              col = col_set_fine[j], lwd = 0.7)
      }
    }
    for (j in seq_along(climate_scenarios)) {
      mean_vals <- apply(final_vals[i, j, time_seg, ], 1, mean, na.rm = TRUE)
      mean_vals <- c(mean(initial_vals[i, j, ]), mean_vals)
      lines(mean_vals ~ c(x_ids[1] - diff(x_ids)[1], x_ids),
            col = col_set[j], lwd = 4)
    }
    mtext(sp_list[i], side = 3, line = 0.5, adj = 0, cex = 1.25)
  }
  
  # reset plot settings
  par(mfrow = old_mfrow)
  
  if (!is.null(file))
    dev.off()
  
  # print results to a CSV
  if (!is.null(table_file)) {
    store_vals <- apply(final_vals[, , dim(final_vals)[3], ], c(1, 2), mean)
    colnames(store_vals) <- c("dry", "ave", "wet")
    rownames(store_vals) <- sp_list
    write.csv(store_vals,
              file = table_file)
  }
  
}