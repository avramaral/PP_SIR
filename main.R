source('libraries.R')
source('utils.R')
source('data_generation.R')
source('intensities.R')
source('SIR_estimation.R')
source('spatioTemporal_estimation.R')
source('error_analysis.R')

source('scenarios.R')

start    <- 1
Terminal <- 100
delta    <- 1
I0       <- c(1, 1)

mult_factor <- 1

nu    <- 1
scale <- 0.05
sig_2 <- 0.2
mu    <- -1 * sig_2 / 2 
sd    <- 0.1
a     <- 0.5


prop_class <- c(0.25, 0.75) # 25% of kids and 75% of adults
C <- matrix(data = c(18, 9, 3, 12), nrow = length(prop_class), ncol = length(prop_class), byrow = TRUE)

S <- length(scenarios)
for (s in 1:S) {
  
  print(paste(sprintf('%02d', s), ' out of ', S, sep = ''))
  
  inv_phi    <- scenarios[[s]]$inv_phi
  beta       <- scenarios[[s]]$beta
  gamma      <- scenarios[[s]]$gamma
  model      <- scenarios[[s]]$model_generation
  null_model <- scenarios[[s]]$null_model
  AR_include <- scenarios[[s]]$AR_include
  
  set.seed(0)
  
  ###################
  # DATA GENERATION #
  ###################
  
  area_pop <- generate_area() # You can safely ignore error and warning messages from this function.
  map <- live_map(area_pop = area_pop, zoom = 16, maptype = 'roadmap')
  
  N_population <- sum(values(area_pop), na.rm = TRUE)
  
  SIR <- simulate_SIR(start = start, Terminal = Terminal, delta = delta, N_population = N_population, prop_class = prop_class, C = C, I0 = I0, beta = beta, gamma = gamma, phi = (1 / inv_phi))
  plot_SIR(SIR = SIR$SIR, s = s) # PLOT
  plot_SIR_sep(SIR_sep = SIR$SIR_sep, prop_class = prop_class, s = s) # PLOT
  
  int_process <- generate_intensity(area_pop = area_pop, SIR_sep = SIR$SIR_sep, prop_class = prop_class, nu = nu, scale = scale, sig_2 = sig_2, mu = mu, sd = sd, a = a, model = model)
  intensities <- list()
  for (i in 1:length(prop_class)) { intensities[[i]] <- int_process[[i]]$intensities }
  
  infect_locations <- simulate_locations(SIR_sep = SIR$SIR_sep, intensities = intensities) 
  # plot_infect_locations(area_pop = area_pop, Terminal = Terminal, map = map, infect_locations = infect_locations, s = s) # PLOT
  
  saveRDS(object = area_pop, file = paste('output/', sprintf('%02d', s), '/rds/area_pop.rds', sep = ''))
  saveRDS(object = SIR, file = paste('output/', sprintf('%02d', s), '/rds/SIR.rds', sep = ''))
  saveRDS(object = intensities, file = paste('output/', sprintf('%02d', s), '/rds/intensities.rds', sep = ''))
  saveRDS(object = infect_locations, file = paste('output/', sprintf('%02d', s), '/rds/infect_locations.rds', sep = ''))
  
  #####################
  # TEMPORAL MODELING #
  #####################
  
  N <- nrow(SIR$SIR)
  N_restricted <- 50
  
  fit <- SIR_NB(start = start, Terminal = Terminal, N = N, N_restricted = N_restricted, N_population = N_population, SIR = SIR, prop_class = prop_class, C = C)
  parameters <- mixing_assessment(fit = fit)
  parameters_sum <- parameters_summary(parameters = parameters)
  Y_hat <- generate_Y_hat(parameters = parameters, SIR_sep = SIR$SIR_sep, prop_class = prop_class, C = C, N_population = N_population)
  
  est_pars_file <- file(paste('output/', sprintf('%02d', s), '/txt/parameters_summary.txt', sep = ''))
  writeLines(as.character(parameters_sum), est_pars_file)
  close(est_pars_file)
  
  plot_estimated_infectious(Y_hat = Y_hat, SIR = SIR, N_restricted = N_restricted, n_classes = length(prop_class), s = s) # PLOT
  
  saveRDS(object = fit, file = paste('output/', sprintf('%02d', s), '/rds/fit_temporal.rds', sep = ''))

  ############################
  # SPATIO-TEMPORAL MODELING #
  ############################

  n_row_count <- nrow(area_pop) * mult_factor
  n_col_count <- ncol(area_pop) * mult_factor
  
  count_cells <- counting_events(area_pop = area_pop, infect_locations = infect_locations, start = start, Terminal = Terminal, delta = delta, n_row = n_row_count, n_col = n_col_count)
    
  n_classes <- length(prop_class)
  
  result <- list()
  processed_result <- list()
  for (k in 1:n_classes) {
    print(paste('Classes: ', sprintf('%02d', k), ' out of ', sprintf('%02d', n_classes), sep = ''))
    area_pop_tmp <- area_pop
    values(area_pop_tmp) <- values(area_pop) * prop_class[k]
    Y_hat_tmp <- Y_hat[, , c(k, (n_classes + k), (n_classes * 2 + k))]
    result[[k]] <- fit_spatioTemporal(area_pop = area_pop_tmp, count_cells = count_cells[[k]], Y_hat = Y_hat_tmp, N_restricted = N_restricted, n_row_count = n_row_count, n_col_count = n_col_count, null_model = null_model, AR_include = AR_include)
    processed_result[[k]] <- process_result(result = result[[k]], area_pop = area_pop_tmp, n_row_count = n_row_count, n_col_count = n_col_count) # Unchanged
  }
  
  saveRDS(object = count_cells, file = paste('output/', sprintf('%02d', s), '/rds/count_cells.rds', sep = ''))
  saveRDS(object = result, file = paste('output/', sprintf('%02d', s), '/rds/result_spatioTemporal.rds', sep = ''))
  saveRDS(object = processed_result, file = paste('output/', sprintf('%02d', s), '/rds/processed_result.rds', sep = '')) 
  
  ##################
  # ERROR ANALYSIS #
  ##################
  
  computed_error_MAPE  <- compute_error(intensities = intensities, processed_result = processed_result, error_type = 'MAPE' )
  computed_error_MAAPE <- compute_error(intensities = intensities, processed_result = processed_result, error_type = 'MAAPE')
  
  for (k in 1:n_classes) {
    # plot_error(computed_error = computed_error_MAPE[[k]],  selected_class = k, s = s, error_type = 'MAPE' ) # PLOT
      plot_error(computed_error = computed_error_MAAPE[[k]], selected_class = k, s = s, error_type = 'MAAPE') # PLOT
  }

  saveRDS(object = computed_error_MAPE, file = paste('output/', sprintf('%02d', s), '/rds/computed_error_MAPE.rds',  sep = ''))
  saveRDS(object = computed_error_MAAPE, file = paste('output/', sprintf('%02d', s), '/rds/computed_error_MAAPE.rds', sep = ''))
}
