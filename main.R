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
I0       <- 1

SIR_generation <- 'NB'
SIR_fitting    <- 'NB'

mult_factor <- 1

nu    <- 1
scale <- 0.05
sig_2 <- 0.2
mu    <- -1 * sig_2 / 2 
sd    <- 0.5
a     <- 0.5

S <- length(scenarios)
for (s in 1:S) {
  
  print(paste(sprintf('%02d', s), ' out of ', S, sep = ''))
  
  inv_phi    <- scenarios[[s]]$inv_phi
  alpha      <- scenarios[[s]]$alpha
  beta       <- scenarios[[s]]$beta
  gamma      <- scenarios[[s]]$gamma
  model      <- scenarios[[s]]$model_generation
  null_model <- scenarios[[s]]$null_model
  AR_include <- scenarios[[s]]$AR_include
  
  set.seed(123)
  
  ###################
  # DATA GENERATION #
  ###################
  
  area_pop <- generate_area() # You can safely ignore error and warning messages from this function.
  map <- live_map(area_pop = area_pop, zoom = 16, maptype = 'roadmap')
  
  N_population <- sum(values(area_pop), na.rm = TRUE)
  
  SIR <- simulate_SIR(start = start, Terminal = Terminal, delta = delta, N_population = N_population, I0 = I0, beta = beta, gamma = gamma, method = SIR_generation, phi = (1 / inv_phi))
  plot_SIR(SIR, s = s) # PLOT
  
  int_process <- generate_intensity(area_pop = area_pop, SIR = SIR, nu = nu, scale = scale, sig_2 = sig_2, mu = mu, sd = sd, a = a, model = model)
  intensities <- int_process$intensities
  
  infect_locations <- simulate_locations(SIR = SIR, intensities = intensities) 
  # plot_infect_locations(area_pop = area_pop, Terminal = Terminal, map = map, infect_locations = infect_locations, s = s) # PLOT
  
  saveRDS(object = area_pop, file = paste('output/', sprintf('%02d', s), '/rds/area_pop.rds', sep = ''))
  saveRDS(object = SIR, file = paste('output/', sprintf('%02d', s), '/rds/SIR.rds', sep = ''))
  saveRDS(object = intensities, file = paste('output/', sprintf('%02d', s), '/rds/intensities.rds', sep = ''))
  saveRDS(object = infect_locations, file = paste('output/', sprintf('%02d', s), '/rds/infect_locations.rds', sep = ''))
  
  #####################
  # TEMPORAL MODELING #
  #####################
  
  N <- nrow(SIR)
  N_restricted <- 50
  
  if (SIR_fitting == 'SDE') {
    
    fit <- SIR_SDE(start = start, Terminal = Terminal, delta = delta, N = N, N_restricted = N_restricted, N_population = N_population, SIR = SIR)
    parameters <- mixing_assessment(fit = fit, scaled_beta = TRUE, N_population = N_population)
    parameters_sum <- parameters_summary(parameters = parameters)
    Y_hat <- rstan::extract(fit, pars = c('Y_hat'))$Y_hat
    
  } else if (SIR_fitting == 'NB') {
    
    par_names <- c('phi', 'beta', 'gamma')
    fit <- SIR_NB(start = start, Terminal = Terminal, N = N, N_restricted = N_restricted, N_population = N_population, SIR = SIR)
    parameters <- mixing_assessment(fit = fit, par_names = par_names, scaled_beta = TRUE, N_population = N_population)
    parameters_sum <- parameters_summary(parameters = parameters, par_names = par_names)
    Y_hat <- generate_Y_hat(parameters = parameters, N_population = N_population)
    
  }
  
  est_pars_file <- file(paste('output/', sprintf('%02d', s), '/txt/parameters_summary.txt', sep = ''))
  writeLines(as.character(parameters_sum), est_pars_file)
  close(est_pars_file)
  
  plot_estimated_infectious(Y_hat = Y_hat, SIR = SIR, N_restricted = N_restricted, s = s) # PLOT
  
  saveRDS(object = fit, file = paste('output/', sprintf('%02d', s), '/rds/fit_temporal.rds', sep = ''))

  ############################
  # SPATIO-TEMPORAL MODELING #
  ############################

  n_row_count <- nrow(area_pop) * mult_factor
  n_col_count <- ncol(area_pop) * mult_factor
  
  count_cells <- counting_events(area_pop = area_pop, infect_locations = infect_locations, start = start, Terminal = Terminal, delta = delta, n_row = n_row_count, n_col = n_col_count)
  
  result <- fit_spatioTemporal(area_pop = area_pop, count_cells = count_cells, Y_hat = Y_hat, N_restricted = N_restricted, n_row_count = n_row_count, n_col_count = n_col_count, null_model = null_model, AR_include = AR_include)

  processed_result <- process_result(result = result, area_pop = area_pop, n_row_count = n_row_count, n_col_count = n_col_count)
  
  saveRDS(object = count_cells, file = paste('output/', sprintf('%02d', s), '/rds/count_cells.rds', sep = ''))
  saveRDS(object = result, file = paste('output/', sprintf('%02d', s), '/rds/result_spatioTemporal.rds', sep = ''))
  saveRDS(object = processed_result, file = paste('output/', sprintf('%02d', s), '/rds/processed_result.rds', sep = '')) 
  
  ##################
  # ERROR ANALYSIS #
  ##################
  
  computed_error_MAPE  <- compute_error(intensities = intensities, results = processed_result$result_r_mean, error_type = 'MAPE' )
  computed_error_MAAPE <- compute_error(intensities = intensities, results = processed_result$result_r_mean, error_type = 'MAAPE')
  
  # plot_error(computed_error = computed_error_MAPE , s = s, error_type = 'MAPE' ) # PLOT
  plot_error(computed_error = computed_error_MAAPE, s = s, error_type = 'MAAPE') # PLOT
  
  saveRDS(object = computed_error_MAPE, file = paste('output/', sprintf('%02d', s), '/rds/computed_error_MAPE.rds',  sep = ''))
  saveRDS(object = computed_error_MAAPE, file = paste('output/', sprintf('%02d', s), '/rds/computed_error_MAAPE.rds', sep = ''))
}
