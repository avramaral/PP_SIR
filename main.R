source('header.R')

for (s in 1:S) {
  
  print(paste(sprintf('%02d', s), ' out of ', S, sep = ''))
  
  inv_phi    <- scenarios[[s]]$inv_phi
  beta       <- scenarios[[s]]$beta
  gamma      <- scenarios[[s]]$gamma
  model      <- scenarios[[s]]$model_generation
  null_model <- scenarios[[s]]$null_model
  AR_include <- scenarios[[s]]$AR_include
  
  set.seed(159)
  
  ###################
  # DATA GENERATION #
  ###################
  
  area_pop <- generate_area() 
  map <- live_map(area_pop = area_pop, zoom = 16, maptype = 'roadmap')
  
  N_population <- sum(values(area_pop), na.rm = TRUE)
  
  if (s %in% seq( 1, 16)) { SIR <- simulate_SIR(start = start, Terminal = Terminal, delta = delta, N_population = N_population, prop_class = prop_class, C = C, I0 = I0, beta = beta, gamma = gamma, phi = (1 / inv_phi)) }
  if (s %in% seq(17, 24)) { SIR <- acnb(gamma) }
  plot_SIR(SIR = SIR$SIR, s = s, save = F) 
  plot_SIR_sep(SIR_sep = SIR$SIR_sep, prop_class = prop_class, s = s, save = F) 
  
  int_process <- generate_intensity(area_pop = area_pop, SIR_sep = SIR$SIR_sep, prop_class = prop_class, nu = nu, scale = scale, sig_2 = sig_2, mu = mu, sd = sd, sd_AR1 = sd_AR1, a = a, model = model)
  intensities <- list()
  for (i in 1:length(prop_class)) { intensities[[i]] <- int_process[[i]]$intensities }
  
  infect_locations <- simulate_locations(SIR_sep = SIR$SIR_sep, intensities = intensities) 
  
  saveRDS(object = area_pop, file = paste('output/', sprintf('%02d', s), '/rds/area_pop.rds', sep = ''))
  saveRDS(object = SIR, file = paste('output/', sprintf('%02d', s), '/rds/SIR.rds', sep = ''))
  saveRDS(object = intensities, file = paste('output/', sprintf('%02d', s), '/rds/intensities.rds', sep = ''))
  saveRDS(object = infect_locations, file = paste('output/', sprintf('%02d', s), '/rds/infect_locations.rds', sep = ''))
  
  #####################
  # TEMPORAL MODELING #
  #####################
  
  SIR <- SIR_obs_gen(SIR = SIR, intensities = intensities, area_pop = area_pop)
  SIR$SIR <- round(SIR$SIR)
  SIR$SIR_sep <- round(SIR$SIR_sep)
  
  N <- nrow(SIR$SIR)
  N_restricted <- 50
  
  fit <- SIR_NB(start = start, Terminal = Terminal, N = N, N_restricted = N_restricted, N_population = N_population, SIR = SIR, prop_class = prop_class, C = C)
  parameters <- mixing_assessment(fit = fit)
  parameters_sum <- parameters_summary(parameters = parameters)
  Y_hat <- generate_Y_hat(parameters = parameters, SIR_sep = SIR$SIR_sep, prop_class = prop_class, C = C, N_population = N_population)
  
  est_pars_file <- file(paste('output/', sprintf('%02d', s), '/txt/parameters_summary.txt', sep = ''))
  writeLines(as.character(parameters_sum), est_pars_file)
  close(est_pars_file)
  
  plot_estimated_infectious(Y_hat = Y_hat, SIR = SIR, N_restricted = N_restricted, n_classes = length(prop_class), s = s, save = F) # PLOT
  
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
    processed_result[[k]] <- process_result(result = result[[k]], area_pop = area_pop_tmp, n_row_count = n_row_count, n_col_count = n_col_count) 
  }
  
  saveRDS(object = count_cells, file = paste('output/', sprintf('%02d', s), '/rds/count_cells.rds', sep = ''))
  saveRDS(object = result, file = paste('output/', sprintf('%02d', s), '/rds/result_spatioTemporal.rds', sep = ''))
  saveRDS(object = processed_result, file = paste('output/', sprintf('%02d', s), '/rds/processed_result.rds', sep = '')) 
  
}
