source('header.R')
library('scoringRules')

log_0 <- function (x, ...) { ifelse(test = x == 0, yes = 0, no = log(x)) }

n_classes <- length(prop_class)
N <- Terminal
N_restricted <- 50

n_row_count <- 10
n_col_count <- 30
n_cells <- n_row_count * n_col_count

# ERRORS

Q <- 100 # Number of simulated intensities

set.seed(1)

# For computing the scoring rules
simulate_posterior <- function () {

  r_s <- list()
  for (s in 1:S) {
    
    # s <- 1
    print(s)
    
    null_model <- FALSE
    if (s %in% seq(1, S, 2)) { null_model <- TRUE }
    
    # Reading files
    result <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/result_spatioTemporal.rds', sep = ''))
    SIR <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/SIR.rds', sep = ''))
    area_pop <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/area_pop.rds', sep = ''))
    
    # Compute necessary quantities
    SIR_I <- SIR$SIR_sep[, (2 + n_classes):(1 + (2 * n_classes))]
    N_population <- sum(values(area_pop), na.rm = TRUE)
    cellarea <- prod(res(area_pop))
    
    r_k <- list()
    for (k in 1:n_classes) {
      
      # k <- 1
      
      simulated_sample <- inla.posterior.sample(n = 100, result = result[[k]])
      
      r_q <- list()
      for (q in 1:Q) {
        
        # q <- 1
        
        r_t <- list()
        for (t in 1:Terminal) {
          
          # t <- 100
          
          sim <- simulated_sample[[q]]$latent[(((t - 1) * n_cells) + 1):(t * n_cells)]
          sim <- convert_result(area_pop = area_pop, val = sim, n_row_count = n_row_count, n_col_count = n_col_count)
          if (null_model) {
            values(sim) <- exp(values(sim)) / cellarea    
          } else {
            values(sim) <- exp(values(sim)) * log((SIR_I[t, k] + 1e-12) * values(area_pop)) / cellarea  
          }
          
          r_t[[t]] <- sim
        }
        r_q[[q]] <- r_t
      }
      r_k[[k]] <- r_q
    }
    r_s[[s]] <- r_k
  }
  
  r_s
}

# sim <- simulate_posterior()
# saveRDS(object = sim, file = 'output/sim_newC.rds')

# sim <- readRDS(file = 'output/sim.rds')


compute_score_es <- function () {
  
  es_s <- list()
  for (s in 1:S) {
    
    print(s)
    
    intensities <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/intensities.rds', sep = ''))
    
    es_k <- list()
    for (k in 1:n_classes) {
      es_t <- c()
      for (t in 1:Terminal) {
        for (q in 1:Q) {
          if (q == 1) {
            sim_int <- values(sim[[s]][[k]][[q]][[t]])
          } else {
            sim_int <- cbind(sim_int, values(sim[[s]][[k]][[q]][[t]]))
          }
        }
        es_t <- c(es_t, es_sample(y = values(intensities[[k]][[t]]), dat = sim_int))
      }
      es_k[[k]] <- es_t
    }
    es_s[[s]] <- es_k
  }
  
  es_s
}
  
es_all <- compute_score_es()



es3k1 <- es_all[[1]][[1]]
es4k1 <- es_all[[2]][[1]]

plot(NA, xlim = c(1, 100), ylim = c(0, max(c(es3k1, es4k1))), main = 'Scenarios 03 and 04', xlab = 'Time', ylab = 'ES')
lines(1:N, es3k1, col = 2)
lines(1:N, es4k1, col = 3)
abline(v = N_restricted)



####################################################
####################################################


# s <- 1
# k <- 1
# 
# 
# result <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/result_spatioTemporal.rds', sep = ''))
# simulated_sample <- inla.posterior.sample(n = Q, result = result[[k]])
# 
# n_classes <- 3
# 
# SIR_I <- SIR$SIR_sep[, (2 + n_classes):(1 + (2 * n_classes))]
# 
# sim <- simulated_sample[[10]]$latent[1:300]
# sim <- convert_result(area_pop = area_pop, val = sim, n_row_count = 10, n_col_count = 30)
# values(sim) <- exp(values(sim)) * log(SIR_I[1, k] * values(area_pop)) / cellarea
# plot(sim)
# 
# sum(values(sim) * cellarea)
# 
# plot(intensities[[k]][[1]])
# plot(processed_result[[k]]$result_r_mean[[1]])
# 
# cellarea <- prod(res(area_pop))
# 
# 
# 
# 
# 
# int <- list()
# for (s in 1:S) {
# 
#   # s = 1
#   print(s)
#   area_pop <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/area_pop.rds', sep = ''))
#   N_population <- sum(values(area_pop), na.rm = TRUE)
# 
#   int_q <- list()
#   for (q in 1:Q) {
#     print(q)
#     inv_phi    <- scenarios[[s]]$inv_phi
#     beta       <- scenarios[[s]]$beta
#     gamma      <- scenarios[[s]]$gamma
#     model      <- scenarios[[s]]$model_generation
#     null_model <- scenarios[[s]]$null_model
#     AR_include <- scenarios[[s]]$AR_include
# 
# 
#     SIR <- simulate_SIR(start = start, Terminal = Terminal, delta = delta, N_population = N_population, prop_class = prop_class, C = C, I0 = I0, beta = beta, gamma = gamma, phi = (1 / inv_phi))
#     int_process <- generate_intensity(area_pop = area_pop, SIR_sep = SIR$SIR_sep, prop_class = prop_class, nu = nu, scale = scale, sig_2 = sig_2, mu = mu, sd = sd, sd_AR1 = sd_AR1, a = a, model = model)
#     int_partial <- list()
#     for (i in 1:length(prop_class)) { int_partial[[i]] <- int_process[[i]]$intensities }
#     int_q[[q]] <- int_partial
#   }
#   int[[s]] <- int_q
# }
# 
# saveRDS(object = int, file = 'output/simulated_intensities.rds')
# 
# k <- 1
# m <- list()
# for (t in 1:Terminal) {
#   for (q in 1:Q) {
#     if (q == 1) {
#       m[[t]] <- values(int[[3]][[q]][[k]][[t]])
#     } else {
#       m[[t]] <- cbind(m[[t]], values(int[[3]][[q]][[k]][[t]]) )
#     }
#   }
# }
# 
# # processed_results: 3
# 
# es <- c()
# for (t in 1:Terminal) {
#   es <- c(es, es_sample(y = values(processed_result[[1]]$result_r_mean[[t]]), dat = m[[t]]))
# }
# 
# es3 <- es
# 
# #######################
# 
# k <- 1
# m <- list()
# for (t in 1:Terminal) {
#   for (q in 1:Q) {
#     if (q == 1) {
#       m[[t]] <- values(int[[4]][[q]][[k]][[t]])
#     } else {
#       m[[t]] <- cbind(m[[t]], values(int[[4]][[q]][[k]][[t]]) )
#     }
#   }
# }
# 
# # processed_results: 4
# 
# es <- c()
# for (t in 1:Terminal) {
#   es <- c(es, es_sample(y = values(processed_result[[1]]$result_r_mean[[t]]), dat = m[[t]]))
# }
# 
# es4 <- es
# 
# plot(es4)
# 
# 
# plot(NA, xlim = c(1, 100), ylim = c(0, max(c(es3, es4))), main = 'Scenarios 03 and 04', xlab = 'Time', ylab = 'RMSE')
# lines(1:N, es3, col = 2)
# lines(1:N, es4, col = 3)
# abline(v = N_restricted)











area_pop <- readRDS(file = paste('output/01/rds/area_pop.rds', sep = ''))
# Compute errors again
for (s in 1:S) {
  print(s)
  #s <- 4
  # Read values
  SIR <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/SIR.rds', sep = ''))
  intensities <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/intensities.rds', sep = ''))
  processed_result <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/processed_result.rds', sep = ''))

  # Compute necessary quantities
  n_inf <- n_inf_func(n_classes = n_classes, N = N, processed_result = processed_result, area_pop = area_pop)
  SIR_I <- SIR$SIR_sep[, (2 + n_classes):(1 + (2 * n_classes))]

  computed_error_RMSE           <- compute_error(intensities = intensities, processed_result = processed_result, error_type = 'RMSE')
  computed_error_RMSE_norm      <- compute_error(intensities = intensities, processed_result = processed_result, normalize = TRUE, n_inf = n_inf, SIR_I = SIR_I, error_type = 'RMSE')
  computed_error_RMSE_count_inf <- compute_error(intensities = intensities, processed_result = processed_result, count_inf = TRUE, error_type = 'RMSE')

  # computed_error_MAAPE           <- compute_error(intensities = intensities, processed_result = processed_result, error_type = 'MAAPE')
  # computed_error_MAAPE_norm      <- compute_error(intensities = intensities, processed_result = processed_result, normalize = TRUE, n_inf = n_inf, SIR_I = SIR_I, error_type = 'MAAPE')
  # computed_error_MAAPE_count_inf <- compute_error(intensities = intensities, processed_result = processed_result, count_inf = TRUE, error_type = 'MAAPE')


  # plot(NA, xlim = c(1, 100), ylim = c(0, 1.58), main = s)
  # lines(1:N, computed_error_MAAPE[[1]], col = 1)
  # lines(1:N, computed_error_MAAPE_norm[[1]], col = 2)
  # lines(1:N, computed_error_MAAPE_count_inf[[1]], col = 3)
  # abline(v = N_restricted)

  plot(NA, xlim = c(1, 100), ylim = c(0, max(computed_error_RMSE_norm[[1]])), main = s)
  lines(1:N, computed_error_RMSE[[1]], col = 1)
  lines(1:N, computed_error_RMSE_norm[[1]], col = 2)
  lines(1:N, computed_error_RMSE_count_inf[[1]], col = 3)
  abline(v = N_restricted)

}

# Plot for the e-mail
#
# r3 <- computed_error_MAAPE # s = 3
# r4 <- computed_error_MAAPE # s = 4
#
# plot(NA, xlim = c(1, 100), ylim = c(0, 1.58), main = 'Scenarios 03 and 04', xlab = 'Time', ylab = 'MAAPE')
# lines(1:N, r3[[1]], col = 2)
# lines(1:N, r3[[2]], col = 3)
# lines(1:N, r3[[3]], col = 4)
# lines(1:N, r4[[1]], col = 2, lty = 2)
# lines(1:N, r4[[2]], col = 3, lty = 2)
# lines(1:N, r4[[3]], col = 4, lty = 2)
# abline(v = N_restricted)
#
# q3 <- computed_error_RMSE # s = 3
# q4 <- computed_error_RMSE # s = 4
#
# plot(NA, xlim = c(1, 100), ylim = c(0, max(c(unlist(q3), unlist(q4)))), main = 'Scenarios 03 and 04 (Error for Intensity Function)', xlab = 'Time', ylab = 'RMSE')
# lines(1:N, q3[[1]], col = 2)
# lines(1:N, q3[[2]], col = 3)
# lines(1:N, q3[[3]], col = 4)
# lines(1:N, q4[[1]], col = 2, lty = 2)
# lines(1:N, q4[[2]], col = 3, lty = 2)
# lines(1:N, q4[[3]], col = 4, lty = 2)
# abline(v = N_restricted)
#
#
# s3 <- computed_error_RMSE_norm # s = 3
# s4 <- computed_error_RMSE_norm # s = 4
#
# plot(NA, xlim = c(1, 100), ylim = c(0, 7), main = 'Scenarios 03 and 04 (Error for the Counting)', xlab = 'Time', ylab = 'RMSE')
# lines(1:N, s3[[1]], col = 2)
# lines(1:N, s3[[2]], col = 3)
# lines(1:N, s3[[3]], col = 4)
# lines(1:N, s4[[1]], col = 2, lty = 2)
# lines(1:N, s4[[2]], col = 3, lty = 2)
# lines(1:N, s4[[3]], col = 4, lty = 2)
# abline(v = N_restricted)


sum(values(processed_result[[1]]$result_r_mean[[25]]) * cellarea)

