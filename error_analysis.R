compute_error <- function (intensities, processed_result, error_type = 'MAPE') {
  
  n_classes <- length(intensities)
  
  errors <- list() 
  for (k in 1:n_classes) {
    results <- processed_result[[k]]$result_r_mean
    
    N <- length(intensities[[k]])
    error <- c()
    pts <- rasterToPoints(x = results[[1]])[, 1:2]
    
    progressbar <- txtProgressBar(min = 1, max = N, initial = 1) 
    for (i in 1:N) {
      int_val <- raster::extract(x = intensities[[k]][[i]], y = pts)
      if (error_type == 'MAPE') {
        error <- c(error,  compute_MAPE(int_val = int_val, res_val = values(results[[i]])))
      } else if (error_type == 'MAAPE') {
        error <- c(error, compute_MAAPE(int_val = int_val, res_val = values(results[[i]])))
      }
      setTxtProgressBar(progressbar, i)
    }
    close(progressbar)
    
    errors[[k]] <- error
  }
  
  errors
}

compute_MAAPE <- function (int_val, res_val) {
  N <- length(int_val)
  r <- atan(abs((int_val - res_val) / (int_val)))
  sum(r) / N
}


compute_MAPE <- function (int_val, res_val) {
  N <- length(int_val)
  r <- (abs((int_val - res_val) / (int_val)))
  sum(r) / N
}

plot_error <- function (computed_error, selected_class, s, error_type = 'MAPE', save = TRUE) {
  if (save) { png(filename = paste('output/', sprintf('%02d', s), '/plots/error_', error_type, '_class_', selected_class, '.png', sep = ''), width = 800, height = 600) }
  par(family = 'LM Roman 10', mfrow = c(1, 1))
  if (error_type == 'MAAPE') { y_max <- 1.58 } else { y_max <- max(computed_error) }
  plot(computed_error, type = 'l', ylim = c(0, y_max), xlab = 'Time', ylab = error_type, main = paste(error_type, ' (', sprintf('%02d', s), ')', sep = ''))
  # lines(computed_error, col = 'red')
  abline(v = N_restricted, lty = 2)
  # legend(x = 'topleft', legend = c('Alternative model', 'Null model'), col = c('black', 'red'), lty = 1)
  if (save) { dev.off() }
}
