compute_error <- function (intensities, results, error_type = 'MAPE') {
  
  N <- length(intensities)
  errors <- c()
  pts <- rasterToPoints(x = results[[1]])[, 1:2]
  
  progressbar <- txtProgressBar(min = 1, max = N, initial = 1) 
  for (i in 1:N) {
    int_val <- raster::extract(x = intensities[[i]], y = pts)
    if (error_type == 'MAPE') {
      errors <- c(errors,  compute_MAPE(int_val = int_val, res_val = values(results[[i]])))
    } else if (error_type == 'MAAPE') {
      errors <- c(errors, compute_MAAPE(int_val = int_val, res_val = values(results[[i]])))
    }
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)
  
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

plot_error <- function (computed_error, s, error_type = 'MAPE') {
  png(filename = paste('output/', sprintf('%02d', s), '/plots/error_', error_type, '.png', sep = ''), width = 800, height = 600)
  par(family = 'LM Roman 10', mfrow = c(1, 1))
  if (error_type == 'MAAPE') { y_max <- 1.58 } else { y_max <- max(computed_error) }
  plot(computed_error, type = 'l', ylim = c(0, y_max), xlab = 'Time', ylab = error_type, main = paste(error_type, ' (', sprintf('%02d', s), ')', sep = ''))
  # lines(computed_error, col = 'red')
  abline(v = N_restricted, lty = 2)
  # legend(x = 'topleft', legend = c('Alternative model', 'Null model'), col = c('black', 'red'), lty = 1)
  dev.off()
}
