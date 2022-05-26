generate_intensity <- function (area_pop, SIR_sep, prop_class, nu, scale, sig_2, mu, sd, sd_AR1, a, model) {
  n_classes <- length(prop_class)
  
  ts <- SIR_sep$time
  n_row_count <- nrow(area_pop) * mult_factor
  n_col_count <- ncol(area_pop) * mult_factor
  
  result <- list()
  for (i in 1:n_classes) {
    ref <- area_pop
    ref <- raster(nrow = n_row_count, ncol = n_col_count)
    extent(ref) <- extent(area_pop)
    inner_cells <- ((n_row_count / nrow(area_pop)) * (n_col_count / ncol(area_pop)))
    new_pop_val <- matrix(data = raster::extract(x = area_pop, y = rasterToPoints(ref)), nrow = n_row_count, ncol = n_col_count, byrow = TRUE) / inner_cells
    values(ref) <- as.vector(t(new_pop_val))
    
    pos <- which(!is.na(values(ref)))
    
    ex <- extent(ref)
    rs <- res(ref)
    x_seq <- seq(from = ex[1] + (rs[1] / 2), to = ex[2] - (rs[1] / 2), by = rs[1])
    y_seq <- seq(from = ex[3] + (rs[2] / 2), to = ex[4] - (rs[2] / 2), by = rs[2])
    
    if (model == 'IID') {
      mu <- -1 * (sig_2 + (sd ** 2)) / 2 
      
      cvModel <- RMmatern(nu = nu, var = sig_2, scale = scale) + RMtrend(mean = 0)
      
      process <- list()
      intensities <- list()
      progressbar <- txtProgressBar(min = 1, max = length(ts), initial = 1) 
      for (t in ts) {
        p <- raster(RFsimulate(model = cvModel, x = x_seq, y = y_seq))
        r <- ref
        values(r)[pos] <- values(p)[pos]
        values(r) <- exp(values(r) + mu + rnorm(n = 1, mean = 0, sd = sd))
        process[[t]] <- r
        
        tmp <- (process[[t]])
        partial <- (ref  * tmp) 
        intensities[[t]] <- ((partial / sum(values(partial), na.rm = TRUE)) / prod(res(area_pop))) * SIR_sep[t, (i + n_classes + 1)]
        setTxtProgressBar(progressbar, t)
      }
      close(progressbar)
    } else if (model == 'AR1') {
      mu <- -1 * (sig_2 + ((sd_AR1 ** 2) / (1 - a ** 2)) + (sd ** 2)) / 2 
      
      cvModel <- RMmatern(nu = nu, var = sig_2, scale = scale) + RMtrend(mean = 0)
      
      process <- list()
      intensities <- list()
      AR1_process <- arima.sim(list(order = c(1, 0, 0), ar = a), n = length(ts), innov = rnorm(n = length(ts), mean = 0, sd = sd_AR1))
      progressbar <- txtProgressBar(min = 1, max = length(ts), initial = 1) 
      for (t in ts) {
        p <- raster(RFsimulate(model = cvModel, x = x_seq, y = y_seq))
        r <- ref
        values(r)[pos] <- values(p)[pos]
        values(r) <- exp(values(r) + mu + rnorm(n = 1, mean = 0, sd = sd) + AR1_process[t])
        process[[t]] <- r
        
        tmp <- (process[[t]])
        partial <- (ref  * tmp)
        intensities[[t]] <- ((partial / sum(values(partial), na.rm = TRUE)) / prod(res(area_pop))) * SIR_sep[t, (i + n_classes + 1)]
        setTxtProgressBar(progressbar, t)
      }
      close(progressbar)
    } else { stop('Choose a valid model.') }
    
    result[[i]] <- list(process = process, intensities = intensities)
  }
  
  result
}
