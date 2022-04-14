generate_intensity <- function (area_pop, SIR, nu, scale, sig_2, mu, sd, a, model) {
  
  ts <- SIR$time
  n_row_count <- nrow(area_pop) * mult_factor
  n_col_count <- ncol(area_pop) * mult_factor
  
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
    
    cvModel <- RMmatern(nu = nu, var = sig_2, scale = scale) + RMtrend(mean = mu)
    
    process <- list()
    intensities <- list()
    for (t in ts) {
      p <- raster(RFsimulate(model = cvModel, x = x_seq, y = y_seq))
      r <- ref
      values(r)[pos] <- values(p)[pos]
      values(r) <- exp(values(r) + rnorm(n = 1, mean = 0, sd = sd))
      process[[t]] <- r
      
      tmp <- (process[[t]] * ref)
      intensities[[t]] <- ((tmp / sum(values(tmp), na.rm = TRUE)) * SIR$I[t]) / prod(res(ref))
    }
  } else if (model == 'AR1') {
    mu <- -1 * ((sig_2 / (1 - a ** 2)) + (sd ** 2)) / 2 
    
    cvModel <- RMmatern(nu = nu, var = sig_2, scale = scale) + RMtrend(mean = 0)
    
    p <- raster(RFsimulate(model = cvModel, x = x_seq, y = y_seq))
    r <- ref
    values(r)[pos] <- values(p)[pos]
    
    process <- list()
    process[[1]] <- r
    intensities <- list()
    for (t in ts[-1]) {
      p <- raster(RFsimulate(model = cvModel, x = x_seq, y = y_seq))
      values(r) <- a * values(process[[(t - 1)]]) + values(p) # Preserve NA's, if needed
      process[[t]] <- r
    }
    
    process <- lapply(X = process, FUN = function (gp) { values(gp) <- exp(values(gp) + mu + rnorm(n = 1, mean = 0, sd = sd)); gp })
    
    for (t in ts) { 
      tmp <- (process[[t]] * ref)
      intensities[[t]] <- ((tmp / sum(values(tmp), na.rm = TRUE)) * SIR$I[t]) / prod(res(ref))
    }
  } else { stop('Choose a valid model.') }
  
  list(process = process, intensities = intensities)
}
