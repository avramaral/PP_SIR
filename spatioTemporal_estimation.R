counting_events <- function (area_pop, infect_locations, start, Terminal, delta, n_row = NULL, n_col = NULL) {
  
  if (is.null(n_row) | is.null(n_col)) {
    n_row <- nrow(area_pop)
    n_col <- ncol(area_pop)
  }
  
  area <- raster(nrow = n_row, ncol = n_col)
  extent(area) <- extent(area_pop)
  values(area) <- 0
  
  b <- st_as_sf(st_union(st_as_sf(rasterToPolygons(area_pop), merge = T)))
  
  area <- raster::rasterize(x = b, y = area, mask = TRUE)
  
  ts <- seq(from = start, to = Terminal, by = delta)
  progressbar <- txtProgressBar(min = 1, max = length(ts), initial = 1) 
  
  r <- list()
  grid_r <- list()
  area_p <- rasterToPolygons(area)
  for (t in ts) {
    r[[t]] <- area
    r[[t]][] <- 0
    if (infect_locations[[t]]$n != 0) {
      pp <- as.SpatialPoints.ppp(infect_locations[[t]])
      tab <- table(cellFromXY(object = area, xy = pp))
      r[[t]][as.numeric(names(tab))] <- tab
    }
    grid_r[[t]] <- rasterToPolygons(r[[t]])
    grid_r[[t]]$id <- 1:nrow(grid_r[[t]])
    
    grid_r_temp <- raster::intersect(x = grid_r[[t]], y = area_p)
    grid_r[[t]] <- grid_r[[t]][grid_r[[t]]$id %in% grid_r_temp$id, ]
    
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  
  grid_r
}

fit_spatioTemporal <- function (area_pop, count_cells, Y_hat, N_restricted, n_row_count = 10, n_col_count = 30, null_model = TRUE, AR_include = FALSE, verbose = FALSE) {
  
  ts <- 1:length(count_cells)
  
  convert_counts <- function (ts, area_pop, count_cells, N_restricted, n_row_count, n_col_count) {
    N <- length(count_cells)
    ts <- ts[1:N_restricted]
    
    progressbar <- txtProgressBar(min = 1, max = length(ts), initial = 1) 
    
    counts <- c() 
    for (t in ts) {
      r <- raster(nrow = n_row_count, ncol = n_col_count)
      extent(r) <- extent(area_pop)
      r <- rasterize(x = count_cells[[t]], y = r, field = 'layer')
      m <- matrix(data = values(r), nrow = n_row_count, ncol = n_col_count, byrow = TRUE)
      counts <- c(counts, na.omit(as.vector(m)))
      setTxtProgressBar(progressbar, t)
    }
    counts <- c(counts, rep(NA, nrow(count_cells[[1]]) * (N - length(ts))))
    
    close(progressbar)
    counts
  }
  
  convert_pop <- function (area_pop, N_tm, n_row_count, n_col_count) {
    r <- raster(nrow = n_row_count, ncol = n_col_count)
    extent(r) <- extent(area_pop)
    mult_factor <- ((n_row_count / nrow(area_pop)) * (n_col_count / ncol(area_pop)))
    m <- matrix(data = raster::extract(x = area_pop, y = rasterToPoints(r)), nrow = n_row_count, ncol = n_col_count, byrow = TRUE) / mult_factor
    v <- na.omit(as.vector(m))
    v <- rep(v, N_tm)
    v / sum(v, na.rm = TRUE)
  }
  
  id_sp <- 1:nrow(count_cells[[1]]); N_sp <- length(id_sp)
  id_tm <- ts; N_tm <- length(id_tm)
  cellarea <- prod(res(area_pop)) / ((n_row_count / nrow(area_pop)) * (n_col_count / ncol(area_pop)))
  
  print('Creating data object...')
  data_INLA <- data.frame(id = 1:(N_sp * N_tm),
                          id_sp = rep(id_sp, N_tm),
                          id_tm = rep(id_tm, each = N_sp),
                          counts = convert_counts(ts = ts, area_pop = area_pop, count_cells = count_cells, N_restricted = N_restricted, n_row_count = n_row_count, n_col_count = n_col_count),
                          infect = rep(apply(Y_hat, c(2, 3), mean)[, 2], each = N_sp),
                          lambda_0 = convert_pop(area_pop = area_pop, N_tm = N_tm, n_row_count = n_row_count, n_col_count = n_col_count) / cellarea)

  if (null_model) {
    if (!AR_include) {
      formula <- counts ~ 1 + f(id,  model = 'iid')  + f(id_sp, model = 'matern2d', nrow = n_row_count, ncol = n_col_count, nu = 1)
    } else {
      formula <- counts ~ 1 + f(id,  model = 'iid')  + f(id_sp, model = 'matern2d', nrow = n_row_count, ncol = n_col_count, nu = 1) + f(id_tm,  model = 'ar1')
    }
  } else {
    if (!AR_include) {
      formula <- counts ~ 0 + offset(log(infect * lambda_0)) + f(id,  model = 'iid')  + f(id_sp, model = 'matern2d', nrow = n_row_count, ncol = n_col_count, nu = 1)
    } else {
      formula <- counts ~ 0 + offset(log(infect * lambda_0)) + f(id,  model = 'iid')  + f(id_sp, model = 'matern2d', nrow = n_row_count, ncol = n_col_count, nu = 1) + f(id_tm,  model = 'ar1') 
    }
  }
  
  result <- inla(formula = formula,
                 family = 'poisson',
                 data = data_INLA,
                 # E = rep(cellarea, (N_sp * N_tm)),
                 control.predictor = list(compute = TRUE, link = 1),
                 control.compute = list(config = TRUE),
                 verbose = verbose, safe = TRUE)

  result
}

process_result <- function (result, area_pop, n_row_count, n_col_count) {
  
  N_sp <- n_row_count * n_col_count
  ts <- 1:(nrow(result$summary.fitted.values) / N_sp)
  cellarea <- prod(res(area_pop)) / ((n_row_count / nrow(area_pop)) * (n_col_count / ncol(area_pop)))
  
    convert_result <- function (area_pop, values, n_row_count, n_col_count) {
      r <- raster(nrow = n_row_count, ncol = n_col_count)
      extent(r) <- extent(area_pop)
      mult_factor <- (n_row_count / nrow(area_pop)) * (n_col_count / ncol(area_pop))
      values(r) <- raster::extract(x = area_pop, y = rasterToPoints(r)) /  mult_factor
      m <- matrix(data = values(r), nrow = n_row_count, ncol = n_col_count, byrow = TRUE)
      v <- as.vector(m)
      v[which(!is.na(values(r)))] <- values
      m <- matrix(data = v, nrow = n_row_count, ncol = n_col_count)
      values(r) <- m
      r
    }
    
  quant <- c('0.025quant', '0.975quant', 'mean')
  fitted_values <- list()
  result_r_0250 <- list()
  result_r_0975 <- list()
  result_r_mean <- list()
  for (t in ts) {
    fitted_values[[t]] <- result$summary.fitted.values[((t - 1) * N_sp + 1):((t) * N_sp), quant] / cellarea
    result_r_0250[[t]] <- convert_result(area_pop = area_pop, values = fitted_values[[t]][, 1], n_row_count = n_row_count, n_col_count = n_col_count)
    result_r_0975[[t]] <- convert_result(area_pop = area_pop, values = fitted_values[[t]][, 2], n_row_count = n_row_count, n_col_count = n_col_count)
    result_r_mean[[t]] <- convert_result(area_pop = area_pop, values = fitted_values[[t]][, 3], n_row_count = n_row_count, n_col_count = n_col_count)
  }
  
  list(result_r_0250 = result_r_0250, result_r_0975 = result_r_0975, result_r_mean = result_r_mean)
}