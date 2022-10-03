generate_area <- function (region = 'Brazil', subregion = '3550308', x_center = 4589832.86, y_center = 7318721.75, x_size = 3000, y_size = 1000, nrow = 10, ncol = 30) {

  area <- st_transform(x = read_municipality(code_muni = as.numeric(subregion), year = 2020), crs = CRS('EPSG:4326'))
  pop <- raster(x = paste('data/population/', tolower(region), '.tif', sep = ''))
  # 'brazil.tif': https://drive.google.com/file/d/1nMhIY-DH852qBZDdRnjhhzfT_q3O0fXF/view?usp=sharing
  pop <- crop(x = pop, y = area)
  pop <- rasterize(x = area, y = pop, mask = TRUE)
  
  center_point <- SpatialPoints(cbind(x = x_center, y = y_center), proj4string = CRS('EPSG:5641'))
  
  rect_around_point <- function(x, x_size, y_size){
    bbox <- st_bbox(x)
    bbox <- bbox + c(-x_size / 2, -y_size / 2, x_size / 2, y_size / 2)
    return(st_as_sf(st_as_sfc(bbox)))
  }
  
  created_area <- rect_around_point(x = center_point, x_size = x_size, y_size = y_size) # 3 km2
  created_area <- st_transform(x = created_area, crs = CRS('EPSG:4326'))

  area_pop <- crop(x = pop, y = created_area)
  area_pop <- rasterize(x = created_area, y = area_pop, mask = TRUE)
  area_pop[is.na(area_pop)] <- 0
  N_population <- sum(values(area_pop))
  
  crt_area <- raster(created_area, nrow = nrow, ncol = ncol)
  area_pop <- resample(area_pop, crt_area, method = 'bilinear')
  values(area_pop) <- sapply(X = values(area_pop), FUN = function (x) { max(x, 0) })
  values(area_pop) <- round(values(area_pop) / sum(values(area_pop)) * N_population, 0)
  
  area_pop
}

live_map <- function (area_pop, zoom = 10, maptype = 'terrain') {
  ll <- extent(area_pop)
  get_map(location = c(ll[1], ll[3] , ll[2], ll[4]), zoom = zoom, maptype = maptype)
}

simulate_SIR <- function (start, Terminal, delta, N_population, prop_class, C, I0, beta, gamma, phi) {
  
  validate_contact_matrix <- function (prop_class, C) {
    v <- TRUE
    for (i in 1:length(prop_class)) {
      for (j in 1:length(prop_class)) {
        if ((prop_class[i] * C[i, j]) != (prop_class[j] * C[j, i])) { v <- FALSE }
      } 
    }
    v
  }
  
  times <- seq(from = start, to = Terminal, by = delta)
  N <- length(times)
  S <- matrix(data = rep(x = 0, times = N * length(prop_class)), nrow = N, ncol = length(prop_class), byrow = FALSE)
  I <- matrix(data = rep(x = 0, times = N * length(prop_class)), nrow = N, ncol = length(prop_class), byrow = FALSE)
  R <- matrix(data = rep(x = 0, times = N * length(prop_class)), nrow = N, ncol = length(prop_class), byrow = FALSE)
  
  N_population <- N_population * prop_class
  
  S[1, ] <- N_population - I0
  I[1, ] <- I0
  R[1, ] <- rep(x = 0, times = length(prop_class))
  
  generated_SIR <- generate_SIR(N_population = N_population, C = C, SIR0 = list(S = S[1, ], I = I[1, ], R = R[1, ]), beta = beta, gamma = gamma, N = N)
  I <- rbind(I0, generated_SIR[-1, (length(prop_class) + 2):(length(prop_class) * 2 + 1)])
  for (i in 1:length(prop_class)) {
    I[, i] <- c(I0[i], rnbinom(n = (nrow(I) - 1), size = phi, mu = I[-1, i]))
  }
  R <- floor(generated_SIR[, (length(prop_class) * 2 + 2):(length(prop_class) * 3 + 1)])
  S <- as.data.frame(matrix(data = N_population, nrow = length(times), ncol = length(prop_class), byrow = TRUE)) - I - R
  
  SIR     <- data.frame(time = times, S = rowSums(S), I = rowSums(I), R = rowSums(R))
  SIR_sep <- data.frame(time = times, S = S, I = I, R = R)
  list(SIR = SIR, SIR_sep = SIR_sep)
}

simulate_locations <- function (SIR_sep, intensities) {
  n_classes <- length(intensities)
  times <- SIR_sep$time
  
  result <- list()
  for (i in 1:n_classes) {
    
    progressbar <- txtProgressBar(min = 1, max = length(times), initial = 1) 
    
    infectious_locations <- list()
    for (t in times) {
      if (t == 1) {
        n_points <- rpoint(n = SIR_sep[t, (i + n_classes + 1)], f = as.im(intensities[[i]][[t]] + 1e-6))
        marks(n_points) <- t
        infectious_locations[[t]] <- n_points
      } else {
        rec_diff <- SIR_sep[t, (i + (2 * n_classes) + 1)] - SIR_sep[(t - 1), (i + (2 * n_classes) + 1)]
        if ((SIR_sep[t, (i + n_classes + 1)] != 0) | (SIR_sep[(t - 1), (i + n_classes + 1)] != 0)) { inf_diff <- SIR_sep[t, (i + n_classes + 1)] - SIR_sep[(t - 1), (i + n_classes + 1)] + rec_diff } else { inf_diff <- 0 } # Deal with cases when no infectious occur, but S and R change
        if (inf_diff < 0) {
          rec_diff <- rec_diff + (-1 * inf_diff)
          inf_diff <- 0
        }
        n_points <- rpoint(n = inf_diff, f = as.im(intensities[[i]][[t]] + 1e-6))
        marks(n_points) <- t
        
        if (rec_diff == 0) { 
          old_marks <- NULL
        } else {
          old_marks <- sort(marks(infectious_locations[[(t - 1)]]))[1:rec_diff]
        }
        old_points <- infectious_locations[[(t - 1)]]
        for (o in old_marks) {
          removed_point <- old_points[marks(old_points) == o][1]
          old_points <- old_points[!((old_points$x == removed_point$x) & (old_points$y == removed_point$y))] 
        }
        infectious_locations[[t]] <- superimpose(n_points, old_points)
      }
      setTxtProgressBar(progressbar, t)
    }
    close(progressbar)
    
    result[[i]] <- infectious_locations
  }
  
  result
}
