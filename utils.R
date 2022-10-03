gplot_data <- function(x, maxpixels = 10e4)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as_tibble(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], by = c('value' = 'ID'))
  }
  dat
}

SIR_solver <- function(beta, gamma, C, SIR0, times) {
  
  SIR_function <- function(time, variables, parameters) {
    n_compart <- 3
    n_classes <- length(variables) / n_compart
    S <- as.matrix(variables[1:n_classes])
    I <- as.matrix(variables[(n_classes + 1):(2 * n_classes)])
    R <- as.matrix(variables[(2 * n_classes + 1):(3 * n_classes)])
    
    with(as.list(parameters), {
      N_population <- S + I + R
      dS <- -1 * as.matrix(beta * S) * (as.matrix(C) %*% as.matrix(I / N_population))
      dI <- -1 * dS - gamma * as.matrix(I)
      dR <- gamma * as.matrix(I)
      return(list(c(dS, dI, dR)))
    })
  }
  
  parameters_values <- c(
    beta  = beta, 
    gamma = gamma, 
    C = C
  )
  
  initial_values <- c(
    S = SIR0$S,  
    I = SIR0$I,  
    R = SIR0$R   
  )
  
  out <- ode(
    y = initial_values,
    times = times,
    func = SIR_function,
    parms = parameters_values
  )
  
  as.data.frame(out)
}

generate_SIR <- function (N_population, C, SIR0, beta, gamma, N = 100) {
  times <- 1:N
  n_classes <- nrow(C)
  
  SIR <- SIR_solver(beta = beta, gamma = gamma, C = C, SIR0 = SIR0, times = times)
  SIR[, 2:(n_classes + 1)] <- floor(SIR[, 2:(n_classes + 1)])
  SIR[, (n_classes + 2):(n_classes * 2 + 1)] <- floor(SIR[, (n_classes + 2):(n_classes * 2 + 1)])
  SIR[, (n_classes * 2 + 2):(n_classes * 3 + 1)] <- as.data.frame(matrix(data = N_population, nrow = length(times), ncol = n_classes, byrow = TRUE)) - SIR[, 2:(n_classes + 1)] - SIR[, (n_classes + 2):(n_classes * 2 + 1)]
  SIR
}

plot_SIR <- function (SIR, s, save = TRUE) {
  if (save) { png(filename = paste('output/', sprintf('%02d', s), '/plots/SIR.png', sep = ''), width = 800, height = 600) }
  par(family = 'LM Roman 10', mfrow = c(1, 1))
  plot(NA, xlim = c(1, tail(SIR$time, 1)), ylim = c(0, SIR$S[1]), xlab = 'Time', ylab = 'Number of individuals')
  lines(SIR$S, col = 'blue',  lwd = 2)
  lines(SIR$I, col = 'red',   lwd = 2)
  lines(SIR$R, col = 'green', lwd = 2)
  if (save) { dev.off() }
}

plot_SIR_sep <- function (SIR_sep, prop_class, s, save = TRUE) {
  n_classes <- length(prop_class)
  if (save) { png(filename = paste('output/', sprintf('%02d', s), '/plots/SIR_sep.png', sep = ''), width = 800, height = 600) }
  par(family = 'LM Roman 10', mfrow = c(1, 1))
  plot(NA, xlim = c(1, tail(SIR_sep$time, 1)), ylim = c(0, max(SIR_sep)), xlab = 'Time', ylab = 'Number of individuals')
  for (i in 1:n_classes) {
    lines(SIR_sep[, (i + 1)], col = 'blue',  lwd = 2, lty = (i + 1))
    lines(SIR_sep[, (i + n_classes + 1)], col = 'red',   lwd = 2, lty = (i + 1))
    lines(SIR_sep[, (i + (2 * n_classes) + 1)], col = 'green', lwd = 2, lty = (i + 1))
  }
  if (save) { dev.off() }
}


plot_infect_locations <- function (area_pop, Terminal, map, infect_locations, s, save = TRUE) {
  mx_legend <- ceiling(max(gplot_data(area_pop)$value) / 100) * 100
  n_classes <- length(infect_locations)
  
  select_pt <- function (selected_window, n_classes) {
    for (i in 1:n_classes) {
      if (i == 1) {
        pts <- as.data.frame(infect_locations[[i]][[selected_window]])
        pts$class <- rep(x = i, times = nrow(pts))
      } else {
        partial <- as.data.frame(infect_locations[[i]][[selected_window]])
        partial$class <- rep(x = i, times = nrow(partial))
        pts <- rbind(pts, partial)
      }
    }
    pts
  }
  
  for (selected_window in 1:Terminal) {
    print(selected_window)
    m <- ggmap(map) +
      geom_tile(data = gplot_data(area_pop), aes(x = x, y = y, fill = value), alpha = 0.5) +
      geom_point(data = select_pt(selected_window = selected_window, n_classes = n_classes), aes(x = x, y = y, color = as.factor(class)), alpha = 0.5) +
      guides(color = 'none') +
      scale_fill_gradientn(colors = rainbow(n = mx_legend / 100, start = 0.1, end = 0.9), name = 'Population', limits = c(0, mx_legend), breaks = seq(0, mx_legend, 100)) +
      labs(x = 'Longitude', y = 'Latitude') +
      theme(text = element_text(family = 'LM Roman 10'), legend.key.width = unit(0.5, 'cm'), legend.key.height = unit(1.29, 'cm'),
            panel.background = element_rect(fill = 'transparent', color = NA), plot.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
    if (save) { ggsave(filename = paste('output/', sprintf('%02d', s), '/plots/maps/plot_', selected_window, '.png', sep = ''), plot = m, width = 3000, height = 1000, units = 'px', dpi = 300, bg = 'transparent') } else { print(m) }
  }
}

plot_estimated_infectious <- function (Y_hat, SIR, N_restricted, n_classes, s, save = TRUE) {
  
  Y_hat_cp <- Y_hat
  
  for (k in 1:n_classes) {
    Y_hat <- Y_hat_cp[, , c(k, (n_classes + k), (n_classes * 2 + k))]
    
    Y_hat_mean <- apply(Y_hat, c(2, 3), mean) 
    Y_hat_.025 <- apply(Y_hat, c(2, 3), quantile, prob = c(0.025))
    Y_hat_.975 <- apply(Y_hat, c(2, 3), quantile, prob = c(0.975))
    M <- Y_hat_mean[, 2]; L <- Y_hat_.025[, 2]; U <- Y_hat_.975[, 2];
    max_y <- max(Y_hat_.975[, 2]) + 100
    
    dfI <- data.frame(t = SIR$SIR$time, M = M, L = L, U = U)
    
    if (save) { png(filename = paste('output/', sprintf('%02d', s), '/plots/Infect_group_', k, '.png', sep = ''), width = 800, height = 600) }
    par(family = 'LM Roman 10', mfrow = c(1, 1))
    plot(NA, xlim = c(0, nrow(SIR$SIR)), ylim = c(0, max_y), main = paste('Group ', k, sep = ''), xlab = 'Time', ylab = 'Number of individuals', xaxs = 'i', yaxs = 'i')
    polygon(c(dfI$t, rev(dfI$t)), c(dfI$L, rev(dfI$U)), col = rgb(1, 0, 0, alpha = 0.1), border = FALSE)
    lines(x = SIR$SIR$time, y = SIR$SIR_sep[, (n_classes + k + 1)], col = 2, lwd = 3)
    lines(x = SIR$SIR$time, y = dfI$M, col = 2, lty = 2, lwd = 3)
    abline(v = N_restricted, lty = 2)
    legend(x = 'topleft', legend = c('Infected', 'Est. Infected'), col = rep('red', 2), lty = c(1, 2), lwd = rep(2, 2), cex = 0.75, bg = 'white')
    if (save) { dev.off() }
  }
}

convert_result <- function (area_pop, val, n_row_count, n_col_count) {
  r <- raster(nrow = n_row_count, ncol = n_col_count)
  extent(r) <- extent(area_pop)
  mult_factor <- (n_row_count / nrow(area_pop)) * (n_col_count / ncol(area_pop))
  values(r) <- raster::extract(x = area_pop, y = rasterToPoints(r)) /  mult_factor
  m <- matrix(data = values(r), nrow = n_row_count, ncol = n_col_count, byrow = TRUE)
  v <- as.vector(m)
  v[which(!is.na(v))] <- val
  m <- matrix(data = v, nrow = n_row_count, ncol = n_col_count, byrow = )
  values(r) <- m
  r
}

n_inf_func <- function (n_classes, N, processed_result, area_pop) {
  
  cellarea <- prod(res(area_pop))
  
  inf_fit <- list()
  for (k in 1:n_classes) {
    part <- c()
    for (t in 1:N) {
      part <- c(part, sum(values(processed_result[[k]]$result_r_mean[[t]]) * cellarea)
      )
    }
    inf_fit[[k]] <- part
  }
  
  as.data.frame(do.call(cbind, inf_fit))
}

SIR_obs_gen <- function (SIR, intensities, area_pop) {
  SIR_obs <- SIR
  
  for (k in 1:3) {
    for (t in 1:Terminal) {
      SIR_obs$SIR_sep[t, ((3 + 1) + k)] <- sum(values(intensities[[k]][[t]]) * prod(res(area_pop)))
    }
  }
  
  SIR_obs$SIR$I <- rowSums(SIR$SIR_sep[, 5:7])
  
  SIR_obs
}
