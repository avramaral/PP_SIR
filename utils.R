gplot_data <- function(x, maxpixels = 10e4)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  # Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as_tibble(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], by = c('value' = 'ID'))
  }
  dat
}

SIR_solver <- function(beta, gamma, SIR0, times, N_population) {
  
  SIR_function <- function(time, variables, parameters, N_population) {
    with(as.list(c(variables, parameters)), {
      dS <- -1 * beta * I * S
      dI <- beta * I * S - gamma * I
      dR <- -1 * dS - dI
      return(list(c(dS, dI, dR)))
    })
  }
  
  parameters_values <- c(
    beta  = beta, 
    gamma = gamma 
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
    parms = parameters_values,
    N_population = N_population
  )
  
  as.data.frame(out)
}

generate_SIR <- function (N_population, beta, gamma, N = 100) {
  SIR0  <- list(S = N_population - 1, I = 1, R = 0)
  times <- 1:N
  
  SIR <- SIR_solver(beta = beta, gamma = gamma, SIR0 = SIR0, times = times, N_population = N_population)
  SIR[, 2] <- floor(SIR[, 2])
  SIR[, 3] <- floor(SIR[, 3])
  SIR[, 4] <- N_population - SIR[, 2] - SIR[, 3]
  SIR
}

plot_SIR <- function (SIR, s) {
  png(filename = paste('output/', sprintf('%02d', s), '/plots/SIR.png', sep = ''), width = 800, height = 600)
  par(family = 'LM Roman 10', mfrow = c(1, 1))
  plot(NA, xlim = c(1, tail(SIR$time, 1)), ylim = c(0, SIR$S[1]), xlab = 'Time', ylab = 'Number of individuals')
  lines(SIR$S, col = 'blue',  lwd = 2)
  lines(SIR$I, col = 'red',   lwd = 2)
  lines(SIR$R, col = 'green', lwd = 2)
  dev.off()
}


plot_infect_locations <- function (area_pop, Terminal, map, infect_locations) {
  max_legend <- ceiling(max(gplot_data(area_pop)$value) / 100) * 100
  for (selected_window in 1:Terminal) {
    print(selected_window)
    m <- ggmap(map) +
           geom_tile(data = gplot_data(area_pop), aes(x = x, y = y, fill = value), alpha = 0.5) +
           geom_point(data = as.data.frame(infect_locations[[selected_window]]), aes(x = x, y = y), alpha = 0.5, color = 'red') +
           scale_fill_gradientn(colors = rainbow(n = max_legend / 100, start = 0.1, end = 0.9), name = 'Population', limits = c(0, max_legend), breaks = seq(0, max_legend, 100)) +
           labs(x = 'Longitude', y = 'Latitude') +
           theme(text = element_text(family = 'LM Roman 10'), legend.key.width = unit(0.5, 'cm'), legend.key.height = unit(1.29, 'cm'),
                 panel.background = element_rect(fill = 'transparent', color = NA), plot.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
    ggsave(filename = paste('images/maps/plot_', selected_window, '.png', sep = ''), plot = m, width = 3000, height = 1000, units = 'px', dpi = 300, bg = 'transparent')
  }
}

plot_estimated_infectious <- function (Y_hat, SIR, N_restricted, s) {
  Y_hat_mean <- apply(Y_hat, c(2, 3), mean) 
  Y_hat_.025 <- apply(Y_hat, c(2, 3), quantile, prob = c(0.025))
  Y_hat_.975 <- apply(Y_hat, c(2, 3), quantile, prob = c(0.975))
  M <- Y_hat_mean[, 2]; L <- Y_hat_.025[, 2]; U <- Y_hat_.975[, 2];
  max_y <- max(Y_hat_.975[, 2]) + 100
 
  dfI <- data.frame(t = SIR$time, M = M, L = L, U = U)
  
  png(filename = paste('output/', sprintf('%02d', s), '/plots/Infect.png', sep = ''), width = 800, height = 600)
  par(family = 'LM Roman 10', mfrow = c(1, 1))
  plot(NA, xlim = c(0, N), ylim = c(0, max_y), main = '', xlab = 'Time', ylab = 'Number of individuals', xaxs = 'i', yaxs = 'i')
  polygon(c(dfI$t, rev(dfI$t)), c(dfI$L, rev(dfI$U)), col = rgb(1, 0, 0, alpha = 0.1), border = FALSE)
  lines(x = SIR$time, y = SIR$I, col = 2, lwd = 3)
  lines(x = SIR$time, y = dfI$M, col = 2, lty = 2, lwd = 3)
  abline(v = N_restricted, lty = 2)
  legend(x = 'topleft', legend = c('Infected', 'Est. Infected'), col = rep('red', 2), lty = c(1, 2), lwd = rep(2, 2), cex = 0.75, bg = 'white')
  dev.off()
}

