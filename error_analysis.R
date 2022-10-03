source('header.R')
source('error_functions.R')

library('tidyverse')
library('geofacet')
library('viridis')
require('gridExtra')

n_classes <- 3

plot_error <- function(area, map, s, t) {
  mx_legend <- 1.6
  n_classes <- 3
  
  m <- ggmap(map) +
    geom_tile(data = gplot_data(area), aes(x = x, y = y, fill = value), alpha = 0.5) +
    guides(color = 'none') +
    scale_fill_gradientn(colors = rainbow(n = 100, start = 0.1, end = 0.9), name = 'MAAPE', limits = c(0, mx_legend), breaks = seq(0, mx_legend, 0.4)) +
    labs(x = 'Longitude', y = 'Latitude', title = paste('Scenario ', sprintf('%02d', s), '. Time window ', sprintf('%03d', t), '.', sep = '')) +
    theme(text = element_text(family = 'LM Roman 10'), legend.key.width = unit(0.5, 'cm'), legend.key.height = unit(0.615, 'cm'),
          panel.background = element_rect(fill = 'transparent', color = NA), plot.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
  
  m    
}

ERRORS <- list()
for (s in 1:24) {
  print(s)
  
  inf_gen <- extract_infectious_gen(s = s, R = T)
  inf_fit <- extract_infectious_fit(s = s, R = T)
  
  errors <- list()
  for (k in 1:n_classes) {
    for (t in 1:Terminal) {
      if (t == 1) {
        errors[[k]] <- compute_AAPE(inf_gen[[k]][[t]], inf_fit[[k]][[t]])
      } else {
        errors[[k]] <- c(errors[[k]], compute_AAPE(inf_gen[[k]][[t]], inf_fit[[k]][[t]]))
      }
    }
  }
  
  ERRORS[[s]] <- errors
}

data <- data.frame(scenario = rep(1:S, each = n_classes * Terminal * 300),
                   class = rep(rep(1:n_classes, each = Terminal * 300), times = S),
                   time = rep(rep(1:Terminal, each = 300), times = n_classes * S) - 1,
                   region = rep(1:300, times = Terminal * n_classes * S),
                   error = unlist(ERRORS))

for (s in seq(1, 23, 2)) {
  print(s)
  custLab <- function (x){ c('0-19', '20-59', '60+') } 
  legn_marks <- c(0, 0.4, 0.8, 1.2, 1.6)
  cell_marks <- c(100, 200, 300)
  time_marks <- c(0, 19, 39, 59, 79, 99)
  a <- ggplot(data = data[(data$scenario == s), ], mapping = aes(x = time, y = region, fill = error)) +
    geom_tile() +
    geom_vline(xintercept = 49, color = 'magenta', size = 1.5) +
    scale_x_continuous(expand = c(0, 0), breaks = time_marks, labels = time_marks) +
    scale_y_continuous(expand = c(0, 0), breaks = cell_marks, labels = cell_marks) +
    scale_fill_viridis(name = "AAPE", option = 'turbo', limits = c(0, 1.6), breaks = legn_marks) +
    labs(x = '', y = 'Cells') +
    facet_wrap(~ class, labeller = as_labeller(custLab)) + 
    theme_bw() +
    theme(text = element_text(family = 'LM Roman 10', size = 16), panel.spacing.x = unit(1.25, 'lines'),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(0.80, 'in'),
          plot.margin = unit(c(0.25, 0.1, -0.5, 0.1), 'cm')) # top, right, bottom, left
  b <- ggplot(data = data[(data$scenario == (s + 1)), ], mapping = aes(x = time, y = region, fill = error)) +
    geom_tile() +
    geom_vline(xintercept = 49, color = 'magenta', size = 1.5) +
    scale_x_continuous(expand = c(0, 0), breaks = time_marks, labels = time_marks) +
    scale_y_continuous(expand = c(0, 0), breaks = cell_marks, labels = cell_marks) +
    scale_fill_viridis(name = "AAPE", option = 'turbo', limits = c(0, 1.6), breaks = legn_marks) +
    labs(x = 'Time', y = 'Cells') +
    facet_wrap(~ class, labeller = as_labeller(custLab)) + 
    theme_bw() +
    theme(text = element_text(family = 'LM Roman 10', size = 16), panel.spacing.x = unit(1.25, 'lines'),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(0.765, 'in'),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
  png(filename = paste("result_processing/aape_", sprintf('%02d', s), "-", sprintf('%02d', (s + 1)), ".png", sep = ""), width = 1200, height = 1200, units = "px", res = 120)
  grid.arrange(a, b, nrow = 2, padding = 0)
  dev.off()
}
