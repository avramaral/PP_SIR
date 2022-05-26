source('libraries.R')
source('utils.R')
source('data_generation.R')
source('intensities.R')
source('SIR_estimation.R')
source('spatioTemporal_estimation.R')
source('error_analysis.R')

source('scenarios.R')

start    <- 1
Terminal <- 100
delta    <- 1

mult_factor <- 1

nu     <- 1
scale  <- 0.05
sig_2  <- 0.2
mu     <- -1 * sig_2 / 2 
sd     <- sqrt(0.1)
sd_AR1 <- sd
a      <- 0.5

prop_class <- readRDS(file = 'data/p_3_classes.rds')
C <- readRDS(file = 'data/m_3_classes.rds')
C <- t(C)

I0 <- rep(x = 1, times = nrow(C))

S <- length(scenarios)