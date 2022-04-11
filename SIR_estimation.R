SIR_SDE <- function (start, Terminal, delta, N, N_restricted, N_population, SIR, n_chains = 4, iter = 4e3, warmup = 2e3) {
  
  ts <- seq(from = start, to = Terminal, by = delta)
  data_stan <- list(N = N,
                    N_restricted = N_restricted,
                    N_population = N_population,
                    delta = delta,
                    t = ts[1:N_restricted],
                    Y = as.matrix(SIR[1:N_restricted, 2:4]))
  
  initial_values <- function(chain_id = 1) { 
    list(
      alpha   = runif(1, 0, 1), 
      beta    = runif(1, 0, 2) / N_population, 
      gamma   = runif(1, 0, 1)
    ) 
  }
  
  fit <- stan(file = 'model_SIR_SDE.stan',
              init = lapply(1:n_chains, function(id) initial_values(chain_id = id)),
              data = data_stan, 
              chains = n_chains, 
              iter = iter, 
              warmup = warmup,
              control = list(adapt_delta = 0.99))
  
  fit
}

SIR_NB <- function (start, Terminal, N, N_restricted, N_population, SIR, n_chains = 4, iter = 4e3, warmup = 2e3) {
  
  ts <- seq(from = start, to = Terminal, by = delta)
  data_stan <- list(N = (N_restricted - 1),
                    Y0 = unname(unlist(SIR[1, 2:4])),
                    t0 = 0,
                    ts = ts[1:(N_restricted - 1)],
                    N_population = N_population,
                    cases = SIR$I[2:N_restricted])
  
  fit <- stan(file = 'model_SIR.stan',
              data = data_stan, 
              chains = n_chains, 
              iter = iter, 
              warmup = warmup,
              control = list(adapt_delta = 0.99))
  
  fit
}


mixing_assessment <- function (fit, par_names = c('alpha', 'beta', 'gamma'), plotting = TRUE, scaled_beta = FALSE, N_population = NULL) {
  est_par1 <- rstan::extract(fit, pars = c(par_names[1]))[[1]]
  est_par2 <- rstan::extract(fit, pars = c(par_names[2]))[[1]]
  est_par3 <- rstan::extract(fit, pars = c(par_names[3]))[[1]]
  
  if (scaled_beta) { est_par2 <- est_par2 / N_population }
  
  n_chains <- length(fit@stan_args)
  iter <- fit@stan_args[[1]]$iter
  warmup <- fit@stan_args[[1]]$warmup
  
  n_samples <- iter - warmup
  
  if (plotting) {
    par(family = 'LM Roman 10', mfrow = c(3, 1))
    plot(NA, xlab = 'Iteration', ylab = par_names[1], xlim = c(0, (n_samples)), ylim = c(min(est_par1), max(est_par1)), xaxs = 'i', yaxs = 'i')
    for (i in 1:n_chains) { lines(est_par1[((i - 1) * n_samples + 1):(i * n_samples)], col = (i + 1)) }
    plot(NA, xlab = 'Iteration', ylab = par_names[2] , xlim = c(0, (n_samples)), ylim = c(min(est_par2),  max(est_par2) ), xaxs = 'i', yaxs = 'i')
    for (i in 1:n_chains) { lines(est_par2[((i - 1) * n_samples + 1):(i * n_samples)], col = (i + 1)) }
    plot(NA, xlab = 'Iteration', ylab = par_names[3], xlim = c(0, (n_samples)), ylim = c(min(est_par3), max(est_par3)), xaxs = 'i', yaxs = 'i')
    for (i in 1:n_chains) { lines(est_par3[((i - 1) * n_samples + 1):(i * n_samples)], col = (i + 1)) }
    legend(x = 'bottomright', legend = sprintf('Chain %s', seq(1:n_chains)), col = (2:(n_chains + 1)), lty = rep(x = 1, times = n_chains), cex = 0.75)
    par(mfrow = c(1, 1))
  }
  
  list(est_par1 = est_par1, est_par2 = est_par2, est_par3 = est_par3)
}

parameters_summary <- function (parameters, par_names = c('alpha', 'beta', 'gamma'), probs = c(0.025, 0.975), round_dig = 6) {
  est_par1_q <- round(quantile(x = parameters$est_par1, probs = probs), round_dig)
  est_par1_m <- round(mean(parameters$est_par1), round_dig)
  est_par2_q <- round(quantile(x = parameters$est_par2, probs = probs), round_dig)
  est_par2_m <- round(mean(parameters$est_par2), round_dig)
  est_par3_q <- round(quantile(x = parameters$est_par3, probs = probs), round_dig)
  est_par3_m <- round(mean(parameters$est_par3), round_dig)
  parSummary <- as.data.frame(rbind(c(est_par1_q, est_par1_m), c(est_par2_q, est_par2_m), c(est_par3_q, est_par3_m)))
  rownames(parSummary) <- par_names
  colnames(parSummary) <- c(paste(probs[1] * 100, '%', sep = ''), paste(probs[2] * 100, '%', sep = ''), 'Mean')
  parSummary
}

generate_Y_hat <- function (parameters, N_population) {
  
  phi   <- parameters$est_par1
  beta  <- parameters$est_par2
  gamma <- parameters$est_par3 
  
  progressbar <- txtProgressBar(min = 1, max = length(beta), initial = 1) 
  
  m <- array(data = NA, dim = c(length(beta), 100, 3))
  for (i in 1:length(beta)) {
    m_partial <- as.matrix(generate_SIR(N_population = N_population, beta = beta[i], gamma = gamma[i])[2:4])
    m[i, , ] <- cbind(m_partial[, 1], rnbinom(n = 100, size = phi[i], mu = m_partial[, 2]), m_partial[, 3])
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  m
}

