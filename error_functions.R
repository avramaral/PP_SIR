extract_infectious_gen <- function (s, R = TRUE) {
  
  SIR <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/SIR.rds', sep = ''))
  intensities <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/intensities.rds', sep = ''))
  
  infectious <- list()
  for (k in 1:n_classes) {
    infectious[[k]] <- list()
    for (t in 1:Terminal) {
      partial <- values(intensities[[k]][[t]]) * prod(res(intensities[[k]][[t]]))
      if (R) {
        infectious[[k]][[t]] <- floor(partial)  
      } else {
        infectious[[k]][[t]] <- partial  
      }
    }
  }
  
  infectious
}

extract_infectious_fit <- function (s, R = TRUE) {
  
  processed_result <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/processed_result.rds', sep = ''))
  
  proc_result <- list()
  for (k in 1:n_classes) {
    proc_result[[k]] <- list()
    for (t in 1:Terminal) {
      proc_result[[k]][[t]] <- processed_result[[k]]$result_r_mean[[t]]
      partial <- values(proc_result[[k]][[t]])  * prod(res(proc_result[[k]][[t]]))
      if (R) {
        proc_result[[k]][[t]] <- floor(partial)  
      } else {
        proc_result[[k]][[t]] <- partial  
      }
    }
  }
  
  proc_result
}

extract_infectious_fit_quant <- function (s, R = TRUE) {
  
  processed_result <- readRDS(file = paste('output/', sprintf('%02d', s), '/rds/processed_result.rds', sep = ''))
  
  proc_result_0025 <- list()
  proc_result_mean <- list()
  proc_result_0975 <- list()
  for (k in 1:n_classes) {
    proc_result_0025[[k]] <- list()
    proc_result_mean[[k]] <- list()
    proc_result_0975[[k]] <- list()
    for (t in 1:Terminal) {
      proc_result_0025[[k]][[t]] <- processed_result[[k]]$result_r_0250[[t]]
      proc_result_mean[[k]][[t]] <- processed_result[[k]]$result_r_mean[[t]]
      proc_result_0975[[k]][[t]] <- processed_result[[k]]$result_r_0975[[t]]
      partial_0025 <- values(proc_result_0025[[k]][[t]]) * prod(res(proc_result_0025[[k]][[t]]))
      partial_mean <- values(proc_result_mean[[k]][[t]]) * prod(res(proc_result_mean[[k]][[t]]))
      partial_0975 <- values(proc_result_0975[[k]][[t]]) * prod(res(proc_result_0975[[k]][[t]]))
      if (R) {
        proc_result_0025[[k]][[t]] <- floor(partial_0025)  
        proc_result_mean[[k]][[t]] <- floor(partial_mean)  
        proc_result_0975[[k]][[t]] <- floor(partial_0975)
      } else {
        proc_result_0025[[k]][[t]] <- partial_0025  
        proc_result_mean[[k]][[t]] <- partial_mean  
        proc_result_0975[[k]][[t]] <- partial_0975
      }
    }
  }
  
  list(proc_result_0025 = proc_result_0025, proc_result_mean = proc_result_mean, proc_result_0975 = proc_result_0975)
}

compute_RMSE <- function (x, y) { # Root Mean Square Error
  N <- length(x)
  r <- (x - y) ** 2 
  sqrt(sum(r) / N)
}

compute_MAPE <- function (x, y) {
  N <- length(x)
  r <- (abs((x - y) / (x)))
  for (i in 1:length(r)) { if (is.nan(r[i])) { r[i] <- 0 } } 
  sum(r) / N
}

compute_MAAPE <- function (x, y) {
  N <- length(x)
  r <- atan(abs((x - y) / (x)))
  for (i in 1:length(r)) { if (is.nan(r[i])) { r[i] <- 0 } }
  sum(r) / N
}

compute_APE <- function (x, y) {
  N <- length(x)
  r <- (abs((x - y) / (x)))
  for (i in 1:length(r)) { if (is.nan(r[i])) { r[i] <- 0 } } 
  r
}

compute_AAPE <- function (x, y) {
  N <- length(x)
  r <- atan(abs((x - y) / (x)))
  for (i in 1:length(r)) { if (is.nan(r[i])) { r[i] <- 0 } }
  r
}

indicator_coverage <- function (inf_gen, inf_fit_0025, inf_fit_0975) {
  
  L <- length(inf_gen)
  
  cum <- 0
  for (l in 1:L) {
    if ((inf_gen[l] >= inf_fit_0025[l]) & (inf_gen[l] <= inf_fit_0975[l])) {
      cum <- cum + 1
    } 
  }
  
  (100 * (1 / L) * cum)
}
