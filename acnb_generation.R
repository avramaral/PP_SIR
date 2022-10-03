acnb <- function (gamma) {
  
  set.seed(666)
  
  C_star <- C
  C_star[1, 2] <- mean(C[1, 2], C[2, 1]) 
  C_star[2, 1] <- mean(C[1, 2], C[2, 1]) 
  C_star[1, 3] <- mean(C[1, 3], C[3, 1]) 
  C_star[3, 1] <- mean(C[1, 3], C[3, 1]) 
  C_star[2, 3] <- mean(C[2, 3], C[3, 2]) 
  C_star[3, 2] <- mean(C[2, 3], C[3, 2]) 
  
  Terminal <- 100
  
  omega <- 2 * pi / 100
  
  trend <- ((1:100) / Terminal)
  sinCV <- cos(omega * (1:100))
  cosCV <- sin(omega * (1:100))
  
  coeff <- c(0, 1, -1, 1)
  alpha <- c(0.01, 0.4)
  
  overdisp <- 0.05 
  size <- 1 / overdisp
  
  inf1 <- matrix(data = 0, nrow = Terminal, ncol = 1)
  lam1 <- matrix(data = 0, nrow = Terminal, ncol = 1)
  inf2 <- matrix(data = 0, nrow = Terminal, ncol = 1)
  lam2 <- matrix(data = 0, nrow = Terminal, ncol = 1)
  inf3 <- matrix(data = 0, nrow = Terminal, ncol = 1)
  lam3 <- matrix(data = 0, nrow = Terminal, ncol = 1)
  inf1[1, ] <- 1
  inf2[1, ] <- 1
  inf3[1, ] <- 1
  
  for (t in 2:Terminal) {
    k <- 1
    lam1[t, 1] <- (prop_class[k] * exp(coeff[1] + trend[t] * coeff[2] + sinCV[t] * coeff[3] + cosCV[t] * coeff[4])) + 
      (exp(alpha[1]) * inf1[(t - 1), 1]) + 
      (exp(alpha[2]) * (C_star[k, 1] * inf1[(t - 1), 1] / (N_population * prop_class[1]) + C_star[k, 2] * inf2[(t - 1), 1] / (N_population * prop_class[2]) + C_star[k, 3] * inf3[(t - 1), 1] / (N_population * prop_class[3])))
    inf1[t, 1] <- rnbinom(n = 1, size = size, prob = (size / (size + lam1[t, 1])))
    
    k <- 2           
    lam2[t, 1] <- (prop_class[k] * exp(coeff[1] + trend[t] * coeff[2] + sinCV[t] * coeff[3] + cosCV[t] * coeff[4])) + 
      (exp(alpha[1]) * inf2[(t - 1), 1]) + 
      (exp(alpha[2]) * (C_star[k, 1] * inf1[(t - 1), 1] / (N_population * prop_class[1]) + C_star[k, 2] * inf2[(t - 1), 1] / (N_population * prop_class[2]) + C_star[k, 3] * inf3[(t - 1), 1] / (N_population * prop_class[3])))                   
    inf2[t, 1] <- rnbinom(n = 1, size = size, prob = (size / (size + lam2[t, 1])))
    
    k <- 3           
    lam3[t, 1] <- (prop_class[k] * exp(coeff[1] + trend[t] * coeff[2] + sinCV[t] * coeff[3] + cosCV[t] * coeff[4])) + 
      (exp(alpha[1]) * inf3[(t - 1), 1]) + 
      (exp(alpha[2]) * (C_star[k, 1] * inf1[(t - 1), 1] / (N_population * prop_class[1]) + C_star[k, 2] * inf2[(t - 1), 1] / (N_population * prop_class[2]) + C_star[k, 3] * inf3[(t - 1), 1] / (N_population * prop_class[3])))                   
    inf3[t, 1] <- rnbinom(n = 1, size = size, prob = (size / (size + lam3[t, 1])))
  }
  
  # par(mfrow = c(3, 1))
  # plot(c(inf1[, 1]), type = 'l', col = 'red')
  # plot(c(inf2[, 1]), type = 'l', col = 'red')
  # plot(c(inf3[, 1]), type = 'l', col = 'red')
  
  S.V1 <- c()
  S.V2 <- c()
  S.V3 <- c()
  I.I1 <- c()
  I.I2 <- c()
  I.I3 <- c()
  R.R1 <- c()
  R.R2 <- c()
  R.R3 <- c()
  recover_days <- 1 / gamma
  for (t in 1:Terminal) {
    I.I1 <- c(I.I1, sum(c(inf1)[max(1, (t - ceiling(recover_days))):t]))
    I.I2 <- c(I.I2, sum(c(inf2)[max(1, (t - ceiling(recover_days))):t]))
    I.I3 <- c(I.I3, sum(c(inf3)[max(1, (t - ceiling(recover_days))):t]))
    R.R1 <- c(R.R1, ifelse(test = (t - ceiling(recover_days)) < 1, yes = 0, no = sum(c(inf1)[1:(t - recover_days)])))
    R.R2 <- c(R.R2, ifelse(test = (t - ceiling(recover_days)) < 1, yes = 0, no = sum(c(inf2)[1:(t - recover_days)])))
    R.R3 <- c(R.R3, ifelse(test = (t - ceiling(recover_days)) < 1, yes = 0, no = sum(c(inf3)[1:(t - recover_days)])))
  }
  S.V1 <- (N_population * prop_class[1]) - I.I1 - R.R1
  S.V2 <- (N_population * prop_class[2]) - I.I2 - R.R2
  S.V3 <- (N_population * prop_class[3]) - I.I3 - R.R3
  res <- cbind(time = 1:Terminal, S = (S.V1 + S.V2 + S.V3), I = (I.I1 + I.I2 + I.I3), R = (R.R1 + R.R2 + R.R3))
  
  SIR <- list()
  SIR$SIR <- as.data.frame(res)
  SIR$SIR_sep <- as.data.frame(cbind(time = 1:Terminal, S.V1 = S.V1, S.V2 = S.V2, S.V3 = S.V3, I.I1 = I.I1, I.I2 = I.I2, I.I3 = I.I3, R.R1 = R.R1, R.R2 = R.R2, R.R3 = R.R3))
  
  SIR
}
