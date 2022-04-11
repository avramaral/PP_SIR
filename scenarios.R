
inv_phi_EP <- 0.005
alpha_EP   <- 0.05
beta_EP    <- 0.000025
gamma_EP   <- 0.2

inv_phi_FC <- 0.005
alpha_FC   <- 0.05
beta_FC    <- 0.000009 
gamma_FC   <- 0.1 

scenarios <- list()

scenarios[[1]] <- list(inv_phi = inv_phi_EP,
                       alpha = alpha_EP,
                       beta = beta_EP,
                       gamma = gamma_EP,
                       null_model = TRUE,
                       AR_include = FALSE)

scenarios[[2]] <- list(inv_phi = inv_phi_EP,
                       alpha = alpha_EP,
                       beta = beta_EP,
                       gamma = gamma_EP,
                       null_model = FALSE,
                       AR_include = FALSE)

scenarios[[3]] <- list(inv_phi = inv_phi_EP,
                       alpha = alpha_EP,
                       beta = beta_EP,
                       gamma = gamma_EP,
                       null_model = TRUE,
                       AR_include = TRUE)

scenarios[[4]] <- list(inv_phi = inv_phi_EP,
                       alpha = alpha_EP,
                       beta = beta_EP,
                       gamma = gamma_EP,
                       null_model = FALSE,
                       AR_include = TRUE)

scenarios[[5]] <- list(inv_phi = inv_phi_FC,
                       alpha = alpha_FC,
                       beta = beta_FC,
                       gamma = gamma_FC,
                       null_model = TRUE,
                       AR_include = FALSE)

scenarios[[6]] <- list(inv_phi = inv_phi_FC,
                       alpha = alpha_FC,
                       beta = beta_FC,
                       gamma = gamma_FC,
                       null_model = FALSE,
                       AR_include = FALSE)

scenarios[[7]] <- list(inv_phi = inv_phi_FC,
                       alpha = alpha_FC,
                       beta = beta_FC,
                       gamma = gamma_FC,
                       null_model = TRUE,
                       AR_include = TRUE)

scenarios[[8]] <- list(inv_phi = inv_phi_FC,
                       alpha = alpha_FC,
                       beta = beta_FC,
                       gamma = gamma_FC,
                       null_model = FALSE,
                       AR_include = TRUE)
