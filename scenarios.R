
inv_phi_EP <- 0.010
beta_EP    <- 0.040
gamma_EP   <- 0.200

inv_phi_FC <- 0.0100
beta_FC    <- 0.0175 
gamma_FC   <- 0.1000 

scenarios <- list()

#######################################################################################

scenarios[[1]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'IID',
                        null_model = TRUE,
                        AR_include = FALSE)

scenarios[[2]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'IID',
                        null_model = FALSE,
                        AR_include = FALSE)

scenarios[[3]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'IID',
                        null_model = TRUE,
                        AR_include = TRUE)

scenarios[[4]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'IID',
                        null_model = FALSE,
                        AR_include = TRUE)

scenarios[[5]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'AR1',
                        null_model = TRUE,
                        AR_include = FALSE)

scenarios[[6]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'AR1',
                        null_model = FALSE,
                        AR_include = FALSE)

scenarios[[7]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'AR1',
                        null_model = TRUE,
                        AR_include = TRUE)

scenarios[[8]]  <- list(inv_phi = inv_phi_EP,
                        beta = beta_EP,
                        gamma = gamma_EP,
                        model_generation = 'AR1',
                        null_model = FALSE,
                        AR_include = TRUE)

#######################################################################################

scenarios[[9]]  <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'IID',
                        null_model = TRUE,
                        AR_include = FALSE)

scenarios[[10]] <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'IID',
                        null_model = FALSE,
                        AR_include = FALSE)

scenarios[[11]] <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'IID',
                        null_model = TRUE,
                        AR_include = TRUE)

scenarios[[12]] <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'IID',
                        null_model = FALSE,
                        AR_include = TRUE)

scenarios[[13]] <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'AR1',
                        null_model = TRUE,
                        AR_include = FALSE)

scenarios[[14]] <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'AR1',
                        null_model = FALSE,
                        AR_include = FALSE)

scenarios[[15]] <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'AR1',
                        null_model = TRUE,
                        AR_include = TRUE)

scenarios[[16]] <- list(inv_phi = inv_phi_FC,
                        beta = beta_FC,
                        gamma = gamma_FC,
                        model_generation = 'AR1',
                        null_model = FALSE,
                        AR_include = TRUE)

#######################################################################################

scenarios[[17]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'IID',
                         null_model = TRUE,
                         AR_include = FALSE)

scenarios[[18]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'IID',
                         null_model = FALSE,
                         AR_include = FALSE)

scenarios[[19]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'IID',
                         null_model = TRUE,
                         AR_include = TRUE)

scenarios[[20]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'IID',
                         null_model = FALSE,
                         AR_include = TRUE)

scenarios[[21]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'AR1',
                         null_model = TRUE,
                         AR_include = FALSE)

scenarios[[22]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'AR1',
                         null_model = FALSE,
                         AR_include = FALSE)

scenarios[[23]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'AR1',
                         null_model = TRUE,
                         AR_include = TRUE)

scenarios[[24]]  <- list(phi = 0.5,
                         gamma = 0.1,
                         model_generation = 'AR1',
                         null_model = FALSE,
                         AR_include = TRUE)
