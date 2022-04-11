functions {
  real[] SIR (int N_population, real delta, real b_noise, real[] values, real[] pars) {

    real new_I = values[2] + (pars[2] / N_population) * values[1] * values[2] * delta - pars[3] * values[2] * delta + pars[1] * values[2] * b_noise; 
    real new_R = values[3] + pars[3] * values[2] * delta; 
    real new_S = N_population - new_I - new_R;
    
    return {new_S, new_I, new_R};
  }
  
  real[, ] solve_SDE (int N, int N_population, real delta, real[] b_noise, real[, ] values, real[] pars) {
    
    real result[N, 3];
    
    result[1, ] = values[1, ];
    for (n in 2:N) {
      result[n, ] = SIR(N_population, delta, b_noise[n], result[(n - 1), ], pars);
    }
    return result;
  }
}

data {
  int<lower=1> N;
  int<lower=1> N_restricted;
  int<lower=1> N_population;
  real delta;
  real t[N_restricted];
  real Y[N_restricted, 3];
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> gamma;
}

transformed parameters {
  real pars[3];
  
  pars[1] = alpha;
  pars[2] = beta;
  pars[3] = gamma;
}

model {
  // Priors
  alpha ~ normal(0.5, 1.);
  beta  ~ normal(0.5, 1.);
  gamma ~ normal(0.5, 1.);
  
  for (n in 2:N_restricted) {
    // Conditioning on the previous step
    Y[n, 2] ~ normal(Y[(n - 1), 2] + ((beta / N_population) * Y[(n - 1), 1] * Y[(n - 1), 2] - gamma * Y[(n - 1), 2]) * delta, sqrt((alpha * Y[(n - 1), 2]) ^ 2 * delta));
  }
}

generated quantities {
  real Y_hat [N, 3];
  real b_noise[N];
  
  b_noise = to_array_1d(multi_normal_rng(rep_vector(0, N), diag_matrix(rep_vector(delta, N))));
  Y_hat = solve_SDE(N, N_population, delta, b_noise, Y, pars);
}
