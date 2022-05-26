functions {
  real[] SIR(real t, real[] Y, real[] pars, real[] x_r, int[] x_i) {

      int n_classes = x_i[1];
      vector[n_classes] N_population = to_vector(x_r[1:n_classes]);
      matrix[n_classes, n_classes] C = to_matrix(x_r[(n_classes + 1):(n_classes + (n_classes * n_classes))], n_classes, n_classes, 0);
      
      vector[n_classes] S = to_vector(Y[1:n_classes]);
      vector[n_classes] I = to_vector(Y[(n_classes + 1):(2 * n_classes)]);
      vector[n_classes] R = to_vector(Y[((2 * n_classes) + 1):(3 * n_classes)]);
      
      real beta  = pars[1];
      real gamma = pars[2];
      
      vector[n_classes] dS_dt = (to_vector((-1 * beta) * S)) .* to_vector(C * to_matrix((I ./ N_population)));
      vector[n_classes] dI_dt = (-1 * dS_dt) - (to_vector(gamma * I));
      vector[n_classes] dR_dt = to_vector(gamma * I);
      
      return append_array(append_array(to_array_1d(dS_dt), to_array_1d(dI_dt)), to_array_1d(dR_dt));
  }
}

data {
  int<lower=1> N;
  int<lower=1> n_classes;
  real Y0[3 * n_classes];
  real t0;
  real ts[N];
  int N_population[n_classes];
  int cases[N, n_classes];
  real C[n_classes, n_classes];
}

transformed data {
  int x_r_length = n_classes + (n_classes * n_classes);
  real x_r[x_r_length];
  int x_i[1];
  int count = 1 + n_classes;
  
  x_i[1] = n_classes;
  // Population size for each group
  for (i in 1:n_classes) {
    x_r[i] = N_population[i];
  }
  // Contact matrix
  for (i in 1:n_classes) {
    for (j in 1:n_classes) {
      x_r[count] = C[i, j]; 
      count = count + 1;
    }
  }
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi; 
}

transformed parameters{
  real Y[N, 3 * n_classes];
  {
    real pars[2];
    pars[1] = beta;
    pars[2] = gamma;

    Y = integrate_ode_rk45(SIR, Y0, t0, ts, pars, x_r, x_i);
  }
}

model {
  // Priors
  beta ~ normal(0.5, 1.00);
  gamma ~ normal(0.5, 1.0);
  phi ~ normal(1, 100);
  
  // Sampling distribution
  for (i in 1:n_classes) {
    for (n in 1:N) {
      cases[n, i] ~ neg_binomial_2(col(to_matrix(Y), (n_classes + i))[n], phi);
    }
  }
}
