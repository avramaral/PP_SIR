functions {
  real[] SIR(real t, real[] Y, real[] pars, real[] x_r, int[] x_i) {

      real S = Y[1];
      real I = Y[2];
      real R = Y[3];
      real N_population = x_i[1];
      
      real beta  = pars[1];
      real gamma = pars[2];
      
      real dS_dt = -1 * (beta / N_population) * I * S;
      real dI_dt = (beta / N_population) * I * S - gamma * I;
      real dR_dt = gamma * I;
      
      return { dS_dt, dI_dt, dR_dt };
  }
}

data {
  int<lower=1> N;
  real Y0[3];
  real t0;
  real ts[N];
  int N_population;
  int cases[N];
}

transformed data {
  real x_r[0];
  int x_i[1] = { N_population };
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi; 
}

transformed parameters{
  real Y[N, 3];
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
  cases ~ neg_binomial_2(col(to_matrix(Y), 2), phi);
}
