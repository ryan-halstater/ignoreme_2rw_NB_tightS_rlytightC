//Code for Lee-Carter model with Poisson or Negative Binomial error

data {
  int<lower = 1> J;                    // number of age categories
  int<lower = 1> T;                    // number of years
  array[J*T] int d;                          // vector of deaths
  vector[J* T] e;                      // vector of exposures
  int<lower = 1> Tfor;                  // number of forecast years
  int<lower = 0> Tval;                  // number of validation years
  array[J*Tval] int dval;                     // vector of deaths for validation
  vector[J* Tval] eval;                 // vector of exposures for validation
  int<lower=0,upper=1> family;          // family = 0 for Poisson, 1 for NB
}
transformed data {
  vector[J * T] input_offset = log(e);        // log exposures
  vector[J * Tval] offset2 = log(eval);     // log exposures for validation
  int<lower = 1> L;                     // size of prediction vector
  L=J*Tfor;
}
parameters {
  array[family > 0] real<lower=0> aux;       // neg. binomial inverse dispersion parameter
  vector[J] a;                          // alpha_x
  unit_vector[J] b1;                        // beta_x, strictly positive and sums to 1
  unit_vector[J] b2;                        // beta_x, strictly positive and sums to 1

  real c1;                           // drift term
  real c2;
  vector[T-1] ks1;                   // vector of kappa
  vector[T-1] ks2;                   // vector of kappa

  real<lower = 0> sigma1; 
  real<lower = 0> sigma2;   // standard deviation of the random walk
}
transformed parameters {        // This block defines a new vector where the first component is zero, this is required for identifiability of the Lee-Carter model. Otherwise, the chains will not converge.
  vector[T] k1;
  vector[T] k2;
  real phi = negative_infinity();       // neg. binomial dispersion parameter
  k1[1] = 0;
  k1[2:T]=ks1;
  k2[1] = 0;
  k2[2:T]=ks2;
  if (family > 0) phi = inv(aux[1]);
}

model {
  vector[J * T] mu;           //force of mortality
  int pos = 1;
  for (t in 1:T) for (x in 1:J) {
    mu[pos] = input_offset[pos]+ a[x] + b1[x] * k1[t] + b2[x] * k2[t];      // Predictor dynamics
    pos += 1;
  }

  target += normal_lpdf(ks1[1]|c1,sigma1);
  target += normal_lpdf(ks1[2:(T-1)]|c1+ks1[1:(T- 2)],sigma1);    // Random walk with drift prior
  
  target += normal_lpdf(ks2[1]|c2,sigma2);
  target += normal_lpdf(ks2[2:(T-1)]|c2+ks2[1:(T- 2)],sigma2);
  
  if (family ==0){
    target += poisson_log_lpmf(d |mu);                // Poisson log model
  }
  else {
    target +=neg_binomial_2_log_lpmf (d|mu,phi);      // Negative-Binomial log model
  }
  target += normal_lpdf(a|0,10);              // Prior on alpha_x
  target += normal_lpdf(c1|0,sqrt(0.001));  
  target += normal_lpdf(c2|0,sqrt(0.001));    // Prior on drift
  target += normal_lpdf(b1|0,sqrt(10)); 
  target += normal_lpdf(b2|0,sqrt(10));     // Prior on betas
  target += exponential_lpdf(sigma1 | 2);         // Exponential prior for sigma
  target += exponential_lpdf(sigma2 | 2); 
  if (family > 0) target += normal_lpdf(aux|0,1)- normal_lcdf(0 | 0, 1);
}

generated quantities {
  vector[Tfor] k_p1;
  vector[Tfor] k_p2;
  vector[L] mufor; // predicted death rates
  vector[J*T] log_lik; // log likelihood on data
  vector[J*Tval] log_lik2; // log likelihood on validation data
  int pos = 1;
  int pos2= 1;
  int pos3= 1;
  k_p1[1] = c1+k1[T]+sigma1 * normal_rng(0,1);
  k_p2[1] = c2+k2[T]+sigma2 * normal_rng(0,1);
  for (t in 2:Tfor) k_p1[t] = c1+k_p1[t - 1] + sigma1 * normal_rng(0,1);
  for (t in 2:Tfor) k_p2[t] = c2+k_p2[t - 1] + sigma2 * normal_rng(0,1);
  if (family==0){
    for (t in 1:Tfor) for (x in 1:J) {
    mufor[pos] = a[x] + b1[x] * k_p1[t] + b2[x] * k_p2[t];
    pos += 1;
  }
  mufor = exp(mufor);
  for (t in 1:T) for (x in 1:J) {
    log_lik[pos2] = poisson_log_lpmf (d[pos2] | input_offset[pos2]+ a[x] + b1[x] * k_p1[t] + b2[x] * k_p2[t]);
     pos2 += 1;
}

 for (t in 1:Tval) for (x in 1:J) {
    log_lik2[pos3] = poisson_log_lpmf(dval[pos3] | offset2[pos3]+  a[x] + b1[x] * k_p1[t] + b2[x] * k_p2[t]);
     pos3 += 1;
}
  }

    }


