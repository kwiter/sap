//Robust regression

data { 
  int<lower=0> N;        // number of data items
  int<lower=0> K;        // number of predictors
  vector[N,K] x;         // predictor matrix
  vector[N] y;           // outcome vector
  real<lower=0> nu;      // degress of freedom for student T
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  //error scale
}
model {
  y ~ student_t(nu,alpha + beta * x, sigma); //likelihood
}