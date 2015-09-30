##Use code to make .stan files to run specific models

#list names of all models in file
modelNames <- function(){
  readLines("stanModels.R")[grep("[\\s]*<-[\\s]* \\'",readLines("stanModels.R"),perl=T)]
}
#writes model to file to 'model.stan'
writeModel <- function(model = linearReg){
  writeLines(model,file("model.stan"))
}

modelNames()
writeModel(robustReg)

linearReg <- '

//linear regression
data { 
  int<lower=0> N;        // number of data items
  int<lower=0> K;        // number of predictors
  vector[N,K] x;         // predictor matrix
  vector[N] y;           // outcome vector
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  //error scale
}
model {
  y ~ normal(alpha + beta * x, sigma); //likelihood
}
'

robustReg <- '

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
'

logisticReg <- '

//logistic regression
data { 
  int<lower=0> N;        // number of data items
  int<lower=0> K;        // number of predictors
  vector[N,K] x;         // predictor matrix
  int<lower=0,upper=1> y[N];           // outcome vector
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  //error scale
}
model {
  y ~ bernoulli_logit(alpha + beta * x); //likelihood
}
'

probitReg <- '

//probit regression
data { 
  int<lower=0> N;        // number of data items
  int<lower=0> K;        // number of predictors
  vector[N,K] x;         // predictor matrix
  int<lower=0,upper=1> y[N];           // outcome vector
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  //error scale
}
model {
  y ~ bernoulli(Phi(alpha + beta * x)); //likelihood faster version y ~ bernoulli(Phi_approx(alpha + beta * x));
}
'

multiLogit <- '

data {
  int<lower=2> K;      // possible outcomes
  int<lower=0> N;      // data items
  int<lower=1> D;      //number of predictors
  int<lower=1,upper=K> y[N];
  vector[D] x[N];
}
transformed data{
  vector[D] zeros;
  zeros <- rep_vector(0,D);
}
parameters {
  matrix[K-1,D] beta_raw;
}
transformed parameters{
  matrix[K,D] beta;
  beta <- append_col(beta_raw,zeros)
}
model {
  to_vector(beta) ~ normal(0,5);
  for (n in 1:N)
    y[n] ~ catagorical_logit(beta * x[n]); // y[n] ~ catagorical(softmax(beta * x[n]));
}

'

OrderedLogistic <- '

data {
  int<lower=2> K;      // possible outcomes
  int<lower=0> N;      // data items
  int<lower=1> D;      //number of predictors
  int<lower=1,upper=K> y[N];
  row_vector[D] x[N];
}
parameters {
  vector[D] beta;
  ordered[K-1] c;
}
model {
  for (n in 1:N)
    y[n] ~ ordered_logistic(x[n] * beta, c);
}

'

OrderedProbit<- '

data {
  int<lower=2> K;      // possible outcomes
  int<lower=0> N;      // data items
  int<lower=1> D;      //number of predictors
  int<lower=1,upper=K> y[N];
  row_vector[D] x[N];
}
parameters {
  vector[D] beta;
  ordered[K-1] c;
}
model {
  vector[K] theta;
  for (n in 1:N) {
    real eta;
    eta <- x[n] * beta;
    theta[1] <- 1 - Phi(eta - c[1]);
    for (k in 2:(K-1))
      theta[k] <- Phi(eta - c[k-1]) = Phi(eta - c[k]);
    theta[K] <- Phi(eta - c[K-1]);
    y[n] ~ catagorical(theta);
  }
}

'

hierLogReg <- '

data {
  int<lower=1> D;
  int<lower=0> N;
  int<lower=1> L;
  int<lower=0,upper=1> y[N];
  int<lower=1,upper=L> ll[N];
  row_vector[D] x[N];
}
parameters {
  real mu[D];
  real<lower=0> sigma[D];
  vector[D] beta[L];
}
model {
  mu ~ normal(0,100);
  for (l in 1:L)
    beta[l] ~ normal(mu,sigma);
  for (n in 1:N)
    y[n] ~ bernoulli_logit(x[n] * beta[ll[n]]);   // y[n] ~ bernoulli(inv_logit(x[n] * beta[ll[n]]));
}

'

multivariateReg <- '

data {
  int<lower=0> N;               //num individuals
  int<lower=1> K;               //num ind predictors
  int<lower=1> J;               // num groups
  int<lower=1> L;               // num group predictors
  int<lower=1,upper=J> jj[N];   // group for ind
  matrix[N,K] x;                // indiv predictors
  row_vector[L] u[J];           // group preditctors
  vector[N] y;                  // outcomes
}
parameters {
  corr_matrix[K] Omega;         // prior correlation
  vector<lower=0>[K] tau;       // prior scale
  matrix[L,K] gamma;            // group coeffs
  vector[K] beta[J];            // indiv coeffs by group
  real<lower=0> sigma;          // prediction error scale
}
model {
  tau ~ cauchy(0,2.5);
  Omega ~ lkj_corr(2);
  to_vector(gamma) ~ normal(0,5);
  {
    row_vector[K] u_gamma[J];
    for (j in 1:J)
      u_gamma[j] <- u[j] * gamma;
    beta ~ multi_normal(u_gamma, quad_form_diag(Omega, tau));
  }
  {
    vector[N] x_beta_jj;
    for (n in 1:N)
      x_beta_jj[n] <- x[n] * beta[jj[n]];
    y ~ normal(x_beta_jj, sigma);
  }
}

'

