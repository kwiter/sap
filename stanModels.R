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



