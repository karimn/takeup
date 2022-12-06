library(rstan)
stan_owen_t_code <-
  '
  functions {
    vector stan_owen_t(vector x, vector y) {
      return owens_t(x, y);
   }
  }
'
expose_stan_functions(stanc(model_code = stan_owen_t_code))