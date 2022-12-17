library(rstan)
if (stan_version() == "2.21.0") {
stan_owen_t_code <-
  '
  functions {
    real stan_owen_t(real x, real y) {
      return owens_t(x, y);
   }
  }
'
} else {
stan_owen_t_code <-
  '
  functions {
    vector stan_owen_t(vector x, vector y) {
      return owens_t(x, y);
   }
  }
'
}
expose_stan_functions(stanc(model_code = stan_owen_t_code))