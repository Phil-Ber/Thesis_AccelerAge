# Source: https://github.com/marije-sluiskes/fitting-accelerage-framework-in-r/blob/main/src/weibull_draw.R
#################################################################
# This script contains functions to draw event times from a Weibull distribution
# which can be parameterized as an AFT and PH model at the same time.
# h_0(t) = lambda * nu * t^(nu - 1)
################################################################

# h(t) = l*v*t^(v-1) (Bender et al. parametrization, https://epub.ub.uni-muenchen.de/1716/1/paper_338.pdf)
weib_baseline_surv <- function(t, lambda, nu){
  return(exp(-lambda*t^nu))
}

inverse_weib_baseline_surv <- function(s, lambda, nu){
  return((-(1/lambda) * log(s))^(1/nu))
}

rweibull_custom <- function(n, lambda, nu, linpred){
  u <- runif(n)
  lambda <- lambda * exp(linpred)
  return(( - (log(u) / lambda) )^(1/nu))
}

rcweibull <- function(n, lambda, nu, linpred, s){
  u <- runif(n)
  lambda <- lambda * exp(linpred)
  return( (s^nu - log(u) / lambda )^(1/nu) )
}

# returns weib_baseline_surv ^ exp(linpred)
weib_PH_surv <- function(t, lambda, nu, linpred){
  weib_baseline_surv(t, lambda, nu) ^ exp(linpred)
}

# Weibull survival with covariate-adjusted scale
S_weibull_cov <- function(t, shape, scale_adj) {
  exp(-(t / scale_adj)^shape)
}


rweib_c <- function(n, shape, scale, linpred) {
  u = runif(n, 0, 150) * exp(linpred) # Point sampled from distr
  return(dweibull(u, shape, scale))
}