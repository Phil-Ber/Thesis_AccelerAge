library(survival)
library(scam)
library(glue)
library(psbcGroup)
library(tidyverse)
# Load custom functions
source("R/bioage_estimate_median.R")
source("R/gompertz_draw.R")
source("R/weibull_draw.R")

generate_population_lifetable = function(
    # Source: https://github.com/marije-sluiskes/fitting-accelerage-framework-in-r/blob/main/AccelerAge-framework-illustration.md
    ################################################################
    #### Create lifetable for population
    ################################################################
    N_pop = 1e5,            # Size of theoretical population
    M = 3,                  # Number of predictors
    betas = c(-1,0,1),      # Vector of coefficients (length should be M)
    a = exp(-9),                              # Gompertz baseline param
    b = 0.085,                                # Gompertz baseline param
    filename = "",          # Optional filename for saving the lifetable
    seed = 123,             # Seed for random number generation (default 123)
    force_recalc = F,       # Force recalculation even if file exists
    X                       # Matrix of samples
    
) {
  sigma <- 1/b # canonical parametrization
  tau <- a/b # canonical parametrization
  path = glue("pop_lifetable_{filename}.rds")
  if (!file.exists(path) || force_recalc) {
    # Ages
    set.seed(seed)
    # X = matrix( rnorm(N_pop*M,mean=0,sd=1), N_pop, M) # matrix of predictors
    # Generate X -> Beta
    
    linpred = (as.matrix(X) %*% as.matrix(betas))
    t = vector(length = N_pop)                                                     
    for (i in 1:N_pop){
      t_i = rgompertz_aft(1, sigma = sigma, tau = tau, linpred = linpred[i]) # vector of ages-of-death
      if (t_i > 150) { # If the age-of-death > 150, we assume 150
        t_i = 150
      }
      t[i] = t_i
    }
    
    print(glue("Range of death: {range(t)}"))
    
    plot(ecdf(t), xlim = c(0,100), main = "Cumulative ")
    
    lifetable_pop = as.data.frame(cbind(X, t))
    lifetable_pop = lifetable_pop[order(lifetable_pop$t),]
    
    mrl = vector(length = N_pop)
    for (j in 1:N_pop){
      mrl[j] = mean(lifetable_pop$t[j:N_pop]) - lifetable_pop$t[j]
    }
    
    lifetable_pop$mrl = mrl
    
    # smoothen 
    fitsmooth = scam(mrl ~ s(t, bs = "mpd"), data = lifetable_pop)
    xx = seq(0,max(lifetable_pop$t), by = 0.11)
    lt = as.data.frame(cbind(t = xx, mrl =  predict(fitsmooth, data.frame(t=xx))))
    
    #M  Montecarlo trials
    #p number of coefficents
    # estBetaBayes=matrix(NA,M,p)
    # estBetaEha=matrix(NA,M,p)
    # for( i in 1:M){
    #   set.seed(i)
    # }
    # Save lifetable to R object
    saveRDS(lt, path)
  } else {
    # If lifetable exists, read it
    lt = readRDS(path)
  }
  ggplot(lt, aes(x = t, y = mrl)) +
    geom_line() +
    labs(
      title = "Mean Residual Life Function",
      x = "Time (t)",
      y = "Mean Residual Life"
    ) +
    theme_minimal()
  return(lt)
}


create_dataset = function(
    # Source: https://github.com/marije-sluiskes/fitting-accelerage-framework-in-r/blob/main/AccelerAge-framework-illustration.md
    ################################################################
    #### Create data set 
    ################################################################
    M,                                        # number of predictors
    n_obs,                                    # number of *observed* subjects
    a,                              # Gompertz baseline param
    b,                                # Gompertz baseline param
    followup,                            # maximum follow‑up time
    G,                                    # number of groups (for grouped)
    gsize,                                # group size (for grouped)
    seednr,                             # seed for reproducibility
    betas,                         # Pop param
    lt,
    X_rho,
    X_plots,
    X_scale
) {
  
  sigma <- 1/b # canonical parametrization
  tau <- a/b # canonical parametrization
  # Check if betas are of length M, and is a numeric/float vector
  
  
  n_gen <- 5 * n_obs # 5 TIMES as many to ensure I generate enough, because for some T < C => not observed

  # Generate it in 
  result = generate_X(n = n_gen, p = M, g = G, rho = X_rho, rho_between = 0, seed = seednr, scale = X_scale, X_plots = X_plots)

  X = result$X # Extract X from the list
  
  cnames <- as.character(glue("x{1:M}"))
  colnames(X) <- cnames
  age_start <- runif(n_gen, 20, 80)
  linpred <- rowSums(sweep(X, 2, betas, "*"))
  
  # Get age of death
  age_death <- vector(length = n_gen)
  for (i in 1:n_gen){
    age_death[i] <- rgompertz_aft(1, sigma = sigma, linpred = linpred[i], tau = tau)
  }
  
  # Remove observations that are left-truncated
  valid_indices <- which(age_start < age_death)
  
  # Check if we have enough valid cases
  if (length(valid_indices) < n_obs) {
    warning(paste("Only", length(valid_indices), "valid cases where age_start < age_death, but", n_obs, "requested."))
    # Use all available valid indices
    indx_obs <- valid_indices
  } else {
    # Use the first n_obs valid indices
    indx_obs <- valid_indices[1:n_obs]
  }
  
  print(glue("Nr of valid indices: {length(valid_indices)}"))
  
  # Create the data frame with valid observations only
  df_sim <- as.data.frame(cbind(X, age_death, age_start, linpred))[indx_obs,]
  
  
  cat("Range of linpred values:", range(linpred), "\n")
  if(any(is.infinite(exp(linpred)))) {
    warning("Some exp(linpred) values are Inf - consider scaling betas down further")
  }
  
  # Get mean residual life
  for (i in 1:nrow(df_sim)){
    mrl_uncon <- integrate(gomp_baseline_surv, lower = (df_sim$age_start[i] * exp(df_sim$linpred[i])), 
                           upper=Inf, a = a, b = b)$value
    s_cond <-  gomp_baseline_surv(df_sim$age_start[i] * exp(df_sim$linpred[i]), a = a, b = b) 
    df_sim$mrl[i] <- (mrl_uncon / s_cond) * exp(-df_sim$linpred[i])
  }
  
  # Get biological age (via population lifetable)
  for (i in 1:nrow(df_sim)){
    df_sim$b_age[i] <- lt$t[ which.min(abs(lt$mrl - df_sim$mrl[i])) ]
  }
  
  # Add censoring 
  df_sim$yrs_rem <- df_sim$age_death - df_sim$age_start
  wh <- which(df_sim$yrs_rem > followup) # censored
  df_sim$status <- 1
  df_sim$status[wh] <- 0
  df_sim$follow_up_time <- df_sim$yrs_rem
  df_sim$follow_up_time[wh] <- followup
  df_sim$age_end <- df_sim$age_start + df_sim$follow_up_time
  
  return(df_sim)
  
}


sim_grouped_betas_pberends = function(p, g) {
  # Function I created to experiment with grouped random variable creation
  # Overly complex for no reason????
  p_g = p/g
  
  # Generate hyperspace
  anchors = matrix(data = runif(g^2,-4,4), nrow = g, ncol = g)
  anchors = anchors %>% apply(MARGIN = 2, FUN = cummean)
  # magnitudes = rnorm(g)
  
  # anchors = anchors %*% diag(magnitudes)
  group_shift = rnorm(g, mean = 0, sd = 40)  # Increase sd for more spread
  
  # Jittered values will be of mean anchor and sd magnitude
  betas = matrix(NA, nrow = g, ncol = p)
  for (group in 0:(g-1)) {
    for (beta in 1:p_g) {
      mu = anchors[, group + 1]
      sigma = 1 #sqrt(sum(mu^2))
      beta_g = rnorm(length(mu), mean = mu, sd = sigma)
      betas[,group*p_g+beta] = beta_g
    }
  }
  # bring down to one dimension (1 x p from g x p)
  betas = as.tibble(t(matrix(1,1,g) %*% betas))
  betas = betas + rep(group_shift, each = p_g)
  
  betas$group = as.factor(rep(1:g, each = p_g))
  colnames(betas)[1] = "beta"
  betas$beta = scale(betas$beta, center = F)
  return(betas)
}



generate_betas = function(p, g, rho, rho_between, seed,
                          mu_u, mu_l, beta_scale, plot,
                          non_zero_groups, active_hazard) {
  
  # # debug
  # p = 200; g = 20; rho = 0.9; rho_between = 0.2; seed = 123; mu_u = 2; mu_l = -2
  # beta_scale = 0.5; plot = T; non_zero_groups = 0.5
  
  if (!is.null(seed)) set.seed(seed)
  library(MASS)
  stopifnot(p %% g == 0)
  
  p_g = p/g
  group_membership = rep(1:g, each = p_g)
  # Sample Beta means from uniform distribution
  mus = runif(g, mu_l, mu_u)
  betas = numeric(p)
  
  # Randomly select which groups will have non-zero coefficients
  active_groups = sample(1:g, size = round(non_zero_groups * g))
  
  sigma_full = matrix(0, p, p)
  
  # Fill each block based on its distance from diagonal
  for (i in 1:g) {
    for (j in 1:g) {
      # Indices for the current block
      i_idx = ((i-1) * p_g + 1):(i * p_g)
      j_idx = ((j-1) * p_g + 1):(j * p_g)
      
      # Fill in with appropriate correlation based on block distance
      block_distance = abs(i - j)
      if (block_distance == 0) {
        # On diagonal - within group correlation
        sigma_full[i_idx, j_idx] = toeplitz(pmax(rho^(0:(p_g-1)), rep(rho_between, p_g)))
      } else {
        # Correlation decreases with group distance 
        sigma_full[i_idx, j_idx] = rho_between#^block_distance
        
      }
    }
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # Create the mean vector
  mu_full = numeric(p)
  for (group in active_groups) {
    idx = ((group-1) * p_g + 1):(group * p_g)
    mu_full[idx] = mus[group]
  }
  
  # Generate all betas using block matrix
  betas = mvrnorm(1, mu = mu_full, Sigma =  sigma_full * beta_scale)
  
  
  # Set inactive groups to exactly zero for Betas
  for (group in setdiff(1:g, active_groups)) {
    idx = ((group-1) * p_g + 1):(group * p_g)
    betas[idx] = 0
  }
  
  beta_df = tibble(beta = betas, group = group_membership)
  print(glue("Mean of betas pre-scaling: {mean(beta_df$beta)}"))
  print(glue("{c('Lower', 'Upper')} range of beta values pre-scaling: {range(beta_df$beta)}"))
  # Scale active betas to make E(betas) = 0.05
  beta_df = beta_df %>%
    mutate(beta = if_else(group %in% active_groups,
                          (beta - mean(beta_df$beta[beta_df$group %in% active_groups], 
                                      na.rm = TRUE)) + active_hazard,
                          beta))
  
  # Adjust beta magnitude to have a linear predictor of E[pred] = 0, var = 1
  # current_var = sum(beta_df$beta^2)
  # current_sd = sqrt(current_var)
  # 
  # print(glue("Mean of betas post-scaling: {mean(beta_df$beta)}"))
  # 
  # 
  # print(glue("{c('Lower', 'Upper')} range of beta values post-scaling: {range(beta_df$beta)}"))
  
  target_variance = 0.2^2  # Target sd as variance
  # Because we assume X ~ N(0,1)
  # We try to set Var(XB) -> B'Var(X)B -> B'1B -> ||B^2||
  
  # Then we set it so some v to rescale this
  # v ||B^2|| = target
  # v = target / ||B^2||
  # sqrt(v) = sqrt(target / ||B^2||)
  current_variance = sum(beta_df$beta^2) * 0.2
  scale_factor = sqrt(target_variance / current_variance)
  beta_df$beta = beta_df$beta * scale_factor
  
  print(glue("Mean of betas post-scaling: {mean(beta_df$beta)}"))
  print(glue("{c('Lower', 'Upper')} range of beta values post-scaling: {range(beta_df$beta)}"))
  
  if (plot) {
    heatmap(sigma_full, Rowv = NA, Colv = NA, scale = "none", main = "Betas sigma")
    ((betas) %*% t(betas)) %>% heatmap(Rowv = NA, Colv = NA, scale = "none", main = "Beta Gram Matrix")
    
    #  boxplot visualization for betas
    beta_by_group <- split(beta_df$beta, rep(1:g, each = p_g))
    boxplot(beta_by_group, 
            main = "Beta Coefficients Distribution by Group",
            xlab = "Group",
            ylab = "Coefficient Value",
            col = rainbow(g),
            border = "black",
            outline = TRUE)
    # Add a horizontal line at y=0 for reference
    abline(h = 0, lty = 2, col = "gray50")
  }
  return(beta_df)
}

generate_X = function(n, p, g, rho, rho_between, seed = NULL, 
                      scale = 0.5, X_plots = T) {
  # https://projecteuclid.org/journals/annals-of-applied-statistics/volume-7/issue-3/A-method-for-generating-realistic-correlation-matrices/10.1214/13-AOAS638.full
  # https://en.wikipedia.org/wiki/Block_matrix#Block_Toeplitz_matrices
  if (!is.null(seed)) set.seed(seed)
  library(MASS)
  stopifnot(p %% g == 0)
  p_g = p/g
  n_g = n/g
  
  # --------------Create a block Toeplitz matrix for X--------------
  sigma_var = matrix(0, p, p)
  
  # Fill each block based on its distance from diagonal
  for (i in 1:g) {
    for (j in 1:g) {
      # Indices for the current block
      i_idx = ((i-1) * p_g + 1):(i * p_g)
      j_idx = ((j-1) * p_g + 1):(j * p_g)
      
      # Fill in with appropriate correlation based on block distance
      block_distance = abs(i - j) # GROUP VALUE
      if (block_distance == 0) { # On diagonal
        # Within group correlation, pmax takes the highest value so within goup can not be lower than between
        sigma_var[i_idx, j_idx] = toeplitz(pmax(rho^(0:(p_g-1)), rep(rho_between, p_g)))
        # sigma_var[i_idx, j_idx] = toeplitz(rho^(0:(p_g-1)))
      } else { # Off Diagonal
        # Correlation decreases with group distance 
        sigma_var[i_idx, j_idx] = rho_between
      }
    }
  }
  
  
  # Generate all covariates using block matrix
  # X = matrix(rnorm(n*p, mean = 0, sd = scale), nrow = n, ncol = p)
  X = mvrnorm(n, mu = rep(0, p), Sigma =  sigma_var * scale)

  # Vector indicating group memberships
  group_membership = rep(1:g, each = p_g)
  
  
  X_df = as.tibble(X)
  
  
  # ----------Plot metrics--------------
  if (X_plots) {
    # Visualize correlation structure
    heatmap(sigma_var, Rowv = NA, Colv = NA, scale = "none", main = "Sigma structure")
    
    # # Visualize X matrix subsample
    # X_sample = X[1:min(100, n), 1:min(100, p)]
    heatmap(X, Rowv = NA, Colv = NA, scale = "none", main = "X Matrix Sample")
    
    X_long = data.frame(
      Value = as.vector(X),
      Group = rep(rep(1:g, each = p_g), n),
      Variable = rep(1:p, each = n)
    )
    
    # boxplot of all variables by group
    boxplot(Value ~ Group, data = X_long,
            main = "Distribution of X Values by Group",
            xlab = "Group", 
            ylab = "Value",
            col = rainbow(g),
            border = "black")
    abline(h = 0, lty = 2, col = "gray50")
    
    # Correlation of X
    cor(X_df) %>% heatmap(Rowv = NA, Colv = NA, scale = "none", main = "X Matrix Correlations")
  }
  
  colnames(X_df) = paste0("X", 1:p)
  return(list(
    X = X_df,
    group_membership = group_membership
  ))
}
