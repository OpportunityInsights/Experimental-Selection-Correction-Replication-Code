############################################################
##  exp_sc_boot : Experimental-Selection Correction
############################################################
# exp_df   – experimental data (randomized)
# obs_df   – observational data (potentially confounded)
# outcome  – name of long-run outcome Y (string)
# z_vars   – vector of surrogate names Z_g (character vector)
# treat    – name of treatment indicator D (string; default "D")
# aux      – optional covariate names X_v (character vector; default none)
# B        – number of bootstrap replications (default 1000)
#
# returns a list with
#   $tau_hat  – ESC point estimate of causal effect
#   $se_hat   – bootstrap standard error
#   $betas    – first-stage β_g estimates
#   $tau_boot – vector of bootstrap τ estimates
############################################################
exp_sc <- function(exp_df, obs_df,
                        outcome, z_vars,
                        treat = "D", aux = character(0),
                        B = 1000) {
  # 1. Point estimate --------------------------------------------------
  # 1a) first-stage in experimental data
  betas <- setNames(numeric(length(z_vars)), z_vars)
  for (g in seq_along(z_vars)) {
    f1 <- reformulate(c(treat, aux), response = z_vars[g])
    betas[g] <- coef(lm(f1, data = exp_df))[treat]
  }
  
  # 1b) build residuals in obs_df
  obs_aug <- obs_df
  for (g in seq_along(z_vars)) {
    sname <- paste0("s_", z_vars[g])
    obs_aug[[sname]] <- obs_aug[[z_vars[g]]] -
      betas[g] * obs_aug[[treat]]
    if (length(aux)) {
      for (v in aux)
        obs_aug[[sname]] <- obs_aug[[sname]] -
          0  # placeholder for β_{v,g}*X_v if needed
    }
  }
  s_vars <- paste0("s_", z_vars)
  
  # 1c) second-stage regression
  f2 <- reformulate(c(treat, s_vars, aux), response = outcome)
  fit2 <- lm(f2, data = obs_aug)
  tau_hat <- coef(fit2)[treat]
  
  # 2. Bootstrap -------------------------------------------------------
  n_exp <- nrow(exp_df)
  n_obs <- nrow(obs_df)
  tau_boot <- numeric(B)
  for (b in seq_len(B)) {
    # 2a) resample both samples
    exp_b <- exp_df[sample(n_exp, replace = TRUE), , drop = FALSE]
    obs_b <- obs_df[sample(n_obs, replace = TRUE), , drop = FALSE]
    
    # 2b) first-stage on exp_b
    betas_b <- numeric(length(z_vars))
    for (g in seq_along(z_vars)) {
      f1b <- reformulate(c(treat, aux), response = z_vars[g])
      betas_b[g] <- coef(lm(f1b, data = exp_b))[treat]
    }
    
    # 2c) build residuals in obs_b
    obs_b_aug <- obs_b
    for (g in seq_along(z_vars)) {
      sname <- paste0("s_", z_vars[g])
      obs_b_aug[[sname]] <- obs_b_aug[[z_vars[g]]] -
        betas_b[g] * obs_b_aug[[treat]]
      if (length(aux)) {
        for (v in aux)
          obs_b_aug[[sname]] <- obs_b_aug[[sname]] -
            0
      }
    }
    
    # 2d) second-stage on obs_b_aug
    f2b <- reformulate(c(treat, s_vars, aux), response = outcome)
    fit2b <- lm(f2b, data = obs_b_aug)
    tau_boot[b] <- coef(fit2b)[treat]
  }
  
  se_hat <- sd(tau_boot, na.rm = TRUE)
  
  list(
    tau_hat  = tau_hat,
    se_hat   = se_hat,
    betas    = betas,
    tau_boot = tau_boot
  )
}



#  Example: Apply ESC with bootstrap to your dataset
# res <- exp_sc(
#   exp_df, obs_df,
#   outcome = "Y",
#   z_vars  = "Z",
#   treat   = "D",
#   B       = 1000  
# )