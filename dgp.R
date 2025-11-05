# Author: Sam Lee
# Disclaimer: ChatGPT-5o helped integrate the cholesky decomposition into
# existing code written by Connor Thompson

get_data <- function(beta_0, beta_1, covariance_block, n_sub, n_obs_per = 5) {
  stopifnot(is.matrix(covariance_block), nrow(covariance_block) == n_obs_per)

  # Cholesky: chol(V) returns R with V = R'R. Use lower L = t(R).
  L_block <- t(chol(covariance_block))
  L_full <- diag(n_sub) %x% L_block # block-diagonal Cholesky
  p <- n_sub * n_obs_per # observations per group

  # Control group
  times_c <- rep(0:(n_obs_per - 1), times = n_sub)
  mu_c <- rep(beta_0, p)
  z_c <- rnorm(p)
  y_c <- as.numeric(mu_c + L_full %*% z_c)
  control_df <- data.frame(
    time = times_c,
    response = y_c,
    treatment = factor(0L, levels = c(0, 1))
  )

  # Treated group
  times_t <- rep(0:(n_obs_per - 1), times = n_sub)
  mu_t <- rep(beta_0, p) + beta_1 * times_t
  z_t <- rnorm(p)
  y_t <- as.numeric(mu_t + L_full %*% z_t)
  treated_df <- data.frame(
    time = times_t,
    response = y_t,
    treatment = factor(1L, levels = c(0, 1))
  )

  dat <- rbind(control_df, treated_df)
  dat$subject <- factor(rep(1:(2 * n_sub), each = n_obs_per))
  rownames(dat) <- NULL
  dat
}

datCS <- get_data(
  1,
  1,
  V_list$CS,
  n_sub = 1,
  n_obs_per = 5
)
datCS$CovStructure <- "Compounds Symmetric"
datRC <- get_data(
  1,
  1,
  V_list$RC,
  n_sub = 1,
  n_obs_per = 5
)
datRC$CovStructure <- "Random Coefficients"
datAR <- get_data(
  1,
  1,
  V_list$AR1,
  n_sub = 1,
  n_obs_per = 5
)
datAR$CovStructure <- "AR(1)"

#Sample for paper
# rbind(datCS, rbind(datRC, datAR))
