# Author: Sam Lee
# Disclaimer: ChatGPT 5-o helped me code up the rank-1 and rank-2 determinant update formulas

#This code creates one block for each specified kind of covariance structure:
# - Compound Symmetric
# - AR(1)
# - Random Coefficients (with a parameters for random intercept and slope)
#These covariance structures are created simultaneously as to fix the marginal variance of each block
# by fixing the determinant (D) and number of observations in each group (n)

#We also fix:
# sigma2_e_cs = 1, # for Compound Symmetry
# phi_ar = 0.5, # for AR(1)

# sigma2_e_rc = 0.6, # for Random Coefficients
#sigma2_int_rc = 0.45

solve_compound_sigma2_block <- function(D, n, sigma2_e) {
  (D / (sigma2_e^(n - 1)) - sigma2_e) / n
}

solve_ar_sigma2 <- function(D, n, phi) {
  (D / ((1 - phi^2)^(n - 1)))^(1 / n)
}

solve_rancoef_sigma2_slp <- function(D, n, sigma2_e, sigma2_int) {
  S11 <- n
  S12 <- n * (n + 1) / 2
  S22 <- n * (n + 1) * (2 * n + 1) / 6
  a <- sigma2_int / sigma2_e
  K <- D / (sigma2_e^n)
  discr <- S22 + a * (S11 * S22 - S12^2)
  if (discr <= 0) {
    stop("Infeasible parameters: denominator <= 0")
  }
  b <- (K - 1 - a * S11) / discr
  if (b < 0) {
    stop("Infeasible: implied slope variance < 0")
  }
  b * sigma2_e
}

# Convenience: choose D from one block, then match others
common_D_from_compound <- function(n, sigma2_e, sigma2_b) {
  sigma2_e^(n - 1) * (sigma2_e + n * sigma2_b)
}

make_cov_structures <- function(
  n,
  D,
  sigma2_e_cs = 1, # for Compound Symmetry
  phi_ar = 0.5, # for AR(1)
  sigma2_e_rc = 0.6, # for Random Coefficients
  sigma2_int_rc = 0.45
) {
  # for Random Coefficients
  # 1) Compound Symmetric
  sigma2_b_cs <- solve_compound_sigma2_block(D, n, sigma2_e_cs)
  V_cs <- sigma2_b_cs * matrix(1, n, n) + sigma2_e_cs * diag(n)

  # 2) AR(1)
  sigma2_ar <- solve_ar_sigma2(D, n, phi_ar)
  V_ar <- sigma2_ar * phi_ar^abs(outer(1:n, 1:n, "-"))

  # 3) Random Coefficients (intercept + slope)
  # Solve for slope variance
  sigma2_slp_rc <- solve_rancoef_sigma2_slp(D, n, sigma2_e_rc, sigma2_int_rc)

  # Construct RC covariance:
  # V = sigma_{int}^2 * 1 1' + sigma_{slp}^2 * t t' + sigma_e^2 I
  t <- 1:n
  V_rc <- sigma2_int_rc *
    outer(rep(1, n), rep(1, n)) +
    sigma2_slp_rc * outer(t, t) +
    sigma2_e_rc * diag(n)

  # Return the structures
  return(list(
    CS = V_cs,
    AR1 = V_ar,
    RC = V_rc
  ))
}
