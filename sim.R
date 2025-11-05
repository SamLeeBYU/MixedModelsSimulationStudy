# Author: Connor Thompson, Daisy Harris, Sam Lee

library(lme4)
library(nlme)
library(lmerTest)

# Requires: lme4, nlme
# Assumes you already defined get_data(beta_0, beta_1, covariance_block, n_sub, n_obs_per)

simulate <- function(
  V_list, # list(CS=..., AR1=..., RC=...) covariance blocks (n_obs_per x n_obs_per)
  n_replications = 100,
  beta_0 = 1,
  beta_1 = 0, # set 0 for size; >0 for power
  n_subjects = 10,
  n_obs_per = 5,
  alpha = 0.05,
  show_progress = TRUE
) {
  stopifnot(
    is.list(V_list),
    all(c("CS", "AR1", "RC") %in% names(V_list)),
    nrow(V_list$CS) == n_obs_per,
    nrow(V_list$AR1) == n_obs_per,
    nrow(V_list$RC) == n_obs_per
  )

  # results matrix
  results <- matrix(NA_real_, nrow = n_replications, ncol = 6)
  colnames(results) <- c(
    "compound_lrt_p",
    "compound_f_p",
    "ar_lrt_p",
    "ar_f_p",
    "rancoef_lrt_p",
    "rancoef_f_p"
  )

  # progress bar
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = n_replications, style = 3)
  }

  for (rep in 1:n_replications) {
    ## ----- 1) Compound symmetry block: fit RI LMM -----
    dat <- get_data(
      beta_0,
      beta_1,
      V_list$CS,
      n_sub = n_subjects,
      n_obs_per = n_obs_per
    )

    # LRT (ML)
    fit_full_ml <- lme4::lmer(
      response ~ -1 + treatment * time + (1 | subject),
      data = dat,
      REML = FALSE
    )
    fit_null_ml <- lme4::lmer(
      response ~ -1 + treatment + time + (1 | subject),
      data = dat,
      REML = FALSE
    )
    lrt_tab <- stats::anova(fit_full_ml, fit_null_ml)
    results[rep, "compound_lrt_p"] <- lrt_tab$`Pr(>Chisq)`[2]

    # Wald (REML) for treatment:time
    fit_reml <- lme4::lmer(
      response ~ -1 + treatment * time + (1 | subject),
      data = dat,
      REML = TRUE
    )
    C <- matrix(c(0, 0, 0, 1), nrow = 1) # contrast for treatment:timetreatment1:time
    bhat <- lme4::fixef(fit_reml)
    Vb <- as.matrix(stats::vcov(fit_reml))
    W <- as.numeric(t(C %*% bhat) %*% solve(C %*% Vb %*% t(C)) %*% (C %*% bhat)) # chi-square(1)
    results[rep, "compound_f_p"] <- stats::pf(
      W,
      df1 = 1,
      df2 = 2 * n_subjects,
      lower.tail = FALSE
    )

    ## ----- 2) AR(1) block: fit GLS with AR(1) -----
    dat <- get_data(
      beta_0,
      beta_1,
      V_list$AR1,
      n_sub = n_subjects,
      n_obs_per = n_obs_per
    )

    # LRT (ML)
    fit_full_gls_ml <- nlme::gls(
      response ~ -1 + treatment * time,
      correlation = nlme::corAR1(form = ~ 1 | subject),
      data = dat,
      method = "ML"
    )
    fit_null_gls_ml <- nlme::gls(
      response ~ -1 + treatment + time,
      correlation = nlme::corAR1(form = ~ 1 | subject),
      data = dat,
      method = "ML"
    )
    lrt_tab <- stats::anova(fit_full_gls_ml, fit_null_gls_ml)
    results[rep, "ar_lrt_p"] <- lrt_tab$`p-value`[2]

    # Wald (REML)
    fit_gls_reml <- nlme::gls(
      response ~ -1 + treatment * time,
      correlation = nlme::corAR1(form = ~ 1 | subject),
      data = dat,
      method = "REML"
    )
    bhat <- stats::coef(fit_gls_reml)
    Vb <- stats::vcov(fit_gls_reml)
    # Coef order: (-1 + treatment*time) -> c(treatment0, treatment1, time, treatment1:time)
    # Same contrast:
    W <- as.numeric(t(C %*% bhat) %*% solve(C %*% Vb %*% t(C)) %*% (C %*% bhat))
    results[rep, "ar_f_p"] <- stats::pf(
      W,
      df1 = 1,
      df2 = 2 * n_subjects,
      lower.tail = FALSE
    )

    ## ----- 3) Random coefficients block: fit RI+RS LMM -----
    dat <- get_data(
      beta_0,
      beta_1,
      V_list$RC,
      n_sub = n_subjects,
      n_obs_per = n_obs_per
    )

    # LRT (ML)
    fit_full_ml <- lme4::lmer(
      response ~ -1 + treatment * time + (time | subject),
      data = dat,
      REML = FALSE
    )
    fit_null_ml <- lme4::lmer(
      response ~ -1 + treatment + time + (time | subject),
      data = dat,
      REML = FALSE
    )
    lrt_tab <- stats::anova(fit_full_ml, fit_null_ml)
    results[rep, "rancoef_lrt_p"] <- lrt_tab$`Pr(>Chisq)`[2]

    # Wald (REML)
    fit_reml <- lme4::lmer(
      response ~ -1 + treatment * time + (time | subject),
      data = dat,
      REML = TRUE
    )
    bhat <- lme4::fixef(fit_reml)
    Vb <- as.matrix(stats::vcov(fit_reml))
    W <- as.numeric(t(C %*% bhat) %*% solve(C %*% Vb %*% t(C)) %*% (C %*% bhat))
    results[rep, "rancoef_f_p"] <- stats::pf(
      W,
      df1 = 1,
      df2 = 2 * n_subjects,
      lower.tail = FALSE
    )

    if (show_progress) utils::setTxtProgressBar(pb, rep)
  }
  if (show_progress) {
    close(pb)
  }

  # Return list: raw p-values and size/power summaries at alpha=0.05
  list(
    results = results,
    reject_rates_alpha = colMeans(results < alpha, na.rm = TRUE)
  )
}
