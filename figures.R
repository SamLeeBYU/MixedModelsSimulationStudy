source("make_cov.R")

library(ggplot2)
library(patchwork)
V_list <- make_cov_structures(n = 5, D = 6)

source("dgp.R")

source("sim.R")

set.seed(537)

data <- get_data(1, 2, V_list$RC, 3)

make_data_plot <- function(data) {
  ggplot(
    data,
    aes(x = time, y = response, group = subject, color = treatment)
  ) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    scale_color_manual(
      values = c("0" = "#0072B2", "1" = "#E69F00"),
      labels = c("Control", "Treated"),
      name = "Treatment"
    ) +
    labs(
      x = "Time (x)",
      y = "Response (y)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

make_power_curve <- function(beta_1_values, n_replications, n_subjects) {
  run_sim <- function(beta_1) {
    out <- simulate(
      V_list = V_list,
      n_replications = n_replications,
      beta_0 = 1,
      beta_1 = beta_1,
      n_subjects = n_subjects,
      n_obs_per = 5
    )
    power <- out$reject_rates_alpha_0_05
    # Return a long-format data frame
    df <- data.frame(
      beta_1 = beta_1,
      Test = names(power),
      Power = as.numeric(power)
    )
  }

  power_results <- do.call(rbind, lapply(beta_1_values, run_sim))

  test_labels <- c(
    ar_wald_p = "AR(1) F-test",
    compound_wald_p = "Compound Symmetry F-test",
    rancoef_wald_p = "Random Coefficients F-test",
    ar_lrt_p = "AR(1) LRT",
    compound_lrt_p = "Compound Symmetry LRT",
    rancoef_lrt_p = "Random Coefficients LRT"
  )

  ggplot(
    power_results,
    aes(x = beta_1, y = Power, color = Test, linetype = Test, group = Test)
  ) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(
      values = c(
        "ar_wald_p" = "#0072B2",
        "compound_wald_p" = "#E69F00",
        "rancoef_wald_p" = "#009E73",
        "ar_lrt_p" = "#0072B2",
        "compound_lrt_p" = "#E69F00",
        "rancoef_lrt_p" = "#009E73"
      ),
      labels = test_labels,
      name = "Test"
    ) +
    scale_linetype_manual(
      values = c(
        "ar_wald_p" = "solid",
        "compound_wald_p" = "solid",
        "rancoef_wald_p" = "solid",
        "ar_lrt_p" = "dashed",
        "compound_lrt_p" = "dashed",
        "rancoef_lrt_p" = "dashed"
      ),
      labels = test_labels,
      name = "Test"
    ) +
    labs(
      title = "Power Curves Across Covariance Structures",
      x = expression(beta[1]),
      y = "Empirical Power"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    xlim(0, 1.6)
}
