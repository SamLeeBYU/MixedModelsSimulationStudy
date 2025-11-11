source("make_cov.R")

library(ggplot2)
V_list <- make_cov_structures(n = 5, D = 6)

source("dgp.R")

source("sim.R")

set.seed(537)

data <- get_data(1,2,V_list$RC,3)

make_data_plot <- function(data){
  ggplot(data, aes(x = time, y = response, group = subject, color = treatment)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    scale_color_manual(values = c("0" = "#0072B2", "1" = "#E69F00"),
                       labels = c("Control", "Treated"),
                       name = "Treatment") +
    labs(
      title = "Simulated Repeated Measures by Subject",
      x = "Time",
      y = "Response"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

make_power_curve <- function(beta_1_values, n_replications = 100, n_subjects = 10){
  run_sim <- function(beta_1){
    out <- simulate(
      V_list = V_list,
      n_replications = n_replications,
      beta_0 = 1,
      beta_1 = beta_1,
      n_subjects = n_subjects,
      n_obs_per = 5
    )
    power <- out$reject_rates_alpha
    # Return a long-format data frame
    df <- data.frame(
      beta_1 = beta_1,
      Test = names(power),
      Power = as.numeric(power)
    )
  }
  
  power_results <- do.call(rbind, lapply(beta_1_values, run_sim))
  
  test_labels <- c(
    ar_f_p = "AR(1) F-test",
    compound_f_p = "Compound Symmetry F-test",
    rancoef_f_p = "Random Coefficients F-test",
    ar_lrt_p = "AR(1) LRT",
    compound_lrt_p = "Compound Symmetry LRT",
    rancoef_lrt_p = "Random Coefficients LRT"
  )
  
  ggplot(power_results, aes(x = beta_1, y = Power, 
                            color = Test, linetype = Test, group = Test)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c(
      "ar_f_p" = "#0072B2",         
      "compound_f_p" = "#E69F00",   
      "rancoef_f_p" = "#009E73",    
      "ar_lrt_p" = "#0072B2",       
      "compound_lrt_p" = "#E69F00", 
      "rancoef_lrt_p" = "#009E73"   
    ), labels = test_labels, name = "Test") +
    scale_linetype_manual(values = c(
      "ar_f_p" = "solid",
      "compound_f_p" = "solid",
      "rancoef_f_p" = "solid",
      "ar_lrt_p" = "dashed",
      "compound_lrt_p" = "dashed",
      "rancoef_lrt_p" = "dashed"
    ), labels = test_labels, name = "Test") +
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


set.seed(537)
run_sim <- function(beta_1, n_subjects, alpha = 0.05) {
  out <- simulate(
    V_list = V_list,
    n_replications = 1000,
    beta_0 = 1,
    beta_1 = beta_1,
    n_subjects = n_subjects,
    n_obs_per = 5
  )
  
  power <- out$reject_rates_alpha
  
  pvals <- out$results
  
  rejects <- (pvals < alpha)
  
  n_rep <- nrow(pvals)
  mcse <- apply(rejects, 2, function(x) {
    p <- mean(x)
    sqrt(p * (1 - p) / n_rep)
  })
  
  df <- data.frame(
    beta_1 = beta_1,
    Test = names(power),
    Power = as.numeric(power),
    MCSE = as.numeric(mcse[names(power)])  # align names just in case
  )
  
  return(df)
}


set.seed(537)
power_results_5 <- do.call(
  rbind,
  lapply(seq(0, 1.625, by = 0.125), function(b) run_sim(beta_1 = b, n_subjects = 5))
)

power_results_10 <- do.call(
  rbind,
  lapply(seq(0, 1.625, by = 0.125), function(b) run_sim(beta_1 = b, n_subjects = 10))
)

power_results_25 <- do.call(
  rbind,
  lapply(seq(0, 1.625, by = 0.125), function(b) run_sim(beta_1 = b, n_subjects = 25))
)

power_results_50 <- do.call(rbind,
  lapply(seq(0, 1.625, by = 0.125), function(b) run_sim(beta_1 = b, n_subjects = 50))
)
power_results_5
power_results_10
power_results_25
power_results_50

power_results <- list(power_results_5 = power_results_5,
                      power_results_10 = power_results_10,
                      power_results_25 = power_results_25,
                      power_results_50 = power_results_50)

saveRDS(power_results, file = "power_results.rds")


test_labels <- c(
  ar_f_p = "AR(1) F-test",
  compound_f_p = "Compound Symmetry F-test",
  rancoef_f_p = "Random Coefficients F-test",
  ar_lrt_p = "AR(1) LRT",
  compound_lrt_p = "Compound Symmetry LRT",
  rancoef_lrt_p = "Random Coefficients LRT"
)

library(patchwork)

power_curve_5 <- ggplot(power_results_5, aes(x = beta_1, y = Power, 
                                             color = Test, linetype = Test, group = Test)) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "ar_f_p" = "#0072B2",         
    "compound_f_p" = "#E69F00",   
    "rancoef_f_p" = "#009E73",    
    "ar_lrt_p" = "#0072B2",       
    "compound_lrt_p" = "#E69F00", 
    "rancoef_lrt_p" = "#009E73"   
  ), labels = test_labels, name = "Test") +
  scale_linetype_manual(values = c(
    "ar_f_p" = "solid",
    "compound_f_p" = "solid",
    "rancoef_f_p" = "solid",
    "ar_lrt_p" = "dashed",
    "compound_lrt_p" = "dashed",
    "rancoef_lrt_p" = "dashed"
  ), labels = test_labels, name = "Test") +
  labs(
    title = "n = 5 Subjects",
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

power_curve_10 <- ggplot(power_results_10, aes(x = beta_1, y = Power, 
                                               color = Test, linetype = Test, group = Test)) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "ar_f_p" = "#0072B2",         
    "compound_f_p" = "#E69F00",   
    "rancoef_f_p" = "#009E73",    
    "ar_lrt_p" = "#0072B2",       
    "compound_lrt_p" = "#E69F00", 
    "rancoef_lrt_p" = "#009E73"   
  ), labels = test_labels, name = "Test") +
  scale_linetype_manual(values = c(
    "ar_f_p" = "solid",
    "compound_f_p" = "solid",
    "rancoef_f_p" = "solid",
    "ar_lrt_p" = "dashed",
    "compound_lrt_p" = "dashed",
    "rancoef_lrt_p" = "dashed"
  ), labels = test_labels, name = "Test") +
  labs(
    title = "n = 10 Subjects",
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

power_curve_25 <- ggplot(power_results_25, aes(x = beta_1, y = Power, 
                                               color = Test, linetype = Test, group = Test)) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "ar_f_p" = "#0072B2",         
    "compound_f_p" = "#E69F00",   
    "rancoef_f_p" = "#009E73",    
    "ar_lrt_p" = "#0072B2",       
    "compound_lrt_p" = "#E69F00", 
    "rancoef_lrt_p" = "#009E73"   
  ), labels = test_labels, name = "Test") +
  scale_linetype_manual(values = c(
    "ar_f_p" = "solid",
    "compound_f_p" = "solid",
    "rancoef_f_p" = "solid",
    "ar_lrt_p" = "dashed",
    "compound_lrt_p" = "dashed",
    "rancoef_lrt_p" = "dashed"
  ), labels = test_labels, name = "Test") +
  labs(
    title = "n = 25 Subjects",
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

power_curve_50 <- ggplot(power_results_50, aes(x = beta_1, y = Power, 
                                               color = Test, linetype = Test, group = Test)) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "ar_f_p" = "#0072B2",         
    "compound_f_p" = "#E69F00",   
    "rancoef_f_p" = "#009E73",    
    "ar_lrt_p" = "#0072B2",       
    "compound_lrt_p" = "#E69F00", 
    "rancoef_lrt_p" = "#009E73"   
  ), labels = test_labels, name = "Test") +
  scale_linetype_manual(values = c(
    "ar_f_p" = "solid",
    "compound_f_p" = "solid",
    "rancoef_f_p" = "solid",
    "ar_lrt_p" = "dashed",
    "compound_lrt_p" = "dashed",
    "rancoef_lrt_p" = "dashed"
  ), labels = test_labels, name = "Test") +
  labs(
    title = "n = 50 Subjects",
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

(power_curve_5 + power_curve_10 + power_curve_25 + power_curve_50) + plot_layout(ncol = 4, guides = "collect") & theme(legend.position = "bottom")


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # for side-by-side layout

# helper function to convert to data frame for ggplot
mat_to_df <- function(mat, name) {
  as.data.frame(mat) |>
    mutate(Row = factor(1:nrow(mat))) |>
    pivot_longer(-Row, names_to = "Col", values_to = "Value") |>
    mutate(Col = factor(gsub("V", "", Col)), Matrix = name)
}

df_all <- bind_rows(
  mat_to_df(V_list$CS, "CS"),
  mat_to_df(V_list$AR1, "AR1"),
  mat_to_df(V_list$RC, "RC")
)

df_all <- df_all %>%
  mutate(Label = sprintf("%.2f", Value))

covariance_heatmaps <- ggplot(df_all, aes(x = Col, y = Row, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), color = "white", size = 3, fontface = "bold") +
  scale_fill_viridis_c(option = "turbo") +
  scale_y_discrete(limits = rev(levels(df_all$Row))) + # reverse y-axis
  coord_equal() +
  facet_wrap(~Matrix, nrow = 1) +
  theme_minimal(base_size = 12) +
  labs(x = "Column", y = "Row", fill = "Value") +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

# display
print(covariance_heatmaps)

# save to PDF
ggsave(
  "covariance_heatmaps.pdf",
  plot = covariance_heatmaps,
  device = "pdf",
  width = 10,
  height = 4
)
