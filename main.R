#This creates covariance blocks for a determined marginal variance (D)
# and observations per sample (n)
source("make_cov.R")

V_list <- make_cov_structures(n = 5, D = 6)

#This creates the data using the above covariance block structures using
#the Cholesky Decomposition
#This defines the get_data function used in sim.R
source("dgp.R")

source("sim.R")

set.seed(537)
out <- simulate(
  V_list = V_list,
  n_replications = 100,
  beta_0 = 1,
  beta_1 = 0,
  n_subjects = 10,
  n_obs_per = 5
)


# This creates the figures for the paper
source("figures.R")

data_rc <- get_data(1,2,V_list$RC,3)
data_plot_rc <- make_data_plot(data)

data_ar1 <- get_data(1,2,V_list$AR1,3)
data_plot_ar1 <- make_data_plot(data_ar1)

data_cs <- get_data(1,2,V_list$CS,3)
data_plot_cs <- make_data_plot(data_cs)

power_curve <- make_power_curve(
  beta_1_values = seq(0, 1.6, by = 0.1),
  n_replications = 100,
  n_subjects = 10
)
