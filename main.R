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
