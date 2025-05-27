# Example data for BayesianPO package
# This script creates example datasets that users can load

# Load the package
library(BayesianPO)

# Small example dataset for quick testing
small_example <- generate_synthetic_data(
  n_items = 4,
  n_observations = 20,
  K = 2,
  rho_true = 0.7,
  prob_noise_true = 0.1,
  random_seed = 42
)

# Medium example dataset for demonstration
medium_example <- generate_synthetic_data(
  n_items = 6,
  n_observations = 50,
  K = 3,
  rho_true = 0.8,
  prob_noise_true = 0.15,
  random_seed = 123
)

# Save the datasets
save(small_example, file = "inst/extdata/small_example.rda")
save(medium_example, file = "inst/extdata/medium_example.rda") 