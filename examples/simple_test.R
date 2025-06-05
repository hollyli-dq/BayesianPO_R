# Simple test of the comprehensive functions
library(yaml)
library(jsonlite)

# Source the functions
source("R/utilities.R")
source("R/mcmc.R") 
source("R/data_generator.R")

# Test data generation
cat("Testing data generation...\n")

# Create a simple config
config <- list(
  generation = list(n = 4, N = 20, prob_noise_true = 0.1),
  mcmc = list(K = 2),
  prior = list(rho_prior = 5.0, mallow_ua = 1.0, noise_beta_prior = 2.0),
  noise = list(noise_option = "queue_jump"),
  covariates = list(p = 2, beta_true = c(0.5, -0.3)),
  min_sub = 3.0
)

# Generate data
data <- generate_data(config)

cat("Data generation successful!\n")
cat("Generated", length(data$observed_orders), "observed orders\n")
cat("True partial order has", sum(data$true_partial_order), "edges\n")

# Save the data
if (!dir.exists("results")) dir.create("results")
jsonlite::write_json(data, "results/simple_test_data.json", pretty = TRUE, auto_unbox = TRUE)
cat("Data saved to results/simple_test_data.json\n") 