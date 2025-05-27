# Simple CLI test bypassing YAML issues
library(jsonlite)

# Source the functions
source("R/utilities.R")
source("R/mcmc.R")
source("R/mcmc_rj.R")
source("R/analysis.R")
source("R/data_generator.R")
source("R/po_inference.R")

cat("=== Testing Comprehensive Functions ===\n\n")

# Test 1: Data Generation
cat("1. Data Generation\n")
cat("------------------\n")

# Create config directly
data_config <- list(
  generation = list(n = 4, N = 20, prob_noise_true = 0.1),
  mcmc = list(K = 2),
  prior = list(rho_prior = 5.0, mallow_ua = 1.0, noise_beta_prior = 2.0),
  noise = list(noise_option = "queue_jump"),
  covariates = list(p = 2, beta_true = c(0.5, -0.3))
)

# Generate data
data <- generate_data(data_config)

# Save data manually
if (!dir.exists("results")) dir.create("results")
jsonlite::write_json(data, "results/cli_test_data.json", pretty = TRUE, auto_unbox = TRUE)

cat("Data generation completed!\n\n")

# Test 2: Fixed Dimension MCMC
cat("2. Fixed Dimension MCMC\n")
cat("-----------------------\n")

mcmc_config <- list(
  generation = list(K = 2),
  mcmc = list(
    num_iterations = 500,
    random_seed = 42,
    drbeta = 0.1,
    sigma_beta = 1.0,
    update_probabilities = list(rho = 0.25, noise = 0.25, U = 0.4, beta = 0.1)
  ),
  rho = list(dr = 1.1),
  noise = list(noise_option = "queue_jump", sigma_mallow = 0.1),
  prior = list(rho_prior = 5.0, noise_beta_prior = 2.0, mallow_ua = 1.0),
  covariates = list(p = 2),
  visualization = list(burn_in = 100)
)

results_fixed <- run_inference(data, mcmc_config, use_rj_mcmc = FALSE)
save_results(results_fixed, "results", "cli_test_fixed")

cat("Fixed MCMC completed!\n")
cat("Inferred edges:", sum(results_fixed$h), "\n")
cat("Acceptance rate:", round(results_fixed$overall_acceptance_rate * 100, 1), "%\n\n")

# Test 3: Reversible Jump MCMC
cat("3. Reversible Jump MCMC\n")
cat("-----------------------\n")

mcmc_config$mcmc$update_probabilities$K <- 0.1
mcmc_config$mcmc$update_probabilities$U <- 0.3
mcmc_config$prior$K_prior <- 2.0

results_rj <- run_inference(data, mcmc_config, use_rj_mcmc = TRUE)
save_results(results_rj, "results", "cli_test_rj")

cat("RJ-MCMC completed!\n")
cat("Inferred edges:", sum(results_rj$h), "\n")
cat("Acceptance rate:", round(results_rj$overall_acceptance_rate * 100, 1), "%\n")

if (!is.null(results_rj$K_trace)) {
  final_K <- results_rj$K_trace[length(results_rj$K_trace)]
  cat("Final K:", final_K, "\n")
}

cat("\n=== All Tests Completed Successfully ===\n") 