# Simple test of the inference function
library(yaml)
library(jsonlite)

# Source the functions
source("R/utilities.R")
source("R/mcmc.R")
source("R/analysis.R")
source("R/po_inference.R")

# Load the generated data
cat("Loading test data...\n")
data <- jsonlite::fromJSON("results/simple_test_data.json", simplifyVector = FALSE)

# Create a simple MCMC config
config <- list(
  generation = list(K = 2),
  mcmc = list(
    num_iterations = 500,
    random_seed = 42,
    drbeta = 0.1,
    sigma_beta = 1.0,
    update_probabilities = list(
      rho = 0.25,
      noise = 0.25,
      U = 0.4,
      beta = 0.1
    )
  ),
  rho = list(dr = 1.1),
  noise = list(
    noise_option = "queue_jump",
    sigma_mallow = 0.1
  ),
  prior = list(
    rho_prior = 5.0,
    noise_beta_prior = 2.0,
    mallow_ua = 1.0
  ),
  covariates = list(p = 2),
  visualization = list(burn_in = 100)
)

# Run inference
cat("Running MCMC inference...\n")
results <- run_inference(data, config, use_rj_mcmc = FALSE)

cat("Inference completed!\n")
cat("Inferred edges:", sum(results$h), "\n")

# Convert true partial order properly
true_po <- do.call(rbind, lapply(data$true_partial_order, function(x) unlist(x)))
cat("True edges:", sum(true_po), "\n")
cat("Acceptance rate:", round(results$overall_acceptance_rate * 100, 1), "%\n")

# Save results
save_results(results, "results", "simple_inference_test")
cat("Results saved!\n")

# Show comparison
cat("\nTrue partial order:\n")
print(true_po)
cat("\nInferred partial order:\n")
print(results$h) 