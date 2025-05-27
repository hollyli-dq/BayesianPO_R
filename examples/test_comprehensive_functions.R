# Test script for comprehensive Bayesian Partial Order Inference functions
# This script demonstrates the usage of the large, self-contained functions

# Load required libraries
library(yaml)
library(jsonlite)
library(mvtnorm)

# Source all required functions
source("R/utilities.R")
source("R/mcmc.R")
source("R/mcmc_rj.R")
source("R/analysis.R")
source("R/data_generator.R")
source("R/po_inference.R")

cat("=== Testing Comprehensive Bayesian Partial Order Inference Functions ===\n\n")

# Test 1: Data Generation
cat("1. Testing Data Generation Function\n")
cat("-----------------------------------\n")

# Load data generation config
data_config <- load_config("config/data_generator_config.yaml")

# Generate synthetic data
cat("Generating synthetic data...\n")
synthetic_data <- generate_data(data_config)

# Save the data
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data_path <- file.path(output_dir, "test_synthetic_data.json")
jsonlite::write_json(synthetic_data, data_path, pretty = TRUE, auto_unbox = TRUE)
cat("Data saved to:", data_path, "\n\n")

# Test 2: Fixed Dimension MCMC Inference
cat("2. Testing Fixed Dimension MCMC Inference\n")
cat("------------------------------------------\n")

# Load MCMC config
mcmc_config <- load_config("config/mcmc_config.yaml")

# Run fixed dimension inference
cat("Running fixed dimension MCMC inference...\n")
results_fixed <- run_inference(synthetic_data, mcmc_config, use_rj_mcmc = FALSE)

# Save results
save_results(results_fixed, output_dir, "test_fixed_mcmc")

cat("Fixed dimension MCMC completed.\n")
cat("Inferred edges:", sum(results_fixed$h), "\n")
cat("True edges:", sum(synthetic_data$true_partial_order), "\n")
cat("Acceptance rate:", round(results_fixed$overall_acceptance_rate * 100, 1), "%\n\n")

# Test 3: Reversible Jump MCMC Inference
cat("3. Testing Reversible Jump MCMC Inference\n")
cat("------------------------------------------\n")

# Update config for RJ-MCMC
mcmc_config$mcmc$update_probabilities$K <- 0.1
mcmc_config$mcmc$update_probabilities$U <- 0.3

cat("Running reversible jump MCMC inference...\n")
results_rj <- run_inference(synthetic_data, mcmc_config, use_rj_mcmc = TRUE)

# Save results
save_results(results_rj, output_dir, "test_rj_mcmc")

cat("Reversible jump MCMC completed.\n")
cat("Inferred edges:", sum(results_rj$h), "\n")
cat("True edges:", sum(synthetic_data$true_partial_order), "\n")
cat("Acceptance rate:", round(results_rj$overall_acceptance_rate * 100, 1), "%\n")

if (!is.null(results_rj$K_trace)) {
  final_K <- results_rj$K_trace[length(results_rj$K_trace)]
  cat("Final K:", final_K, "\n")
}

cat("\n")

# Test 4: Generate Plots and Analysis
cat("4. Testing Plot Generation and Analysis\n")
cat("---------------------------------------\n")

cat("Generating plots and analysis for fixed dimension MCMC...\n")
generate_plots(results_fixed, synthetic_data, mcmc_config, output_dir, "test_fixed_mcmc")

cat("Generating plots and analysis for reversible jump MCMC...\n")
generate_plots(results_rj, synthetic_data, mcmc_config, output_dir, "test_rj_mcmc")

# Test 5: Compare Results
cat("\n5. Comparing Results\n")
cat("--------------------\n")

true_po <- if (is.list(synthetic_data$true_partial_order)) {
  do.call(rbind, synthetic_data$true_partial_order)
} else {
  as.matrix(synthetic_data$true_partial_order)
}

# Calculate accuracies
accuracy_fixed <- sum(results_fixed$h == true_po) / length(true_po)
accuracy_rj <- sum(results_rj$h == true_po) / length(true_po)

cat("True partial order:\n")
print(true_po)
cat("\nFixed dimension inferred partial order:\n")
print(results_fixed$h)
cat("\nReversible jump inferred partial order:\n")
print(results_rj$h)

cat("\nAccuracy Comparison:\n")
cat("Fixed dimension MCMC accuracy:", round(accuracy_fixed * 100, 1), "%\n")
cat("Reversible jump MCMC accuracy:", round(accuracy_rj * 100, 1), "%\n")

cat("\nEdge Count Comparison:\n")
cat("True edges:", sum(true_po), "\n")
cat("Fixed dimension inferred edges:", sum(results_fixed$h), "\n")
cat("Reversible jump inferred edges:", sum(results_rj$h), "\n")

cat("\n=== All Tests Completed Successfully ===\n")
cat("Results saved in:", output_dir, "\n")
cat("Files generated:\n")
cat("- test_synthetic_data.json\n")
cat("- test_fixed_mcmc_results.json\n")
cat("- test_fixed_mcmc_partial_order.rds\n")
cat("- test_rj_mcmc_results.json\n")
cat("- test_rj_mcmc_partial_order.rds\n") 