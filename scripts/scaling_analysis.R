#!/usr/bin/env Rscript

# =============================================================================
# Scaling Analysis for MCMC Partial Order Inference
# =============================================================================
# This script evaluates computational performance and accuracy across different
# numbers of nodes (items) in the partial order.

# Load required packages and functions
source("scripts/install_packages.R")
source("R/utilities.R")
source("R/mcmc.R")
source("R/mcmc_rj.R")
source("R/analysis.R")

library(mvtnorm)
library(ggplot2)
library(igraph)
library(graph)

# =============================================================================
# Configuration
# =============================================================================

# Node sizes to test
node_sizes <- c(5, 6, 7, 8, 9, 10, 15, 20, 25)

# Fixed parameters across all runs
N <- 30                    # number of total orders
K <- 3                     # latent dimensions (for fixed MCMC)
p <- 2                     # number of covariates
min_sub <- 3               # minimum subset size
rho_prior <- 0.16667       # prior mean for rho
prob_noise_true <- 0.10    # queue-jump noise
beta_true <- c(0.5, -0.3)  # regression coefficients
rho_true <- 0.9            # true correlation
rng_seed <- 42             # reproducibility
noise_option <- "queue_jump"

# MCMC Configuration
mcmc_config <- list(
  iters = 50000,             # Reduced for faster runs across many sizes
  mcmc_pt = c(rho = 0.20, noise = 0.20, U = 0.40, beta = 0.20),
  dr = 0.10,
  drbeta = 0.10,
  sigma_mal = 0.10,
  sigma_b = 1.00,
  noise_beta_prior = 9
)

# RJ-MCMC Configuration  
rj_config <- list(
  iters = 50000,             # Reduced for faster runs
  mcmc_pt = c(rho = 0.15, noise = 0.15, U = 0.30, beta = 0.15, K = 0.25),
  dr = 0.10,
  drbeta = 0.10,
  sigma_mal = 0.10,
  sigma_b = 1.00,
  K_prior = 3.0,
  noise_beta_prior = 9
)

burn_in <- 5000            # Burn-in period
burn_in_index <- burn_in / 100

# =============================================================================
# Helper Functions
# =============================================================================

# No accuracy computation needed - focus on timing only

# =============================================================================
# Main Scaling Analysis
# =============================================================================

# Initialize results storage - focusing on timing and performance metrics only
results <- data.frame(
  n_nodes = integer(),
  data_gen_time = numeric(),
  fixed_mcmc_time = numeric(),
  rj_mcmc_time = numeric(),
  total_time = numeric(),
  fixed_acceptance_rate = numeric(),
  rj_acceptance_rate = numeric(),
  rj_mean_K = numeric(),
  stringsAsFactors = FALSE
)

cat("Starting scaling analysis for node sizes:", paste(node_sizes, collapse = ", "), "\n")
cat("Each run includes data generation + Fixed MCMC + RJ-MCMC\n")
cat("Progress:\n")

for (i in seq_along(node_sizes)) {
  n <- node_sizes[i]
  cat(sprintf("[%d/%d] Processing n = %d nodes...\n", i, length(node_sizes), n))
  
  set.seed(rng_seed)  # Ensure reproducibility
  
  # -------------------------------------------------------------------------
  # 1. Data Generation
  # -------------------------------------------------------------------------
  cat("  - Generating data...")
  data_gen_start <- Sys.time()
  
  # Generate design matrix
  X <- matrix(rnorm(p * n), nrow = p, ncol = n)
  
  # Generate synthetic data
  synthetic_data <- generate_synthetic_data(
    n_items = n,
    n_observations = N,
    min_sub = min_sub,
    K = K,
    rho_true = rho_true,
    prob_noise_true = prob_noise_true,
    beta_true = beta_true,
    X = t(X),
    random_seed = rng_seed
  )
  
  data_gen_end <- Sys.time()
  data_gen_time <- as.numeric(difftime(data_gen_end, data_gen_start, units = "secs"))
  cat(sprintf(" %.2f sec\n", data_gen_time))
  
  # Extract components
  items <- synthetic_data$items
  observed_orders <- synthetic_data$observed_orders
  choice_sets <- synthetic_data$choice_sets
  h_true <- synthetic_data$h_true
  
  # -------------------------------------------------------------------------
  # 2. Fixed Dimension MCMC
  # -------------------------------------------------------------------------
  cat("  - Running Fixed MCMC...")
  fixed_mcmc_start <- Sys.time()
  
  args_fixed <- list(
    observed_orders = observed_orders,
    choice_sets = choice_sets,
    num_iterations = mcmc_config$iters,
    K = K,
    X = X,
    dr = mcmc_config$dr,
    drbeta = mcmc_config$drbeta,
    sigma_mallow = mcmc_config$sigma_mal,
    sigma_beta = mcmc_config$sigma_b,
    mcmc_pt = mcmc_config$mcmc_pt,
    rho_prior = rho_prior,
    noise_option = noise_option,
    noise_beta_prior = mcmc_config$noise_beta_prior,
    random_seed = rng_seed
  )
  
  mcmc_fixed <- do.call(mcmc_partial_order, args_fixed)
  
  fixed_mcmc_end <- Sys.time()
  fixed_mcmc_time <- as.numeric(difftime(fixed_mcmc_end, fixed_mcmc_start, units = "secs"))
  cat(sprintf(" %.2f sec\n", fixed_mcmc_time))
  
  # -------------------------------------------------------------------------
  # 3. Reversible Jump MCMC
  # -------------------------------------------------------------------------
  cat("  - Running RJ-MCMC...")
  rj_mcmc_start <- Sys.time()
  
  args_rj <- list(
    observed_orders = observed_orders,
    choice_sets = choice_sets,
    num_iterations = rj_config$iters,
    X = X,
    dr = rj_config$dr,
    drbeta = rj_config$drbeta,
    sigma_mallow = rj_config$sigma_mal,
    sigma_beta = rj_config$sigma_b,
    mcmc_pt = rj_config$mcmc_pt,
    rho_prior = rho_prior,
    K_prior = rj_config$K_prior,
    noise_option = noise_option,
    noise_beta_prior = rj_config$noise_beta_prior,
    random_seed = rng_seed
  )
  
  mcmc_rj <- do.call(mcmc_partial_order_k, args_rj)
  
  rj_mcmc_end <- Sys.time()
  rj_mcmc_time <- as.numeric(difftime(rj_mcmc_end, rj_mcmc_start, units = "secs"))
  cat(sprintf(" %.2f sec\n", rj_mcmc_time))
  
  # -------------------------------------------------------------------------
  # 4. Extract Basic Performance Metrics
  # -------------------------------------------------------------------------
  cat("  - Computing metrics...")
  
  # Extract K estimate from RJ-MCMC
  rj_K_trace_post_burnin <- mcmc_rj$K_trace[(burn_in_index + 1):length(mcmc_rj$K_trace)]
  rj_mean_K <- mean(rj_K_trace_post_burnin)
  
  # Calculate total time
  total_time <- data_gen_time + fixed_mcmc_time + rj_mcmc_time
  
  # Store results - focusing on timing and basic performance
  results[i, ] <- list(
    n_nodes = n,
    data_gen_time = data_gen_time,
    fixed_mcmc_time = fixed_mcmc_time,
    rj_mcmc_time = rj_mcmc_time,
    total_time = total_time,
    fixed_acceptance_rate = mcmc_fixed$overall_acceptance_rate,
    rj_acceptance_rate = mcmc_rj$overall_acceptance_rate,
    rj_mean_K = rj_mean_K
  )
  
  cat(" Done!\n")
  
  # Optional: Save intermediate results
  if (i %% 3 == 0) {  # Save every 3 iterations
    write.csv(results[1:i, ], file = "scaling_analysis_intermediate.csv", row.names = FALSE)
  }
}

# =============================================================================
# Save and Display Results
# =============================================================================

# Save complete results
write.csv(results, file = "scaling_analysis_results.csv", row.names = FALSE)

# Display summary table
cat("\n" + paste(rep("=", 80), collapse = "") + "\n")
cat("SCALING ANALYSIS RESULTS\n")
cat(paste(rep("=", 80), collapse = "") + "\n\n")

# Round numeric columns for display
display_results <- results
numeric_cols <- sapply(display_results, is.numeric)
display_results[numeric_cols] <- lapply(display_results[numeric_cols], function(x) round(x, 3))

print(display_results)

# =============================================================================
# Create Timing Analysis Plots
# =============================================================================

cat("\nGenerating timing analysis plots...\n")

# Timing Analysis
pdf("scaling_timing_analysis.pdf", width = 15, height = 10)
par(mfrow = c(2, 3))

# Data generation time
plot(results$n_nodes, results$data_gen_time, type = "b", pch = 19, col = "blue",
     xlab = "Number of Nodes", ylab = "Time (seconds)", 
     main = "Data Generation Time")

# Fixed MCMC time
plot(results$n_nodes, results$fixed_mcmc_time, type = "b", pch = 19, col = "red",
     xlab = "Number of Nodes", ylab = "Time (seconds)", 
     main = "Fixed MCMC Time")

# RJ-MCMC time
plot(results$n_nodes, results$rj_mcmc_time, type = "b", pch = 19, col = "green",
     xlab = "Number of Nodes", ylab = "Time (seconds)", 
     main = "RJ-MCMC Time")

# Total time
plot(results$n_nodes, results$total_time, type = "b", pch = 19, col = "purple",
     xlab = "Number of Nodes", ylab = "Time (seconds)", 
     main = "Total Runtime")

# MCMC time comparison
plot(results$n_nodes, results$fixed_mcmc_time, type = "b", pch = 19, col = "red",
     xlab = "Number of Nodes", ylab = "Time (seconds)", 
     main = "MCMC Time Comparison", ylim = range(c(results$fixed_mcmc_time, results$rj_mcmc_time)))
lines(results$n_nodes, results$rj_mcmc_time, type = "b", pch = 17, col = "green")
legend("topleft", legend = c("Fixed MCMC", "RJ-MCMC"), 
       col = c("red", "green"), pch = c(19, 17), lty = 1)

# Dimension selection (RJ-MCMC)
plot(results$n_nodes, results$rj_mean_K, type = "b", pch = 17, col = "purple",
     xlab = "Number of Nodes", ylab = "Estimated K", 
     main = "RJ-MCMC Dimension Selection")
abline(h = K, col = "black", lty = 2)
legend("topright", legend = c("Estimated K", "True K"), 
       col = c("purple", "black"), lty = c(1, 2), pch = c(17, NA))

dev.off()

cat("Analysis complete! Results saved to:\n")
cat("- scaling_analysis_results.csv\n")
cat("- scaling_timing_analysis.pdf\n") 