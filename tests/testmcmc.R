# Simple CLI test bypassing YAML issues
library(jsonlite)

# Source the functions
source("R/utilities.R")
source("R/mcmc.R")
source("R/mcmc_rj.R")
source("R/analysis.R")
source("R/po_inference.R")

cat("=== load data ===\n\n")

# Load the JSON data
data <- jsonlite::fromJSON("results/partial_order_data.json")


# Print completion message
cat("Data loading completed!\n\n")

# Optional: Print structure of loaded data
cat("Data structure:\n")
str(data, max.level = 2)


h_true<-data$true_partial_order
items<-data$items

# Test 2: Fixed Dimension MCMC
cat("2. Fixed Dimension MCMC\n")
cat("-----------------------\n")

# config.R


fixed_cfg <- list(
  iters            = 10000,                      # total MH iterations
  K=3.0,
  random_seed =42,
  mcmc_pt          = c(rho  = 0.10,             # update probabilities
                       noise = 0.20,
                       U     = 0.30,
                       beta  = 0.20),
  dr               = 0.10,                      # rho step
  drbeta           = 0.80,                      # beta step
  sigma_mal        = 0.10,                      # σ for Mallows θ proposal
  sigma_b          = 0.5,                      # prior sd for β
  noise_option     = "queue_jump",              # or "mallows_noise"
  rho_prior        = 0.16667,
  noise_beta_prior = 9                        # queue-jump prior (omit for Mallows)
  
)

# 2.  Common arguments --------------------------------------------------------
args_fix <- list(
  observed_orders = data$observed_orders,
  choice_sets     = data$choice_sets,
  num_iterations  = fixed_cfg$iters,
  K               = fixed_cfg$K,
  X               = data$X,
  dr              = fixed_cfg$dr,
  drbeta          = fixed_cfg$drbeta,
  sigma_mallow    = fixed_cfg$sigma_mal, 
  sigma_beta      = fixed_cfg$sigma_b,
  mcmc_pt         = fixed_cfg$mcmc_pt,
  rho_prior       = fixed_cfg$rho_prior,
  random_seed     = fixed_cfg$rng_seed
)

# 3.  Attach the correct noise-specific prior ---------------------------------
if (fixed_cfg$noise_option == "queue_jump") {
  args_fix$noise_option     <- "queue_jump"
  args_fix$noise_beta_prior <- fixed_cfg$noise_beta_prior   # keep
} else {  # "mallows_noise"
  args_fix$noise_option <- "mallows_noise"
  args_fix$mallow_ua    <- fixed_cfg$mallow_ua              # keep
}

# 4.  Run the sampler ---------------------------------------------------------
results_fixed <- do.call(mcmc_partial_order, args_fix)

cat(
  "Fixed-dimension MCMC completed — acceptance rate:",
  sprintf("%.2f %%", 100 * results_fixed$overall_acceptance_rate), "\n"
)


save_results(results_fixed, "results", "cli_test_fixed")

cat("Fixed MCMC completed!\n")
cat("Acceptance rate:", round(results_fixed$overall_acceptance_rate * 100, 1), "%\n\n")


par(mfrow = c(1, 1))

# Check if we have valid data before plotting
if (length(results_fixed $log_likelihood_currents) > 0 && 
    all(is.finite(results_fixed $log_likelihood_currents))) {
  plot(results_fixed$log_likelihood_currents, type = "l", 
       main = "Log-Likelihood (Fixed K)", xlab = "Iteration", ylab = "Log-Likelihood")
} else {
  plot(1, type = "n", main = "Log-Likelihood (Fixed K) - No Valid Data", 
       xlab = "Iteration", ylab = "Log-Likelihood")
}


burn_in<-4000
burn_in_index <- burn_in / 100 
burn_indices <- (burn_in/100 + 1):length(results_fixed$rho_trace)

# Fixed dimension estimates
rho_est_fixed <- mean(results_fixed$rho_trace[burn_indices])
beta_est_fixed <- colMeans(do.call(rbind, results_fixed$beta_trace[burn_indices]))


true_params <- list(
  rho_true =  data$parameter$rho_true,
  prob_noise_true = data$parameter$prob_noise_true,
  beta_true= data$parameters$beta. 
)

config <- list(
  prior = list(
    rho_prior = args_fix$prior$rho_prior,
    noise_beta_prior =args_fix $prior$noise_beta_prior
  ),
  noise = list(
    noise_option = "queue_jump"
  )
)

plot_mcmc_results(results_fixed, 
                  true_param = true_params,
                  config = config,
                  burn_in = 2)



threshold <- 0.5
# Extract post-burn-in partial order traces
h_trace_fixed_post_burnin <- results_fixed$h_trace[burn_indices]

# Compute posterior mean partial orders
cat("Computing posterior mean partial orders...\n")

# Fixed K MCMC: Average across all post-burn-in iterations
if (length(h_trace_fixed_post_burnin) > 0) {
  # Convert list of matrices to 3D array for easier averaging
  h_array_fixed <- array(0, dim = c(nrow(h_trace_fixed_post_burnin[[1]]), 
                                    ncol(h_trace_fixed_post_burnin[[1]]), 
                                    length(h_trace_fixed_post_burnin)))
  
  for (i in seq_along(h_trace_fixed_post_burnin)) {
    h_array_fixed[,,i] <- h_trace_fixed_post_burnin[[i]]
  }
  
  # Compute mean and apply threshold
  mean_h_fixed <- apply(h_array_fixed, c(1,2), mean)
  final_h_fixed <- (mean_h_fixed >= threshold) * 1
  
  # Apply transitive reduction to get minimal representation¯
  final_h_fixed <- transitive_reduction(final_h_fixed)
  
  cat("Fixed K MCMC: Averaged", length(h_trace_fixed_post_burnin), "post-burn-in samples\n")
} else {
  final_h_fixed <- matrix(0, nrow = nrow(h_true), ncol = ncol(h_true))
  cat("Warning: No post-burn-in samples for Fixed K MCMC\n")
}


par(mfrow = c(1, 2))

plot_partial_order(
  h_true, items,
  title         = "True Partial Order",
  vertex_size   = 38,
  edge_arrow_size = 0.6,
  label_cex     = 1.3,  # even bigger text
  frame_width   = 2.3   # bolder outlines
)

plot_partial_order(
  final_h_fixed, items, title = "MCMC Estimate\n(Posterior Mean)",
  vertex_size   = 38,
  edge_arrow_size = 0.6,
  label_cex     = 1.3,  # even bigger text
  frame_width   = 2.3   # bolder outlines
)


cat("The Rho estimated is",  rho_est_fixed ,"Compared with the true rho is",data$parameter$rho_true)
cat("The beta estimated is", beta_est_fixed,"Compared with the true beta is",data$beta_true)
cat("Fixed h is",final_h_fixed, "Compared with the true h is",h_true)


final_h_fixed <- transitive_reduction(final_h_fixed)
h_true<-transitive_reduction(h_true)
h_true
