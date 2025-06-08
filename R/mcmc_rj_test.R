# 1. Load libraries and source files
library(jsonlite)
library(mvtnorm)  # Ensure this is installed

source("R/utilities.R")
source("R/mcmc_rj.R")
source("R/analysis.R")
source("R/po_inference.R")

# 2. Load data
cat("=== load data ===\n\n")
data_rj <- jsonlite::fromJSON("results/partial_order_data_rj.json")
cat("Data loading completed!\n\n")

# 3. Prepare configuration
rj_cfg <- list(
  iters = 10000,
  mcmc_pt = c(rho = 0.1, noise = 0.2, U = 0.3, beta = 0.2, K = 0.2),
  dr = 0.8,
  drbeta = 0.8,
  sigma_mal = 0.1,
  sigma_b = 0.5,
  K_prior = 3.0,
  noise_option = "queue_jump",
  rho_prior = 0.1667,
  noise_beta_prior = 9,
  mallow_ua = 1
)

# 4. Prepare arguments
args_common <- list(
  observed_orders = data_rj$observed_orders,
  choice_sets = data_rj$choice_sets,
  num_iterations = rj_cfg$iters,
  X = data_rj$X,
  dr = rj_cfg$dr,
  drbeta = rj_cfg$drbeta,
  sigma_mallow = rj_cfg$sigma_mal,
  sigma_beta = rj_cfg$sigma_b,
  mcmc_pt = rj_cfg$mcmc_pt,  # Fixed from mcmc_pt_rj
  rho_prior = rj_cfg$rho_prior,
  K_prior = rj_cfg$K_prior,
  random_seed = 42,
  noise_option = rj_cfg$noise_option
)

# Add appropriate noise prior
if (rj_cfg$noise_option == "queue_jump") {
  args_common$noise_beta_prior <- rj_cfg$noise_beta_prior
} else {
  args_common$mallow_ua <- rj_cfg$mallow_ua
}

# 5. Run MCMC
rj_results <- do.call(mcmc_partial_order_k, args_common)



# Plot log likelihood
if (length(rj_results$log_likelihood_currents) > 0) {
  plot(rj_results$log_likelihood_currents, type = "l",
       main = "Log-Likelihood Trace", xlab = "Iteration", ylab = "Log-Likelihood")
}


# 6. Diagnostics and plotting
burn_in <- 2000
burn_in_index <- burn_in / 100  # Since traces are stored every 100 iterations
# Prepare true parameters (adjust according to your actual data structure)
rj_true_params <- list(
  rho_true = data_rj$rho_true,  # Adjust if needed
  prob_noise_true = data_rj$prob_noise_true,  # Adjust if needed
  K_true = data_rj$K_true  # Adjust if needed


  
plot_mcmc_results(rj_results,
                  true_param = rj_true_params,
                  config = list(
                    prior = list(
                      rho_prior = rj_cfg$rho_prior,
                      noise_beta_prior = rj_cfg$noise_beta_prior,
                      K_prior = rj_cfg$K_prior
                    )
                  ),
                  burn_in = 20)  # More reasonable burn-in value


# 7. Estimate final partial order
threshold <- 0.5
h_trace_rj_post_burnin <- rj_results$h_trace[(burn_in_index + 1):length(rj_results$h_trace)]

if (length(h_trace_rj_post_burnin) > 0 && exists("transitive_reduction")) {
  h_array_rj <- array(unlist(h_trace_rj_post_burnin), 
                      dim = c(nrow(h_trace_rj_post_burnin[[1]]), 
                              ncol(h_trace_rj_post_burnin[[1]]), 
                              length(h_trace_rj_post_burnin)))
  
  mean_h_rj <- apply(h_array_rj, c(1,2), mean)
  final_h_rj <- (mean_h_rj >= threshold) * 1
  final_h_rj <- transitive_reduction(final_h_rj)
} else {
  final_h_rj <- matrix(0, nrow = nrow(data_rj$true_partial_order), 
                       ncol = ncol(data_rj$true_partial_order))
}

# 8. Plot comparison
if (exists("plot_partial_order")) {
  par(mfrow = c(1, 2))
  plot_partial_order(
    data_rj$true_partial_order, data_rj$items,
    title = "True Partial Order",
    vertex_size = 38,
    edge_arrow_size = 0.6,
    label_cex = 1.3,
    frame_width = 2.3
  )
  plot_partial_order(
    final_h_rj, data_rj$items, 
    title = "RJ-MCMC Estimate\n(Posterior Mean)",
    vertex_size = 38,
    edge_arrow_size = 0.6,
    label_cex = 1.3,
    frame_width = 2.3
  )
}

