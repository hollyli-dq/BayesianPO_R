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
  iters            = 100000,                      # total MH iterations
  K=3.0,
  random_seed =42,
  mcmc_pt = c(
    rho   = 0.125,   # 0.10 / 0.80
    noise = 0.250,   # 0.20 / 0.80
    U     = 0.375,   # 0.30 / 0.80
    beta  = 0.250    # 0.20 / 0.80
  ),
  dr               = 0.10,                      # rho step
  drbeta           = 0.80,                      # beta step
  sigma_mal        = 0.10,                      # σ for Mallows θ proposal
  sigma_b          = 0.5,                      # prior sd for β
  noise_option     = "queue_jump",              # or "mallows_noise"
  rho_prior        = 0.16667,
  noise_beta_prior = 9,
  random_seed =42# queue-jump prior (omit for Mallows)
  
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


results_fixed <- do.call(mcmc_partial_order, args_fix)

rj_cfg <- list(
  iters = 100000,
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
args_common_rj <- list(
  observed_orders = data$observed_orders,
  choice_sets = data$choice_sets,
  num_iterations = rj_cfg$iters,
  X =data$X,
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
  args_common_rj$noise_beta_prior <- rj_cfg$noise_beta_prior
} else {
  args_common_rj$mallow_ua <- rj_cfg$mallow_ua
}

rj_results <- do.call(mcmc_partial_order_k, args_common_rj)

# 4.  Run the sampler ---------------------------------------------------------


cat(
  "Fixed-dimension MCMC completed — acceptance rate:",
  sprintf("%.2f %%", 100 * results_fixed$overall_acceptance_rate), "\n"
)


save_results(results_fixed, "results", "test_fixed")
save_results(rj_results, "results", "test_rj")


mcmc_results_python <-  jsonlite::fromJSON("results/mcmc_results.json")
mcmc_rj_results_python <-  jsonlite::fromJSON("results/mcmc_results_rj_k.json")


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


# Plot log likelihood
if (length(rj_results$log_likelihood_currents) > 0) {
  plot(rj_results$log_likelihood_currents, type = "l",
       main = "Log-Likelihood Trace(Queue Jump)", xlab = "Iteration", ylab = "Log-Likelihood")
}



burn_in<-20000
burn_in_index <- burn_in / 100 
burn_indices <- (burn_in/100 + 1):length(results_fixed$rho_trace)

true_params <- list(
  rho_true =  data$parameter$rho_true,
  prob_noise_true = data$parameter$prob_noise_true,
  beta_true= data$beta_true. 
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

plot_mcmc_results_dual(results_fixed,
                       mcmc_results_python,
                       true_param = list(rho_true = 0.7),
                       burn_in = 500)


plot_mcmc_results_dual(results_fixed,
                       mcmc_results_python,
                       true_param = list(rho_true = 0.7),
                       burn_in = 500)
par(mfrow = c(1, 2))
beta_true<-as.numeric((data$beta_true))
plot_beta_parameters_dual(results_fixed,
                          mcmc_results_python,
                                      true_param  =beta_true,
                                      burn_in     = 500,
                                      out_dir     = ".",
                                      prefix      = "beta")

  
  
  
threshold <- 0.5
# Extract post-burn-in partial order traces
h_trace_fixed_post_burnin <- results_fixed$h_trace[burn_indices]
h_trace_kj_post_burnin <- rj_results$h_trace[burn_indices]
h_trace_fixed_post_burnin_py <- mcmc_results_python$h_trace[burn_indices]
h_trace_kj_post_burnin_py <- mcmc_rj_results_python$h_trace[burn_indices]

threshold <- 0.5
cat("Computing posterior-mean partial orders …\n")

###############################################################################
# 1.  Fixed-K   (R implementation)   →  final_h_fixed
###############################################################################
if (length(h_trace_fixed_post_burnin) > 0) {
  h_array_fixed <- array(
    unlist(h_trace_fixed_post_burnin, use.names = FALSE),
    dim = c(nrow(h_trace_fixed_post_burnin[[1]]),
            ncol(h_trace_fixed_post_burnin[[1]]),
            length(h_trace_fixed_post_burnin))
  )
  mean_h_fixed  <- apply(h_array_fixed, c(1, 2), mean)
  final_h_fixed <- (mean_h_fixed >= threshold) * 1L
  final_h_fixed <- transitive_reduction(final_h_fixed)
  
  cat("  • Fixed-K (R):", length(h_trace_fixed_post_burnin),
      "post-burn-in samples averaged\n")
} else {
  final_h_fixed <- matrix(0, nrow = nrow(h_true), ncol = ncol(h_true))
  cat("  • Warning: no post-burn-in samples for Fixed-K (R)\n")
}

###############################################################################
# 2.  RJ-K      (R implementation)   →  final_h_kj
###############################################################################
if (length(h_trace_kj_post_burnin) > 0) {
  h_array_kj <- array(
    unlist(h_trace_kj_post_burnin, use.names = FALSE),
    dim = c(nrow(h_trace_kj_post_burnin[[1]]),
            ncol(h_trace_kj_post_burnin[[1]]),
            length(h_trace_kj_post_burnin))
  )
  mean_h_kj  <- apply(h_array_kj, c(1, 2), mean)
  final_h_kj <- (mean_h_kj >= threshold) * 1L
  final_h_kj <- transitive_reduction(final_h_kj)
  
  cat("  • RJ-K (R):", length(h_trace_kj_post_burnin),
      "post-burn-in samples averaged\n")
} else {
  final_h_kj <- matrix(0, nrow = nrow(h_true), ncol = ncol(h_true))
  cat("  • Warning: no post-burn-in samples for RJ-K (R)\n")
}

###############################################################################
# 3.  Fixed-K   (Python trace)       →  final_h_fixed_py
###############################################################################
library(abind)            # install.packages("abind")  once

to_array3d <- function(trace_list) {
  mats <- lapply(trace_list, function(x) {
    m <- as.matrix(x)                 # vector → 1 × p  ; matrix stays matrix
    if (length(dim(m)) != 2L) stop("each h must be 2-D")
    m
  })
  r <- nrow(mats[[1]]);  c <- ncol(mats[[1]])
  if (!all(vapply(mats, function(m) all(dim(m) == c(r, c)), FALSE)))
    stop("not all h-matrices have the same shape")
  abind(mats, along = 3)              # r × c × n_iter
}

###############################################################################
# 3.  Fixed-K   (Python trace)   →  final_h_fixed_py
###############################################################################
library(abind)            # install.packages("abind")  once

to_array3d <- function(trace_list) {
  mats <- lapply(trace_list, function(x) {
    m <- as.matrix(x)                 # vector → 1 × p  ; matrix stays matrix
    if (length(dim(m)) != 2L) stop("each h must be 2-D")
    m
  })
  r <- nrow(mats[[1]]);  c <- ncol(mats[[1]])
  if (!all(vapply(mats, function(m) all(dim(m) == c(r, c)), FALSE)))
    stop("not all h-matrices have the same shape")
  abind(mats, along = 3)              # r × c × n_iter
}

if (length(h_trace_fixed_post_burnin_py) > 0) {
  
  h_array_fixed_py <- to_array3d(h_trace_fixed_post_burnin_py)
  mean_h_fixed_py  <- apply(h_array_fixed_py, c(1, 2), mean)
  
  final_h_fixed_py <- transitive_reduction( (mean_h_fixed_py >= threshold) * 1L )
  
  cat("  • Fixed-K (Py):", dim(h_array_fixed_py)[3],
      "post-burn-in samples averaged\n")
  
} else {
  
  final_h_fixed_py <- matrix(0, nrow = nrow(h_true), ncol = ncol(h_true))
  cat("  • Warning: no post-burn-in samples for Fixed-K (Py)\n")
  
}
###############################################################################
# 4.  RJ-K      (Python trace)       →  final_h_kj_py
###############################################################################

## 4. RJ-K  (Python trace)  → final_h_kj_py ───────────────────────────────
if (length(h_trace_kj_post_burnin_py) > 0) {
  h_array_kj_py <- to_array3d(h_trace_kj_post_burnin_py)     # r×c×n
  mean_h_kj_py  <- apply(h_array_kj_py, c(1, 2), mean)       # r×c
  final_h_kj_py <- transitive_reduction((mean_h_kj_py >= threshold) * 1L)
  
  cat("  • RJ-K (Py):", dim(h_array_kj_py)[3], "post-burn-in samples averaged\n")
} else {
  final_h_kj_py <- matrix(0, nrow = nrow(h_true), ncol = ncol(h_true))
  cat("  • Warning: no post-burn-in samples for RJ-K (Py)\n")
}
###############################################################################
# At this point you have:
#   final_h_fixed      – fixed-K,   R trace
#   final_h_kj         – RJ-K,      R trace
#   final_h_fixed_py   – fixed-K,   Python trace
#   final_h_kj_py      – RJ-K,      Python trace
###############################################################################

op <- par(mfrow = c(2, 3), mar = c(2, 2, 3, 1))  # save old par; tweak margins

plot_partial_order(
  transitive_reduction(h_true),  items,
  title           = "True Partial Order",
  vertex_size     = 38,
  edge_arrow_size = 0.6,
  label_cex       = 1.3,
  frame_width     = 2.3
)

plot_partial_order(
  final_h_fixed,  items,
  title           = "Est. Fixed-K (R)",
  vertex_size     = 38,
  edge_arrow_size = 0.6,
  label_cex       = 1.3,
  frame_width     = 2.3
)

plot_partial_order(
  final_h_kj,     items,
  title           = "Est. RJ-K (R)",
  vertex_size     = 38,
  edge_arrow_size = 0.6,
  label_cex       = 1.3,
  frame_width     = 2.3
)

plot_partial_order(
  final_h_fixed_py, items,
  title           = "Est. Fixed-K (Py)",
  vertex_size     = 38,
  edge_arrow_size = 0.6,
  label_cex       = 1.3,
  frame_width     = 2.3
)

plot_partial_order(
  final_h_kj_py,    items,
  title           = "Est. RJ-K (Py)",
  vertex_size     = 38,
  edge_arrow_size = 0.6,
  label_cex       = 1.3,
  frame_width     = 2.3
)

par(op)  # restore the user’s old graphics settings