# Bayesian Partial Order Inference - MCMC Functions
# Proper likelihood computation and MCMC simulation matching Python implementation
# Author: Converted for Prof Geoff Nicholls, University of Oxford

# ============================================================================
# MCMC SIMULATION (PROPER IMPLEMENTATION)
# ============================================================================
# Create a cache environment at package level (persistent across calls)
.likelihood_cache <- new.env(hash = TRUE)

# Helper function to create cache keys
create_cache_key <- function(adj_matrix) {
  digest::digest(adj_matrix)
}

# Cache-aware number of linear extensions calculation
get_nle <- function(adj_matrix) {
  key <- create_cache_key(adj_matrix)
  
  if (exists(key, envir = .likelihood_cache)) {
    return(get(key, envir = .likelihood_cache))
  }
  
  val <- nle(adj_matrix)
  assign(key, val, envir = .likelihood_cache)
  return(val)
}

# Cache-aware number of extensions with first item calculation
get_nle_first <- function(adj_matrix, local_idx) {
  key <- paste0(create_cache_key(adj_matrix), "_", local_idx)
  
  if (exists(key, envir = .likelihood_cache)) {
    return(get(key, envir = .likelihood_cache))
  }
  
  val <-  num_extensions_with_first(adj_matrix, local_idx)
  assign(key, val, envir = .likelihood_cache)
  return(val)
}



# ============================================================================
# LIKELIHOOD COMPUTATION (PROPER IMPLEMENTATION)
# ============================================================================


#' Calculate log likelihood for queue jump noise model (proper implementation)
log_likelihood_queue_jump <- function(h, observed_orders_idx, choice_sets, item_to_index, prob_noise) {
  log_likelihood <- 0.0
  
  for (idx in seq_along(observed_orders_idx)) {
    y_i <- observed_orders_idx[[idx]]
    O_i <- choice_sets[[idx]]
    O_i_indices <- sort(sapply(O_i, function(item) item_to_index[[as.character(item)]]))
    m <- length(y_i)
    
    # Process each position in the sequence
    for (j in 1:(m-1)) {  # Only go to m-1 like Python
      y_j <- y_i[j]
      remaining_indices <- y_i[j:length(y_i)]
      
      # Extract submatrix for remaining items
      h_Z_remaining <- h[remaining_indices, remaining_indices, drop = FALSE]
      
      # Apply transitive reduction during likelihood calculation (like Python)
      tr_remaining <- transitive_reduction(h_Z_remaining)
      
      # Count total linear extensions
      num_le <- get_nle(tr_remaining)
      
      # Count extensions where y_j is first
      local_idx <- which(remaining_indices == y_j)
      if (length(local_idx) > 0 && num_le > 0) {
        # Use the proper algorithm to count extensions starting with y_j
        num_first_item <- get_nle_first(tr_remaining, local_idx)
      } else {
        num_first_item <- 0
      }
      
      # Calculate probabilities (matching Python logic)
      if (num_le > 0) {
        prob_no_jump <- (1 - prob_noise) * (num_first_item / num_le)
      } else {
        prob_no_jump <- 0
      }
      prob_jump <- prob_noise * (1 / (m - j + 1))
      prob_observed <- prob_no_jump + prob_jump
      
      # Ensure positive probability
      prob_observed <- max(prob_observed, 1e-20)
      log_likelihood <- log_likelihood + log(prob_observed)
    }
  }
  
  return(log_likelihood)
}


mcmc_partial_order <- function(observed_orders, 
                               choice_sets,
                               num_iterations, 
                               K, 
                               dr = 0.1,
                               rho_prior = 1/6, 
                               noise_option = "queue_jump",
                               noise_beta_prior = NULL,  
                               random_seed = 123) {
  
  set.seed(random_seed)
 
  # ----------------------------------------------------------------
  # 1. Setup: Map items to indices, initialize states, etc.
  # ----------------------------------------------------------------
  
  # Setup: Map items to indices
  items <- sort(unique(unlist(choice_sets)))
  n <- length(items)
  item_to_index <- setNames(seq_along(items), items)
  index_to_item <- setNames(items, seq_along(items))
  
  # Convert observed orders to index form
  observed_orders_idx <- lapply(observed_orders, function(order) {
    sapply(order, function(item) item_to_index[[as.character(item)]])
  })
  
  # ----------------------------------------------------------------
  # 2. Prepare Storage for MCMC results
  # ----------------------------------------------------------------
  
  # Initialize MCMC state  
  Z <- rmvnorm(n, mean = rep(0, K), sigma = diag(K))
  eta <- Z
  h_Z <- generate_partial_order(eta)
  
  # Initialize parameters
  rho <- sample_rho_prior(rho_prior)
  prob_noise <- sample_noise_prior(noise_beta_prior)
  
  # Storage for results
  Z_trace <- list()
  rho_trace <- numeric()
  prob_noise_trace <- numeric()
  h_trace <- list()
  log_likelihood_trace <- numeric()
  acceptance_rates <- data.frame(list(rho=0,prob_noise=0,Z=0))
  
  # Progress intervals
  progress_intervals <- seq(0.1, 1.0, 0.1) * num_iterations
  
  cat("Starting MCMC with", num_iterations, "iterations...\n")
  cat("Initial h_Z has", sum(h_Z), "edges\n") # Main MCMC loop
  
  llk_current <- log_likelihood_queue_jump(h_Z, observed_orders_idx, choice_sets, item_to_index, prob_noise)
  
  # ----------------------------------------------------------------
  # 3. Main MCMC Loop (MATCH PYTHON EXACTLY)
  # ----------------------------------------------------------------
  
  for (iteration in 1:num_iterations) {
    
    # Update rho
    delta <- runif(1, dr, 1/dr)  # Correct range
    rho_prime <- 1.0 - (1.0 - rho) * delta
    
    if (!is.na(rho_prime) && rho_prime > 0 && rho_prime < 1) {
      log_prior_current <- log_rho_prior(rho, rho_prior) + log_U_prior(Z, rho, K)
      log_prior_proposed <- log_rho_prior(rho_prime, rho_prior) + log_U_prior(Z, rho_prime, K)
      
      log_acceptance_ratio <- (log_prior_proposed) - (log_prior_current)- log(delta)
      
      acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
      if (runif(1) < acceptance_probability) {
        rho <- rho_prime
        acceptance_rates$rho <- acceptance_rates$rho + 1
      }
    }
    
    # Update noise parameter
    prob_noise_prime <- sample_noise_prior(noise_beta_prior)
    
    log_prior_current <- log_noise_prior(prob_noise, noise_beta_prior)
    log_prior_proposed <- log_noise_prior(prob_noise_prime, noise_beta_prior)
    
    llk_prime <- log_likelihood_queue_jump(h_Z, observed_orders_idx, choice_sets, item_to_index, prob_noise_prime)
    
    log_acceptance_ratio <- llk_prime- llk_current
    acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
    if (runif(1) < acceptance_probability) {
      prob_noise <- prob_noise_prime
      llk_current <- llk_prime
      acceptance_rates$prob_noise <- acceptance_rates$prob_noise + 1
    }
    
    
    # Update U (latent matrix Z) 
    i <- sample(1:n, 1)
    
    Z_prime <- Z
    Sigma_row <- build_sigma_rho(K, rho)
    Z_prime[i,] <- mvrnorm(1, Z[i,], Sigma_row)
    
    eta_prime <- Z_prime 
    h_Z_prime <- generate_partial_order(eta_prime)
    
    log_prior_current <- log_U_prior(Z, rho, K)
    log_prior_proposed <- log_U_prior(Z_prime, rho, K)
    
    llk_prime <- log_likelihood_queue_jump(h_Z_prime, observed_orders_idx, choice_sets, item_to_index, prob_noise)
    
    log_acceptance_ratio <- (log_prior_proposed + llk_prime) - (log_prior_current + llk_current)
    acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
    if (runif(1) < acceptance_probability) {
      Z <- Z_prime
      h_Z <- h_Z_prime
      llk_current <- llk_prime
      acceptance_rates$Z <- acceptance_rates$Z + 1
    }
    
    
    # Store current state every 10 iterations
    if (iteration %% 10 == 0) {
      Z_trace <- append(Z_trace, list(Z))
      h_trace <- append(h_trace, list(h_Z))
      rho_trace <- c(rho_trace, rho)
      prob_noise_trace <- c(prob_noise_trace, prob_noise)
      log_likelihood_trace <- c(log_likelihood_trace,llk_current)
    }
    
    if (iteration %in% progress_intervals) {
      cat(sprintf("Iteration %d/%d - Accept Rate: %.2f%% - Edges: %d\n", 
                  iteration, num_iterations, median(unlist(acceptance_rates/iteration * 100)), sum(h_Z)))
    }
  }
  
  acceptance_rates=acceptance_rates/num_iterations  
  
  return(list(
    Z_trace = Z_trace,
    h_trace = h_trace,
    index_to_item = index_to_item,
    item_to_index = item_to_index,
    rho_trace = rho_trace,
    prob_noise_trace = prob_noise_trace,
    acceptance_rates = acceptance_rates,
    log_likelihood_trace = log_likelihood_trace,
    final_h = h_Z,
    final_Z = Z,
    final_rho = rho,
    final_prob_noise = prob_noise
  ))
} 
