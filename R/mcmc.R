# Bayesian Partial Order Inference - MCMC Functions
# Proper likelihood computation and MCMC simulation matching Python implementation
# Author: Converted for Prof Geoff Nicholls, University of Oxford

# Required libraries
library(mvtnorm)
library(MASS)

# ============================================================================
# LIKELIHOOD COMPUTATION (PROPER IMPLEMENTATION)
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

#' Calculate log likelihood for queue jump noise model (proper implementation)
#' 
#' @param h Partial order matrix
#' @param observed_orders_idx List of observed total orders (as indices)
#' @param choice_sets List of choice sets for each observation
#' @param item_to_index Mapping from items to indices
#' @param prob_noise Noise probability
#' @return Log likelihood value
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

#' Calculate log likelihood for Mallows noise model (proper implementation)
#' 
#' @param h Partial order matrix
#' @param observed_orders_idx List of observed total orders (as indices)
#' @param choice_sets List of choice sets for each observation
#' @param item_to_index Mapping from items to indices
#' @param mallow_theta Mallows parameter
#' @return Log likelihood value
log_likelihood_mallows <- function(h, observed_orders_idx, choice_sets, item_to_index, mallow_theta) {
  log_likelihood <- 0.0
  
  for (idx in seq_along(observed_orders_idx)) {
    y_i <- observed_orders_idx[[idx]]
    O_i <- choice_sets[[idx]]
    O_i_indices <- sort(sapply(O_i, function(item) item_to_index[[as.character(item)]]))
    
    # Extract submatrix for choice set
    h_Z_Oi <- h[O_i_indices, O_i_indices, drop = FALSE]
    
    # Apply transitive reduction like Python
    tr_Oi <- transitive_reduction(h_Z_Oi)
    
    # Count inversions relative to the transitive reduction
    inversions <- 0
    n_items <- length(y_i)
    
    for (i in 1:(n_items - 1)) {
      for (j in (i + 1):n_items) {
        item1_idx <- which(O_i_indices == y_i[i])
        item2_idx <- which(O_i_indices == y_i[j])
        
        if (length(item1_idx) > 0 && length(item2_idx) > 0) {
          # If transitive reduction says item2 should come before item1, count as inversion
          if (!is.na(tr_Oi[item2_idx, item1_idx]) && tr_Oi[item2_idx, item1_idx] == 1) {
            inversions <- inversions + 1
          }
        }
      }
    }
    
    # Mallows probability with proper normalization
    mallows_prob <- exp(-mallow_theta * inversions)
    mallows_prob <- max(mallows_prob, 1e-20)
    log_likelihood <- log_likelihood + log(mallows_prob)
  }
  
  return(log_likelihood)
}

#' Main likelihood calculation function
#' 
#' @param Z Latent variables matrix
#' @param h Partial order matrix
#' @param observed_orders_idx List of observed total orders (as indices)
#' @param choice_sets List of choice sets
#' @param item_to_index Item to index mapping
#' @param prob_noise Noise probability (for queue jump)
#' @param mallow_theta Mallows parameter (for Mallows model)
#' @param noise_option Either "queue_jump" or "mallows_noise"
#' @return Log likelihood value
#' @export
calculate_log_likelihood <- function(Z, h, observed_orders_idx, choice_sets, item_to_index, 
                                   prob_noise, mallow_theta, noise_option) {
  if (noise_option == "queue_jump") {
    return(log_likelihood_queue_jump(h, observed_orders_idx, choice_sets, item_to_index, prob_noise))
  } else if (noise_option == "mallows_noise") {
    return(log_likelihood_mallows(h, observed_orders_idx, choice_sets, item_to_index, mallow_theta))
  } else {
    stop("Unknown noise option: ", noise_option)
  }
}

# ============================================================================
# MCMC SIMULATION (PROPER IMPLEMENTATION)
# ============================================================================

#' Main MCMC function for partial order inference (fixed implementation)
#' 
#' @param observed_orders List of observed total orders
#' @param choice_sets List of choice sets for each observation
#' @param num_iterations Number of MCMC iterations
#' @param K Number of latent dimensions
#' @param X Covariate matrix
#' @param dr Multiplicative step size for rho
#' @param drbeta Multiplicative step size for beta
#' @param sigma_mallow Standard deviation for Mallows proposals
#' @param sigma_beta Prior standard deviation for beta
#' @param noise_option Either "queue_jump" or "mallows_noise"
#' @param mcmc_pt Vector of update probabilities [rho, noise, U, beta]
#' @param rho_prior Prior parameter for rho
#' @param noise_beta_prior Prior parameter for noise
#' @param mallow_ua Prior parameter for Mallows
#' @param random_seed Random seed
#' @return List containing MCMC results
#' @export
mcmc_partial_order <- function(observed_orders, 
                               choice_sets,
                               num_iterations, 
                               K, 
                               X,
                              dr = 0.1, 
                              drbeta = 0.1, 
                              sigma_mallow = 0.1, 
                              sigma_beta = 1.0,
                              noise_option = "queue_jump", 
                              mcmc_pt = c(0.25, 0.25, 0.25, 0.25),
                              rho_prior = 1/6, 
                              noise_beta_prior = NULL,   # <- now nullable
                              mallow_ua        = NULL,   # <- now nullable
                              random_seed = 123) {
  noise_option <- match.arg(noise_option)
  
  ## -- enforce mutually-exclusive hyper-parameters --------------------------
  if (noise_option == "queue_jump") {
    if (is.null(noise_beta_prior))
      stop("`noise_beta_prior` must be provided when noise_option = 'queue_jump'")
    if (!is.null(mallow_ua))
      message("Ignoring `mallow_ua` because queue-jump noise is selected.")
  } else {                              # mallows_noise
    if (is.null(mallow_ua))
      stop("`mallow_ua` must be provided when noise_option = 'mallows_noise'")
    if (!is.null(noise_beta_prior))
      message("Ignoring `noise_beta_prior` because Mallows noise is selected.")
  }
  
  set.seed(random_seed)
  
  # Setup: Map items to indices
  items <- sort(unique(unlist(choice_sets)))
  n <- length(items)
  item_to_index <- setNames(seq_along(items), items)
  index_to_item <- setNames(items, seq_along(items))
  
  # Convert observed orders to index form
  observed_orders_idx <- lapply(observed_orders, function(order) {
    sapply(order, function(item) item_to_index[[as.character(item)]])
  })
  
  # Initialize MCMC state (PROPER INITIALIZATION)
  # Start with random initialization like Python, not zeros
  Z <- matrix(0, nrow = n, ncol = K)
  
  p <- nrow(X)
  beta <- rnorm(p, 0, sigma_beta)
  alpha <- as.vector(t(X) %*% beta)
  print(Z)
  print(alpha
  eta <- transform_U_to_eta(Z, alpha)
  if (anyNA(eta)) {
    stop("transform_U_to_eta() produced NA values")
  }
  h_Z <- generate_partial_order(eta)
  
  # Initialize parameters
  rho <- sample_rho_prior(rho_prior)
  prob_noise <- sample_noise_prior(noise_beta_prior)
  mallow_theta <- sample_theta_prior(mallow_ua)
  
  # Storage for results
  Z_trace <- list()
  rho_trace <- numeric()
  beta_trace <- list()
  prob_noise_trace <- numeric()
  mallow_theta_trace <- numeric()
  h_trace <- list()
  log_likelihood_currents <- numeric()
  acceptance_rates <- numeric()
  num_acceptances <- 0
  
  # Progress intervals
  progress_intervals <- seq(0.1, 1.0, 0.1) * num_iterations
  
  # Unpack update probabilities
  rho_pct <- mcmc_pt[1]
  noise_pct <- mcmc_pt[2]
  U_pct <- mcmc_pt[3]
  beta_pct <- mcmc_pt[4]

  cat("Starting MCMC with", num_iterations, "iterations...\n")
  cat("Initial h_Z has", sum(h_Z), "edges\n") # Main MCMC loop
  for (iteration in 1:num_iterations) {
    llk_current <- calculate_log_likelihood(Z, h_Z, observed_orders_idx, choice_sets,
                                            item_to_index, prob_noise, mallow_theta, noise_option)
    r <- runif(1)
    
    # Update rho
    if (r < rho_pct) {
      delta <- runif(1, dr, 1/dr)  # Correct range
      rho_prime <- 1.0 - (1.0 - rho) * delta
      
      if (!is.na(rho_prime) && rho_prime > 0 && rho_prime < 1) {
        log_prior_current <- log_rho_prior(rho, rho_prior) + log_U_prior(Z, rho, K)
        log_prior_proposed <- log_rho_prior(rho_prime, rho_prior) + log_U_prior(Z, rho_prime, K)
        log_likelihood_proposed <- llk_current  # Z unchanged
        
        log_acceptance_ratio <- (log_prior_proposed) - (log_prior_current)- log(delta)
        
        acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
        if (runif(1) < acceptance_probability) {
          rho <- rho_prime
          num_acceptances <- num_acceptances + 1
          llk_current <- log_likelihood_proposed
        }
      }
      

      log_likelihood_currents <- c(log_likelihood_currents, llk_current)
    }
    
    # Update noise parameter
    else if (r < (rho_pct + noise_pct)) {
      if (noise_option == "mallows_noise") {
        epsilon <- rnorm(1)
        mallow_theta_prime <- mallow_theta * exp(sigma_mallow * epsilon)
        
        log_prior_current <- log_theta_prior(mallow_theta, mallow_ua)
        log_prior_proposed <- log_theta_prior(mallow_theta_prime, mallow_ua)
        

        llk_prime <- calculate_log_likelihood(Z, h_Z, observed_orders_idx, choice_sets,
                                            item_to_index, prob_noise, mallow_theta_prime, noise_option)
        
        log_acceptance_ratio <- (log_prior_proposed + llk_prime) - (log_prior_current + llk_current) + 
                               log(mallow_theta/ mallow_theta_prime)  # Correct Jacobian
        
        acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
        if (runif(1) < acceptance_probability) {
          mallow_theta <- mallow_theta_prime
          num_acceptances <- num_acceptances + 1
          llk_current <- llk_prime
        }
      } else if (noise_option == "queue_jump") {
        prob_noise_prime <- sample_noise_prior(noise_beta_prior)
        
        log_prior_current <- log_noise_prior(prob_noise, noise_beta_prior)
        log_prior_proposed <- log_noise_prior(prob_noise_prime, noise_beta_prior)
        

        llk_prime <- calculate_log_likelihood(Z, h_Z, observed_orders_idx, choice_sets,
                                            item_to_index, prob_noise_prime, mallow_theta, noise_option)
        
        log_acceptance_ratio <- llk_prime- llk_current
        acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
        if (runif(1) < acceptance_probability) {
          prob_noise <- prob_noise_prime
          num_acceptances <- num_acceptances + 1
          llk_current <- llk_prime
        }
      }
      
      # Store log-likelihood every iteration

      log_likelihood_currents <- c(log_likelihood_currents, llk_current)
    }
    
    # Update U (latent matrix Z) - PROPER IMPLEMENTATION
    else if (r <= (rho_pct + noise_pct + U_pct)) {
      # Update single element like Python
      i <- sample(1:n, 1)
      
      Z_prime <- Z
      Sigma_row <- build_sigma_rho(K, rho)
      Z_prime[i,] <- mvrnorm(1, Z[i,], Sigma_row)
      
      eta_prime <- transform_U_to_eta(Z_prime, alpha)
      h_Z_prime <- generate_partial_order(eta_prime)
      
      log_prior_current <- log_U_prior(Z, rho, K)
      log_prior_proposed <- log_U_prior(Z_prime, rho, K)
      

      llk_prime <- calculate_log_likelihood(Z_prime, h_Z_prime, observed_orders_idx, choice_sets,
                                          item_to_index, prob_noise, mallow_theta, noise_option)
      
      log_acceptance_ratio <- (log_prior_proposed + llk_prime) - (log_prior_current + llk_current)
      acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
      if (runif(1) < acceptance_probability) {
        Z <- Z_prime
        h_Z <- h_Z_prime
        llk_current <- llk_prime
        num_acceptances <- num_acceptances + 1
      }
      
      # Store log-likelihood every iteration
      log_likelihood_currents <- c(log_likelihood_currents, llk_current)
    }
    
    # Update beta
    else {
      j <- ((iteration-1) %% p) + 1  # Random component like Python
      epsilon <- rnorm(1, 0, drbeta * sigma_beta) 
      beta_prime <- beta
      beta_prime[j] <- beta_prime[j] + epsilon
      alpha_prime <- as.vector(t(X) %*% beta_prime)
      eta_prime <- transform_U_to_eta(Z, alpha_prime)
      h_Z_prime <- generate_partial_order(eta_prime)
      
      lp_current <- log_beta_prior(beta, sigma_beta)
      lp_proposed <- log_beta_prior(beta_prime, sigma_beta)

      llk_prime <- calculate_log_likelihood(Z, h_Z_prime, observed_orders_idx, choice_sets,
                                          item_to_index, prob_noise, mallow_theta, noise_option)
      
      log_acceptance_ratio <- (lp_proposed + llk_prime) - (lp_current + llk_current)
      acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
      if (runif(1) < acceptance_probability) {
        beta <- beta_prime
        llk_current <- llk_prime
        alpha <- alpha_prime
        h_Z <- h_Z_prime
        num_acceptances <- num_acceptances + 1
      }
      
      # Store log-likelihood every iteration
      log_likelihood_currents <- c(log_likelihood_currents, llk_current)
    }
    
    # Store current state every 100 iterations
    if (iteration %% 100 == 0) {
      Z_trace[[length(Z_trace) + 1]] <- Z
      h_trace[[length(h_trace) + 1]] <- h_Z
      beta_trace[[length(beta_trace) + 1]] <- beta
      rho_trace <- c(rho_trace, rho)
      prob_noise_trace <- c(prob_noise_trace, prob_noise)
      mallow_theta_trace <- c(mallow_theta_trace, mallow_theta)
    }
    
    current_acceptance_rate <- num_acceptances / iteration
    acceptance_rates <- c(acceptance_rates, current_acceptance_rate)
    
    if (iteration %in% progress_intervals) {
      cat(sprintf("Iteration %d/%d - Accept Rate: %.2f%% - Edges: %d\n", 
                  iteration, num_iterations, current_acceptance_rate * 100, sum(h_Z)))
    }
  }
  
  overall_acceptance_rate <- num_acceptances / num_iterations
  cat(sprintf("\nOverall Acceptance Rate after %d iterations: %.2f%%\n", 
              num_iterations, overall_acceptance_rate * 100))
  cat(sprintf("Final h_Z has %d edges\n", sum(h_Z)))
  
  return(list(
    Z_trace = Z_trace,
    h_trace = h_trace,
    index_to_item = index_to_item,
    item_to_index = item_to_index,
    rho_trace = rho_trace,
    beta_trace = beta_trace,
    prob_noise_trace = prob_noise_trace,
    mallow_theta_trace = mallow_theta_trace,
    log_likelihood_currents = log_likelihood_currents,
    acceptance_rates = acceptance_rates,
    overall_acceptance_rate = overall_acceptance_rate,
    final_h = h_Z,
    final_Z = Z,
    final_beta = beta,
    final_rho = rho,
    final_prob_noise = prob_noise,
    final_mallow_theta = mallow_theta
  ))
} 