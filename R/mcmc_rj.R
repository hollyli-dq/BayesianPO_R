mcmc_partial_order_k <- function(
    observed_orders,
    choice_sets,
    num_iterations,
    dr = 0.1,
    rho_prior = 1/6,
    K_prior = 3.0,
    noise_option = "queue_jump",
    noise_beta_prior = NULL, 
    random_seed = 123
){
  
  # ----------------------------------------------------------------
  # 1. Setup: Map items to indices, initialize states, etc.
  # ----------------------------------------------------------------
  
  
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
  
  # Initialize MCMC state (MATCH PYTHON EXACTLY)
  K <- 1  # Start with K=1 like Python
  Z <- rmvnorm(n, mean = rep(0, K), sigma = diag(K))
  
  eta <- Z
  h_Z <- generate_partial_order(eta)
  
  # Initialize parameters using proper priors (MATCH PYTHON)
  rho <- sample_rho_prior(rho_prior)
  prob_noise <- sample_noise_prior(noise_beta_prior)
  llk_current <- log_likelihood_queue_jump(h_Z, observed_orders_idx, choice_sets, item_to_index, prob_noise)
  
  Z_trace <- list()
  h_trace <- list()
  K_trace <- numeric()
  rho_trace <- numeric()
  prob_noise_trace <- numeric()
  log_likelihood_trace <- numeric()
  acceptance_rates <- data.frame(list(rho=0,prob_noise=0,Z=0,K=0))
  
  # Progress intervals (10% increments)
  progress_intervals <- round(seq(0.1, 1.0, by = 0.1) * num_iterations)
  
  # ----------------------------------------------------------------
  # 3. Main MCMC Loop (MATCH PYTHON EXACTLY)
  # ----------------------------------------------------------------
  
  for (iteration in 1:num_iterations) {
    
    delta <- runif(1, dr, 1/dr)  # Correct range
    rho_prime <- 1.0 - (1.0 - rho) * delta
    
    if (!is.na(rho_prime) && rho_prime > 0 && rho_prime < 1) {
      log_prior_current <- log_rho_prior(rho, rho_prior) + log_U_prior(Z, rho, K)
      log_prior_proposed <- log_rho_prior(rho_prime, rho_prior) + log_U_prior(Z, rho_prime, K)
      log_acceptance_ratio <- (log_prior_proposed) - (log_prior_current)- log(delta)
      
      acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
      if (runif(1) < acceptance_probability) {
        rho <- rho_prime
        acceptance_rates$rho = acceptance_rates$rho+1
      }
    }
    
    ################## noise prob update ################
    
    prob_noise_prime <- sample_noise_prior(noise_beta_prior)
    
    log_prior_current <- log_noise_prior(prob_noise, noise_beta_prior)
    log_prior_proposed <- log_noise_prior(prob_noise_prime, noise_beta_prior)
    
    llk_prime <- log_likelihood_queue_jump(h_Z, observed_orders_idx, choice_sets, item_to_index, prob_noise_prime)
    
    # MATCH PYTHON: only likelihood ratio for queue jump
    log_acceptance_ratio <- llk_prime - llk_current
    acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
    
    if (runif(1) < acceptance_probability) {
      prob_noise <- prob_noise_prime
      llk_current <- llk_prime
      acceptance_rates$prob_noise  = acceptance_rates$prob_noise+1
    }
    
    ############# Z-update ################
    
    i <- sample(0:(n-1), 1)  # 0-indexed like Python
    current_row <- Z[i+1, ]  # Convert to 1-indexed for R
    
    # Build proposal covariance matrix (MATCH PYTHON)
    Sigma <- build_sigma_rho(K, rho)
    
    proposed_row <- rmvnorm(1, current_row, Sigma)[1, ]
    Z_prime <- Z
    Z_prime[i+1, ] <- proposed_row  # Convert to 1-indexed for R
    
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
      acceptance_rates$Z = acceptance_rates$Z+1
    } 
    
    #################### K-update ###################333
    
    if (K == 1) {
      move <- "up"
    } else {
      move <- ifelse(runif(1) < 0.5, "up", "down")
    }
    
    if (move == "up") {
      # Birth move: K -> K+1
      K_prime <- K + 1
      col_ins <- sample(0:K, 1)  
      
      b_col <- sample_conditional_column(Z, rho)
      
      # Insert column using np.insert logic (MATCH PYTHON)
      if (col_ins == 0) {
        Z_prime <- cbind(b_col, Z)
      } else if (col_ins == K) {
        Z_prime <- cbind(Z, b_col)
      } else {
        Z_prime <- cbind(Z[, 1:col_ins, drop = FALSE], 
                         b_col, 
                         Z[, (col_ins+1):K, drop = FALSE])
      }
      
      eta_prime <- Z_prime
      h_Z_prime <- generate_partial_order(eta_prime)
      
      log_prior_K <- log_K_prior(K, K_prior)
      log_prior_K_prime <- log_K_prior(K_prime, K_prior)
      
      llk_prime <- log_likelihood_queue_jump(h_Z_prime, observed_orders_idx, choice_sets, item_to_index, prob_noise)
      
      # Jacobian terms  
      rho_fk <- ifelse(K == 1, 1.0, 0.5)
      rho_bk <- 0.5
      
      log_acc <- (log_prior_K_prime + llk_prime) - 
        (log_prior_K + llk_current) + 
        log(rho_bk) - log(rho_fk)
      
      accept_prob <- min(1.0, exp(log_acc))
      if (runif(1) < accept_prob) {
        Z <- Z_prime
        K <- K_prime
        h_Z <- h_Z_prime
        llk_current <- llk_prime
        acceptance_rates$K = acceptance_rates$K+1
      } 
      
    } 
    if (move == "down") {
      # Death move: K -> K-1 (MATCH PYTHON EXACTLY)
      K_prime <- K - 1
      col_del <- sample(0:(K-1), 1)  # MATCH PYTHON: 0 to K-1
      Z_prime <- Z[, -(col_del+1), drop = FALSE]  # Convert to 1-indexed for R
      eta_prime <- Z_prime
      h_Z_prime <- generate_partial_order(eta_prime)
      
      log_prior_K <- log_K_prior(K, K_prior)
      log_prior_K_prime <- log_K_prior(K_prime, K_prior)
      
      llk_prime <- log_likelihood_queue_jump(h_Z_prime, observed_orders_idx, choice_sets, item_to_index, prob_noise)
      
      
      # Jacobian terms for death move
      rho_fk <- 0.5
      rho_bk <- ifelse(K_prime == 1, 1.0, 0.5)
      
      log_acc <- (log_prior_K_prime + llk_prime) - (log_prior_K + llk_current) + 
        log(rho_fk) - log(rho_bk)
      accept_prob <- min(1.0, exp(log_acc))
      
      if (runif(1) < accept_prob) {
        Z <- Z_prime
        K <- K_prime
        h_Z <- h_Z_prime
        llk_current <- llk_prime
        acceptance_rates$K = acceptance_rates$K+1
      } 
    }
    
    # Store current state every 10 iterations
    if (iteration %% 10 == 0) {
      Z_trace <- append(Z_trace, list(Z))
      h_trace <- append(h_trace, list(h_Z))
      K_trace <- c(K_trace, K)
      rho_trace <- c(rho_trace, rho)
      prob_noise_trace <- c(prob_noise_trace, prob_noise)
      log_likelihood_trace <- c(log_likelihood_trace,llk_current)
    }
    
    if (iteration %in% progress_intervals) {
      cat(sprintf("Iteration %d/%d - Accept Rate: %.2f%% - K: %d - Edges: %d\n", 
                  iteration, num_iterations, median(unlist(acceptance_rates/iteration * 100)), K, sum(h_Z)))
    }
  }
  
  acceptance_rates=acceptance_rates/num_iterations
  
  return(list(
    Z_trace = Z_trace,
    h_trace = h_trace,
    K_trace = K_trace,
    index_to_item = index_to_item,
    item_to_index = item_to_index,
    rho_trace = rho_trace,
    prob_noise_trace = prob_noise_trace,
    acceptance_rates = acceptance_rates,
    log_likelihood_trace = log_likelihood_trace,
    final_h = h_Z,
    final_Z = Z,
    final_K = K,
    final_rho = rho,
    final_prob_noise = prob_noise
  ))
}
