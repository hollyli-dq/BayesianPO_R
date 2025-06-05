# Bayesian Partial Order Inference - Reversible Jump MCMC
# Variable dimension MCMC for automatic model selection (CORRECTED TO MATCH PYTHON)
# Author: Converted for Prof Geoff Nicholls, University of Oxford

# Required libraries
library(mvtnorm)

#' Reversible Jump MCMC for Partial Order Inference (CORRECTED IMPLEMENTATION)
#' 
#' Perform MCMC sampling to infer the partial order h, plus parameters (rho, prob_noise, mallow_theta).
#' This function includes reversible jump moves to automatically determine the optimal number of latent dimensions K.
#' 
#' @param observed_orders List of observed total orders
#' @param choice_sets List of choice sets for each observation
#' @param num_iterations Number of MCMC iterations
#' @param X Design matrix for covariates
#' @param dr Multiplicative step size for rho
#' @param drbeta Multiplicative step size for beta
#' @param sigma_mallow Parameter controlling random-walk for mallow_theta
#' @param sigma_beta Standard deviation for beta prior
#' @param noise_option Noise model ("queue_jump" or "mallows_noise")
#' @param mcmc_pt Vector of update probabilities (rho, noise, U, beta, K)
#' @param rho_prior Prior parameter for rho
#' @param noise_beta_prior Prior parameter for noise
#' @param mallow_ua Prior parameter for Mallows model
#' @param K_prior Prior parameter for K (Poisson rate)
#' @param random_seed Random seed for reproducibility
#' @return List containing MCMC traces and diagnostics
#' @export
mcmc_partial_order_k <- function(
    observed_orders,
    choice_sets,
    num_iterations,
    X,
    dr = 0.1,
    drbeta = 0.1,
    sigma_mallow = 0.1,
    sigma_beta = 1.0,
    noise_option = "queue_jump",
    mcmc_pt = c(0.15, 0.15, 0.3, 0.15, 0.25),
    rho_prior = 1/6,
    noise_beta_prior = 9,
    mallow_ua = 1,
    K_prior = 3.0,
    random_seed = 123
) {

    
    noise_option <- match.arg(noise_option)
    
    ## -------------------------------------------------------------------------
    ## Enforce that only the relevant prior is supplied for the chosen noise
    ## model.  (Same logic as in mcmc_partial_order.)
    ## -------------------------------------------------------------------------
    if (noise_option == "queue_jump") {
        if (is.null(noise_beta_prior))
            stop("`noise_beta_prior` must be provided when noise_option = 'queue_jump'")
        if (!is.null(mallow_ua))
            message("Ignoring `mallow_ua` because queue-jump noise is selected.")
    } else {  # mallows_noise
        if (is.null(mallow_ua))
            stop("`mallow_ua` must be provided when noise_option = 'mallows_noise'")
        if (!is.null(noise_beta_prior))
            message("Ignoring `noise_beta_prior` because Mallows noise is selected.")
    }
    
    # Set random seed
    set.seed(random_seed)
    
    # ----------------------------------------------------------------
    # 1. Setup: Map items to indices, initialize states, etc.
    # ----------------------------------------------------------------
    
    items <- sort(unique(unlist(choice_sets)))
    n <- length(items)
    item_to_index <- setNames(0:(length(items)-1), items)
    index_to_item <- setNames(items, 0:(length(items)-1))
    
    # Convert observed orders to index form
    observed_orders_idx <- lapply(observed_orders, function(order) {
        sapply(order, function(item) item_to_index[[as.character(item)]])
    })
    
    # Initialize MCMC state (MATCH PYTHON EXACTLY)
    K <- 1  # Start with K=1 like Python
    Z <- matrix(0, nrow = n, ncol = K)  # Start with zeros like Python
    
    p <- nrow(X)
    beta <- rnorm(p, mean = 0, sd = sigma_beta)
    alpha <- as.vector(t(X) %*% beta)
    eta <- transform_U_to_eta(Z, alpha)
    h_Z <- generate_partial_order(eta)
    
    # Initialize parameters using proper priors (MATCH PYTHON)
    rho <- sample_rho_prior(rho_prior)
    prob_noise <- sample_noise_prior(noise_beta_prior)
    mallow_theta <- sample_theta_prior(mallow_ua)
    
    # ----------------------------------------------------------------
    # 2. Prepare Storage for MCMC results
    # ----------------------------------------------------------------
    
    Z_trace <- list()
    h_trace <- list()
    K_trace <- numeric()
    beta_trace <- list()
    update_records <- list()
    
    rho_trace <- numeric()
    prob_noise_trace <- numeric()
    mallow_theta_trace <- numeric()
    
    proposed_rho_vals <- numeric()
    proposed_prob_noise_vals <- numeric()
    proposed_mallow_theta_vals <- numeric()
    proposed_beta_vals <- list()
    proposed_Zs <- list()
    acceptance_decisions <- numeric()
    acceptance_rates <- numeric()
    log_likelihood_currents <- numeric()
    log_likelihood_primes <- numeric()
    
    num_acceptances <- 0
    
    # Progress intervals (10% increments)
    progress_intervals <- round(seq(0.1, 1.0, by = 0.1) * num_iterations)
    
    # Unpack update probabilities
    rho_pct <- mcmc_pt[1]
    noise_pct <- mcmc_pt[2]
    U_pct <- mcmc_pt[3]
    beta_pct <- mcmc_pt[4]
    K_pct <- mcmc_pt[5]
    
    llk_current <- -Inf
    llk_prime <- -Inf
    
    cat("Starting RJ-MCMC with", num_iterations, "iterations...\n")
    cat("Initial K =", K, ", h_Z has", sum(h_Z), "edges\n")
    
    # ----------------------------------------------------------------
    # 3. Main MCMC Loop (MATCH PYTHON EXACTLY)
    # ----------------------------------------------------------------
    
    for (iteration in 1:num_iterations) {
        r <- runif(1)
        accepted_this_iter <- FALSE
        update_category <- NULL
        
        # Calculate current likelihood (LIKE PYTHON)
        llk_current <- calculate_log_likelihood(
            Z, h_Z, observed_orders_idx, choice_sets,
            item_to_index, prob_noise, mallow_theta, noise_option
        )
        
        # ---- A) Update rho (MATCH PYTHON EXACTLY) ----
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
            else if (r < (rho_pct + noise_pct)) {
            update_category <- "noise"
            
            if (noise_option == "mallows_noise") {
                epsilon <- rnorm(1, 0, 1)
                mallow_theta_prime <- mallow_theta * exp(sigma_mallow * epsilon)
                
                log_prior_current <- log_theta_prior(mallow_theta, mallow_ua)
                log_prior_proposed <- log_theta_prior(mallow_theta_prime, mallow_ua)
                
                llk_prime <- calculate_log_likelihood(
                    Z, h_Z, observed_orders_idx, choice_sets, item_to_index,
                    prob_noise, mallow_theta_prime, noise_option
                )
                
                log_acceptance_ratio <- (log_prior_proposed + llk_prime) - 
                                     (log_prior_current + llk_current) + 
                                     log(mallow_theta / mallow_theta_prime)  # MATCH PYTHON
                
                acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
                if (runif(1) < acceptance_probability) {
                    mallow_theta <- mallow_theta_prime
                    num_acceptances <- num_acceptances + 1
                    acceptance_decisions <- c(acceptance_decisions, 1)
                    accepted_this_iter <- TRUE
                    llk_current <- llk_prime
                } else {
                    acceptance_decisions <- c(acceptance_decisions, 0)
                }
                proposed_mallow_theta_vals <- c(proposed_mallow_theta_vals, mallow_theta_prime)
                
            } else if (noise_option == "queue_jump") {
                prob_noise_prime <- sample_noise_prior(noise_beta_prior)
                
                log_prior_current <- log_noise_prior(prob_noise, noise_beta_prior)
                log_prior_proposed <- log_noise_prior(prob_noise_prime, noise_beta_prior)
                
                llk_prime <- calculate_log_likelihood(
                    Z, h_Z, observed_orders_idx, choice_sets, item_to_index,
                    prob_noise_prime, mallow_theta, noise_option
                )
                
                # MATCH PYTHON: only likelihood ratio for queue jump
                log_acceptance_ratio <- llk_prime - llk_current
                acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
                
                if (runif(1) < acceptance_probability) {
                    prob_noise <- prob_noise_prime
                    num_acceptances <- num_acceptances + 1
                    accepted_this_iter <- TRUE
                    acceptance_decisions <- c(acceptance_decisions, 1)
                    llk_current <- llk_prime
                } else {
                    acceptance_decisions <- c(acceptance_decisions, 0)
                }
                proposed_prob_noise_vals <- c(proposed_prob_noise_vals, prob_noise_prime)
            }
            
        # ---- C) Update U (latent matrix Z) - MATCH PYTHON EXACTLY ----
        } else if (r <= (rho_pct + noise_pct + U_pct)) {
            update_category <- "U"
            i <- sample(0:(n-1), 1)  # 0-indexed like Python
            current_row <- Z[i+1, ]  # Convert to 1-indexed for R
            
            # Build proposal covariance matrix (MATCH PYTHON)
            Sigma <- build_sigma_rho(K, rho)
            
            proposed_row <- rmvnorm(1, current_row, Sigma)[1, ]
            Z_prime <- Z
            Z_prime[i+1, ] <- proposed_row  # Convert to 1-indexed for R
            
            eta_prime <- transform_U_to_eta(Z_prime, alpha)
            h_Z_prime <- generate_partial_order(eta_prime)
            
            log_prior_current <- log_U_prior(Z, rho, K)
            log_prior_proposed <- log_U_prior(Z_prime, rho, K)
            
            llk_prime <- calculate_log_likelihood(
                Z_prime, h_Z_prime, observed_orders_idx, choice_sets, item_to_index,
                prob_noise, mallow_theta, noise_option
            )
            
            log_acceptance_ratio <- (log_prior_proposed + llk_prime) - (log_prior_current + llk_current)
            acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
            
            if (runif(1) < acceptance_probability) {
                Z <- Z_prime
                h_Z <- h_Z_prime
                num_acceptances <- num_acceptances + 1
                acceptance_decisions <- c(acceptance_decisions, 1)
                accepted_this_iter <- TRUE
                llk_current <- llk_prime
            } else {
                acceptance_decisions <- c(acceptance_decisions, 0)
            }
            proposed_Zs <- append(proposed_Zs, list(Z_prime))
            
        # ---- D) Update beta (MATCH PYTHON EXACTLY) ----
        } else if (r <= (rho_pct + noise_pct + U_pct + beta_pct)) {
            update_category <- "beta"
            j <- ((iteration - 1) %% p) + 1  # MATCH PYTHON: cycle through components
            epsilon <- rnorm(1, 0, drbeta * sigma_beta)  # MATCH PYTHON scaling
            beta_prime <- beta
            beta_prime[j] <- beta_prime[j] + epsilon
            alpha_prime <- as.vector(t(X) %*% beta_prime)
            eta_prime <- transform_U_to_eta(Z, alpha_prime)
            h_Z_prime <- generate_partial_order(eta_prime)
            
            lp_current <- log_beta_prior(beta, sigma_beta)
            lp_proposed <- log_beta_prior(beta_prime, sigma_beta)
            
            llk_prime <- calculate_log_likelihood(
                Z, h_Z_prime, observed_orders_idx, choice_sets, item_to_index,
                prob_noise, mallow_theta, noise_option
            )
            
            log_acceptance_ratio <- (lp_proposed + llk_prime) - (lp_current + llk_current)
            acceptance_probability <- min(1.0, exp(log_acceptance_ratio))
            
            if (runif(1) < acceptance_probability) {
                beta <- beta_prime
                alpha <- alpha_prime
                h_Z <- h_Z_prime
                num_acceptances <- num_acceptances + 1
                acceptance_decisions <- c(acceptance_decisions, 1)
                accepted_this_iter <- TRUE
                llk_current <- llk_prime
            } else {
                acceptance_decisions <- c(acceptance_decisions, 0)
            }
            proposed_beta_vals <- append(proposed_beta_vals, list(beta_prime))
            
        # ---- E) Update K (Reversible Jump) - MATCH PYTHON EXACTLY ----
        } else {
            update_category <- "K"
            
            if (K == 1) {
                move <- "up"
            } else {
                move <- ifelse(runif(1) < 0.5, "up", "down")
            }
            
            if (move == "up") {
                # Birth move: K -> K+1
                K_prime <- K + 1
                col_ins <- sample(0:K, 1)  # MATCH PYTHON: 0 to K positions
                
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
                
                eta_prime <- transform_U_to_eta(Z_prime, alpha)
                h_Z_prime <- generate_partial_order(eta_prime)
                
                log_prior_K <- log_K_prior(K, K_prior)
                log_prior_K_prime <- log_K_prior(K_prime, K_prior)
                
                
                llk_prime <- calculate_log_likelihood(
                    Z_prime, h_Z_prime, observed_orders_idx, choice_sets,
                    item_to_index, prob_noise, mallow_theta, noise_option
                )
                
                # Jacobian terms (MATCH PYTHON EXACTLY)
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
                    num_acceptances <- num_acceptances + 1
                    acceptance_decisions <- c(acceptance_decisions, 1)
                    accepted_this_iter <- TRUE
                    llk_current <- llk_prime
                } else {
                    acceptance_decisions <- c(acceptance_decisions, 0)
                }
                
            } else {
                # Death move: K -> K-1 (MATCH PYTHON EXACTLY)
                K_prime <- K - 1
                col_del <- sample(0:(K-1), 1)  # MATCH PYTHON: 0 to K-1
                Z_prime <- Z[, -(col_del+1), drop = FALSE]  # Convert to 1-indexed for R
                eta_prime <- transform_U_to_eta(Z_prime, alpha)
                h_Z_prime <- generate_partial_order(eta_prime)
                
                log_prior_K <- log_K_prior(K, K_prior)
                log_prior_K_prime <- log_K_prior(K_prime, K_prior)
                
                llk_prime <- calculate_log_likelihood(
                    Z_prime, h_Z_prime, observed_orders_idx, choice_sets,
                    item_to_index, prob_noise, mallow_theta, noise_option
                )
                
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
                    num_acceptances <- num_acceptances + 1
                    acceptance_decisions <- c(acceptance_decisions, 1)
                    accepted_this_iter <- TRUE
                    llk_current <- llk_prime
                } else {
                    acceptance_decisions <- c(acceptance_decisions, 0)
                }
            }
        }
        
        # Store current state every 100 iterations (MATCH PYTHON)
        if (iteration %% 100 == 0) {
            Z_trace <- append(Z_trace, list(Z))
            h_trace <- append(h_trace, list(h_Z))
            K_trace <- c(K_trace, K)
            beta_trace <- append(beta_trace, list(beta))
            rho_trace <- c(rho_trace, rho)
            prob_noise_trace <- c(prob_noise_trace, prob_noise)
            mallow_theta_trace <- c(mallow_theta_trace, mallow_theta)
            update_records <- append(update_records, list(c(iteration, update_category, accepted_this_iter)))
        }
        
        log_likelihood_currents <- c(log_likelihood_currents, llk_current)
        log_likelihood_primes <- c(log_likelihood_primes, llk_prime)
        current_acceptance_rate <- num_acceptances / iteration
        acceptance_rates <- c(acceptance_rates, current_acceptance_rate)
        
        if (iteration %in% progress_intervals) {
            cat(sprintf("Iteration %d/%d - Accept Rate: %.2f%% - K: %d - Edges: %d\n", 
                       iteration, num_iterations, current_acceptance_rate * 100, K, sum(h_Z)))
        }
    }
    
    overall_acceptance_rate <- num_acceptances / num_iterations
    cat(sprintf("\nOverall Acceptance Rate after %d iterations: %.2f%%\n", 
               num_iterations, overall_acceptance_rate * 100))
    cat(sprintf("Final K = %d, h_Z has %d edges\n", K, sum(h_Z)))
    
    # Create update dataframe
    update_df <- data.frame(
        iteration = sapply(update_records, function(x) x[1]),
        category = sapply(update_records, function(x) x[2]),
        accepted = as.logical(sapply(update_records, function(x) x[3]))
    )
    
    return(list(
        Z_trace = Z_trace,
        h_trace = h_trace,
        K_trace = K_trace,
        beta_trace = beta_trace,
        index_to_item = index_to_item,
        item_to_index = item_to_index,
        rho_trace = rho_trace,
        prob_noise_trace = prob_noise_trace,
        mallow_theta_trace = mallow_theta_trace,
        proposed_rho_vals = proposed_rho_vals,
        proposed_prob_noise_vals = proposed_prob_noise_vals,
        proposed_mallow_theta_vals = proposed_mallow_theta_vals,
        proposed_beta_vals = proposed_beta_vals,
        proposed_Zs = proposed_Zs,
        acceptance_rates = acceptance_rates,
        acceptance_decisions = acceptance_decisions,
        log_likelihood_currents = log_likelihood_currents,
        log_likelihood_primes = log_likelihood_primes,
        overall_acceptance_rate = overall_acceptance_rate,
        update_df = update_df,
        final_h = h_Z,
        final_Z = Z,
        final_K = K,
        final_beta = beta,
        final_rho = rho,
        final_prob_noise = prob_noise,
        final_mallow_theta = mallow_theta
    ))
} 