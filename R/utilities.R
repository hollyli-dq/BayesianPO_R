#' Bayesian Partial Order Inference Package
#' 
#' This package provides functions for Bayesian inference of strong partial orders
#' from noisy observations using Markov Chain Monte Carlo (MCMC) methods.
#' 
#' @description
#' A strong partial order is a binary relation that satisfies:
#' - Irreflexivity: ¬(a ≺ a)
#' - Antisymmetry: If a ≺ b then ¬(b ≺ a)
#' - Transitivity: If a ≺ b and b ≺ c then a ≺ c
#' 
#' The model uses a latent space representation where each item j has a 
#' K-dimensional latent position U_j ∈ R^K, with correlation controlled by ρ.
 

# ============================================================================
# BASIC UTILITIES
# ============================================================================

#' Build correlation matrix with off-diagonal correlation rho
build_sigma_rho <- function(K, rho) {
  Sigma <- matrix(rho, nrow = K, ncol = K)
  diag(Sigma) <- 1.0
  return(Sigma)
}

#' Generate partial order matrix from latent positions
generate_partial_order <- function(eta) {
  n <- nrow(eta)
  h <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        # Item i precedes item j if all dimensions of eta_i are greater than eta_j
        h[i, j] <- as.integer(all(eta[i, ] > eta[j, ]))
      }
    }
  }
  return(h)
}

#' Compute transitive closure of a binary matrix
transitive_closure <- function(h) {
  n <- nrow(h)
  tc <- h
  
  # Replace any NAs with 0s for computation
  tc[is.na(tc)] <- 0
  
  for (k in 1:n) {
    for (i in 1:n) {
      for (j in 1:n) {
        tc[i, j] <- max(tc[i, j], tc[i, k] * tc[k, j])
      }
    }
  }
  return(tc)
}

#' Compute transitive reduction of a binary matrix
transitive_reduction <- function(h) {
  n <- nrow(h)
  tr <- h
  
  for (k in 1:n) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (tr[i, k] == 1 && tr[k, j] == 1) {
          tr[i, j] <- 0
        }
      }
    }
  }
  return(tr)
}

#' Count number of linear extensions (exact algorithm)
count_linear_extensions <- function(h) {
  # Convert to transitive reduction first
  tr <- transitive_reduction(h)
  return(nle(tr))
}

#' Count number of linear extensions of transitive reduction (exact algorithm)
nle <- function(tr) {
  # Base cases
  if (length(tr) == 0 || nrow(tr) <= 1) {
    return(1)
  }
  
  n <- nrow(tr)
  
  # Calculate column sums (incoming edges) and row sums (outgoing edges)
  cs <- colSums(tr)
  csi <- (cs == 0)  # nodes with no incoming edges (potential tops)
  bs <- rowSums(tr)
  bsi <- (bs == 0)  # nodes with no outgoing edges (potential bottoms)
  
  # Find free nodes (isolated nodes with no incoming or outgoing edges)
  free <- which(bsi & csi)
  k <- length(free)
  
  # If all nodes are free, return n!
  if (k == n) {
    return(factorial(n))
  }
  
  # Handle free nodes
  if (k > 0) {
    # Remove free nodes and calculate factorial contribution
    tr <- tr[-free, -free, drop = FALSE]
    fac <- factorial(n) / factorial(n - k)
  } else {
    fac <- 1
  }
  
  # Recompute after removing free nodes
  if (nrow(tr) == 0) {
    return(fac)
  }
  
  cs <- colSums(tr)
  csi <- (cs == 0)
  bs <- rowSums(tr)
  bsi <- (bs == 0)
  
  tops <- which(csi)  # nodes with no incoming edges
  bots <- which(bsi)  # nodes with no outgoing edges
  
  # Special case: if only 2 nodes remain, return fac
  if (nrow(tr) == 2) {
    return(fac)
  }
  
  # Special case: unique top and bottom
  if (length(tops) == 1 && length(bots) == 1) {
    i <- tops[1]
    j <- bots[1]
    
    # Check bounds
    if (i <= nrow(tr) && j <= ncol(tr)) {
      # Remove both top and bottom nodes
      remaining_indices <- setdiff(1:nrow(tr), c(i, j))
      if (length(remaining_indices) > 0) {
        trr <- tr[remaining_indices, remaining_indices, drop = FALSE]
        return(fac * nle(trr))
      } else {
        return(fac)
      }
    } else {
      return(0)
    }
  }
  
  # General case: iterate over all top elements
  count <- 0
  for (i in tops) {
    if (i <= nrow(tr)) {
      # Remove node i
      remaining_indices <- setdiff(1:nrow(tr), i)
      if (length(remaining_indices) > 0) {
        trr <- tr[remaining_indices, remaining_indices, drop = FALSE]
        count <- count + nle(trr)
      } else {
        count <- count + 1
      }
    }
  }
  
  return(fac * count)
}

#' Find top elements (nodes with no incoming edges)
find_tops <- function(tr) {
  incoming <- colSums(tr)
  tops <- which(incoming == 0)
  return(tops)
}

#' Count linear extensions starting with specific element
num_extensions_with_first <- function(tr, first_item_idx) {
  # Identify top elements of the current poset
  tops <- find_tops(tr)
  
  # If first_item_idx is not a top element, no linear extension can start with it
  if (!(first_item_idx %in% tops)) {
    return(0)
  }
  
  # If it is top, remove it from tr and count the nle of the reduced poset
  remaining_indices <- setdiff(1:nrow(tr), first_item_idx)
  if (length(remaining_indices) > 0) {
    tr_reduced <- tr[remaining_indices, remaining_indices, drop = FALSE]
    return(nle(tr_reduced))
  } else {
    return(1)
  }
}

#' Restrict partial order to subset
restrict_partial_order <- function(h, subset) {
  return(h[subset, subset, drop = FALSE])
}

#' Log prior for latent variables U
log_U_prior <- function(Z, rho, K) {
  Sigma <- build_sigma_rho(K, rho)
  n <- nrow(Z)
  log_prior <- 0
  
  for (i in 1:n) {
    log_prior <- log_prior + dmvnorm(Z[i, ], mean = rep(0, K), sigma = Sigma, log = TRUE)
  }
  return(log_prior)
}

#' Sample conditional column for reversible jump MCMC
sample_conditional_column <- function(Z, rho) {
  n <- nrow(Z)
  K <- ncol(Z)
  
  if (K == 0) {
    return(rnorm(n, 0, 1))
  }
  
  # Build the (K+1)x(K+1) covariance matrix
  Kplus1 <- K + 1
  Sigma_full <- build_sigma_rho(Kplus1, rho)
  
  # Partition the covariance matrix
  Sigma_gg <- Sigma_full[1:K, 1:K, drop = FALSE]
  Sigma_dg <- Sigma_full[Kplus1, 1:K, drop = FALSE]
  Sigma_dd <- Sigma_full[Kplus1, Kplus1]
  
  # Invert Sigma_gg once for all rows
  Sigma_gg_inv <- solve(Sigma_gg)
  
  bridging_col <- numeric(n)
  for (i in 1:n) {
    x_i <- Z[i, ]  # existing coordinates
    # Conditional mean
    mu_cond <- as.numeric(Sigma_dg %*% Sigma_gg_inv %*% x_i)
    # Conditional variance
    var_cond <- Sigma_dd - as.numeric(Sigma_dg %*% Sigma_gg_inv %*% t(Sigma_dg))
    # Sample
    bridging_col[i] <- rnorm(1, mu_cond, sqrt(var_cond))
  }
  
  return(bridging_col)
}

#' Sample conditional element for latent matrix update
sample_conditional_z <- function(Z, rZ, cZ, rho) {
  K <- ncol(Z)
  
  # Special case for K=1: just sample from standard normal
  if (K == 1) {
    return(rnorm(1, 0, 1))
  }
  
  # Build correlation matrix
  Sigma <- build_sigma_rho(K, rho)
  
  dependent_ind <- cZ
  given_inds <- setdiff(1:K, cZ)
  
  # Handle case where there are no given indices (shouldn't happen with K>1, but just in case)
  if (length(given_inds) == 0) {
    return(rnorm(1, 0, 1))
  }
  
  Sigma_dd <- Sigma[dependent_ind, dependent_ind]  # scalar
  Sigma_dg <- Sigma[dependent_ind, given_inds, drop = FALSE]  # row vector
  Sigma_gg <- Sigma[given_inds, given_inds, drop = FALSE]
  
  # X_g is the given values - ensure it's a column vector
  X_g <- as.matrix(Z[rZ, given_inds, drop = FALSE])
  if (nrow(X_g) == 1) {
    X_g <- t(X_g)  # Make it a column vector
  }
  
  # Means are 0
  mu_d <- 0.0
  mu_g <- matrix(0, nrow = length(given_inds), ncol = 1)
  
  # Invert Sigma_gg
  if (length(given_inds) == 1) {
    # Special case for single given variable
    Sigma_gg_inv <- 1 / (Sigma_gg + 1e-8)
    mu_cond <- mu_d + as.numeric(Sigma_dg * Sigma_gg_inv * (X_g - mu_g))
    var_cond <- Sigma_dd - as.numeric(Sigma_dg * Sigma_gg_inv * Sigma_dg)
  } else {
    # General case for multiple given variables
    Sigma_gg_inv <- solve(Sigma_gg + diag(1e-8, nrow(Sigma_gg)))
    mu_cond <- mu_d + as.numeric(Sigma_dg %*% Sigma_gg_inv %*% (X_g - mu_g))
    var_cond <- Sigma_dd - as.numeric(Sigma_dg %*% Sigma_gg_inv %*% t(Sigma_dg))
  }
  
  var_cond <- max(var_cond, 1e-8)
  Z_new <- rnorm(1, mu_cond, sqrt(var_cond))
  
  return(Z_new)
}

# ============================================================================
# PRIOR DISTRIBUTIONS
# ============================================================================

#' Log prior for rho parameter (Beta distribution with truncation)
log_rho_prior <- function(rho, fac = 1/6, tol = 1e-4) {
  if (rho > 1 - tol) {
    return(-Inf)
  }
  log_pdf <- dbeta(rho, 1, fac, log = TRUE)
  log_cdf_trunc <- pbeta(1 - tol, 1, fac, log.p = TRUE)
  return(log_pdf - log_cdf_trunc)
}

#' Sample from rho prior (Beta distribution with truncation)
sample_rho_prior <- function(fac = 1/6, tol = 1e-4) {
  repeat {
    rho <- rbeta(1, 1, fac)
    if (1 - rho >= tol) {
      return(rho)
    }
  }
}

#' Log prior for noise probability (Beta distribution)
log_noise_prior <- function(prob_noise, beta_param) {
  if (prob_noise <= 0 || prob_noise >= 1) {
    return(-Inf)
  }
  return(dbeta(prob_noise, 1, beta_param, log = TRUE))
}

#' Sample from noise probability prior
sample_noise_prior <- function(beta_param) {
  return(rbeta(1, 1, beta_param))
}

#' Log prior for beta coefficients (multivariate normal)
log_beta_prior <- function(beta, sigma_beta) {
  p <- length(beta)
  
  if (length(sigma_beta) == 1) {
    # Scalar case
    log_det_part <- -0.5 * p * log(2 * pi) - p * log(sigma_beta)
    quad_part <- -0.5 * sum(beta^2) / (sigma_beta^2)
  } else {
    # Vector case
    if (length(sigma_beta) != p) {
      stop("sigma_beta must be scalar or have same length as beta")
    }
    log_det_part <- -0.5 * p * log(2 * pi) - sum(log(sigma_beta))
    quad_part <- -0.5 * sum(beta^2 / (sigma_beta^2))
  }
  
  return(log_det_part + quad_part)
}

#' Sample from beta prior
sample_beta_prior <- function(sigma_beta, p) {
  if (length(sigma_beta) == 1) {
    return(rnorm(p, 0, sigma_beta))
  } else {
    if (length(sigma_beta) != p) {
      stop("sigma_beta must be scalar or have length p")
    }
    return(rnorm(p, 0, sigma_beta))
  }
}

#' Log prior for K (truncated Poisson)
log_K_prior <- function(K, K_prior) {
  if (K < 1) {
    return(-Inf)
  }
  # log(k!) using lgamma(k+1)
  log_k_fact <- lgamma(K + 1)
  # normalizing constant for truncation
  norm_const <- -log(1 - exp(-K_prior))
  val <- -K_prior + K * log(K_prior) - log_k_fact + norm_const
  return(val)
}

#' Sample from K prior (truncated Poisson)
sample_K_prior <- function(K_prior) {
  repeat {
    candidate <- rpois(1, K_prior)
    if (candidate >= 1) {
      return(candidate)
    }
  }
}

# ============================================================================
# Convert matrix to data list
# ============================================================================
# nrow(X) is a number of lists
# 0 indicates missing.
mx2list <- function(X){
  observed_orders=list()
  for(i in 1:nrow(X)){
    observed_orders[[i]]=paste0("Item_",X[i,which(X[i,]>0)])
  }
  choice_sets = observed_orders 
  items = paste0("Item_",c(1:max(X)))
  return(list(observed_orders = observed_orders, 
              choice_sets = choice_sets, items=items))
}


# ============================================================================
# DATA GENERATION FUNCTIONS
# ============================================================================

#' Generate synthetic data for testing
generate_synthetic_data <- function(n_items = 6, n_observations = 50, min_sub =2 ,K = 3, 
                                   rho_true = 0.8, prob_noise_true = 0.1, 
                                   random_seed = 123) {
  set.seed(random_seed)
  
  # Generate true latent variables
  Sigma_true <- build_sigma_rho(K, rho_true)
  Z_true <- rmvnorm(n_items, mean = rep(0, K), sigma = Sigma_true)
  
  # Transform to eta and generate true partial order
  eta_true <- Z_true
  
  # First generate the full partial order
  h_full <- generate_partial_order(eta_true)
  
  # Then compute its transitive reduction to remove redundant edges
  h_true <- transitive_closure(h_full)
  
  # Generate item names
  items <- paste0("Item_", 1:n_items)

  subsets <- list()
  if(min_sub!=n){
    for (i in 1:N) {
      # Random choice set (subset of items)
      choice_size <- sample(min_sub:n, 1)
      choice_set <- sample(items, choice_size)
      subsets[[i]] <- choice_set
    }
  } else {subsets = rep(list(items),N)} 
  
  # Generate observations with noise
  observed_orders <- list()
  choice_sets <- list()
  
  for (i in 1:n_observations) {
    choice_set <- subsets[[i]]
    choice_sets[[i]] <-choice_set
    
    if (noise_option == "queue_jump") {
      order <- generate_total_order_queue_jump(
          subset     = choice_set,
          items_all  = items,
          h_global   = h_true,
          prob_noise = prob_noise_true
      
      )
    }
    observed_orders[[i]] <- order   # one assignment is enough
  }                                # closes the *for* loop
 
  return(list(
    observed_orders = observed_orders,
    choice_sets = choice_sets,
    items = items,
    h_true = h_true,
    Z_true = Z_true,
    rho_true = rho_true,
    prob_noise_true = prob_noise_true
  ))
} 

#' Queue-jump total-order generator (one choice-set)
generate_total_order_queue_jump <- function(subset,
                                            items_all,
                                            h_global,
                                            prob_noise = 0.1) {
  subset <- sort(unique(subset))
  if (length(subset) == 0) return(integer(0))
  
  ## 1.  Extract the sub-matrix for this subset 
  idx_map <- match(subset, items_all)                # global → row/col index
  h_sub   <- h_global[idx_map, idx_map, drop = FALSE]
  
  ## 2.  Work with LOCAL indices 
  remaining   <- seq_along(subset)                   # local indices
  order_local <- integer(0)
  
  while (length(remaining) > 0) {
    m <- length(remaining)
    if (m == 1) {                      # only one left – append and stop
      order_local <- c(order_local, remaining)
      break
    }
    
    # (a) current residual partial order
    h_rem <- h_sub[remaining, remaining, drop = FALSE]
    tr_rem <- transitive_reduction(h_rem)
    
    # (b) total number of linear extensions
    N_total <- nle(tr_rem)
    if (N_total == 0)
      stop("Residual graph is cyclic – cannot draw a linear extension")
    
    # (c) probability of each candidate being first, without jump
    p_no_jump <- numeric(m)
    for (k in seq_len(m)) {
      N_first      <- num_extensions_with_first(tr_rem, k)
      p_no_jump[k] <- (1 - prob_noise) * (N_first / N_total)
    }
    
    # (d) add the jump component (uniform across m candidates)
    p_total <- p_no_jump + prob_noise / m
    p_total <- p_total / sum(p_total)                # normalise
    
    # (e) sample the next element
    chosen_pos   <- sample(seq_len(m), size = 1, prob = p_total)
    chosen_local <- remaining[chosen_pos]
    
    order_local <- c(order_local, chosen_local)
    remaining   <- remaining[-chosen_pos]
  }
  
  ## 3.  Convert local indices back to *global* 
  subset[order_local]
}
