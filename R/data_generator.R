# Bayesian Partial Order Inference - Data Generator
# Comprehensive data generation function for partial order inference
# Author: Converted for Prof Geoff Nicholls, University of Oxford

# Required libraries
library(yaml)
library(jsonlite)
library(mvtnorm)

#' Load configuration from YAML file
#' 
#' @param config_path Path to configuration file
#' @return List containing configuration
load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(paste("Config file not found at:", config_path))
  }
  
  tryCatch({
    yaml::read_yaml(config_path)
  }, error = function(e) {
    stop(paste("Error parsing YAML file:", e$message))
  })
}

#' Generate synthetic data for partial order inference
#' 
#' @param config Configuration list containing generation parameters
#' @return List containing synthetic data
#' @export
generate_data <- function(config) {
  tryCatch({
    # 1. Set up parameters for partial order generation
    n <- config$generation$n  # Number of nodes/items in the partial order
    N <- config$generation$N  # Number of total orders to sample from the partial order
    
    # 2. Load configuration parameters
    K <- config$mcmc$K  # Number of dimensions for latent positions
    rho_prior <- config$prior$rho_prior  # Prior parameter for correlation
    noise_option <- config$noise$noise_option  # Noise model specification
 #   mallow_ua <- config$prior$mallow_ua  # Mallows model parameter
    
    items <- 0:(n-1)  # Create list of item indices to represent the items
    
    # 3. Generate true correlation parameter from prior
    rho_true <- rbeta(1, 1, rho_prior)  # Sample from Beta(1, rho_prior) distribution
    cat("True correlation parameter (rho):", sprintf("%.4f", rho_true), "\n")
    
    # 4. Set up covariates and regression parameters
    # X: p × n design matrix (p is number of covariates)
    # β: p × 1 vector of regression coefficients
    p <- config$covariates$p  # Number of covariates
    beta_true <- config$covariates$beta_true   # True regression coefficients
    X <- matrix(rnorm(p * n), nrow = p, ncol = n)  # Example design matrix
    
    # 5. Compute assessor-specific effects
    # α = X^T β is an n × 1 vector of assessor effects
    alpha <- as.vector(t(X) %*% beta_true)
    cat("\nThe covariates effects (alpha):\n")
    print(alpha)
    
    # 6. Generate latent positions for each assessor
    # U matrix represents the base latent positions in K-dimensional space
    Sigma_true <- build_sigma_rho(K, rho_true)
    U <- rmvnorm(n, mean = rep(0, K), sigma = Sigma_true)
    cat("\nBase U matrix (latent positions):\n")
    print(U)
    
    # 7. Generate latent positions with covariates effects
    # η = transform_U_to_eta(U, α) 
    eta <- transform_U_to_eta(U, alpha)
    cat("\nAdjusted latent positions (eta):\n")
    print(eta)
    
    cat("\nRegression Information:\n")
    cat("Design Matrix (X):\n")
    print(X)
    cat("\nTrue Regression Coefficients (beta):\n")
    print(beta_true)
    cat("\nCovariates Effects (alpha = X^T β):\n")
    print(alpha)
    
    # 8. Generate partial order from adjusted latent positions
    # First generate the full partial order
    h <- generate_partial_order(eta)
    # Then compute its transitive reduction to remove redundant edges
    h_true <- transitive_reduction(h)
    cat("\nPartial Order (adjacency matrix):\n")
    print(h_true)
    
    # 9. Print descriptive statistics of the generated partial order
    cat("\nPartial Order Statistics:\n")
    cat("Number of items:", n, "\n")
    cat("Number of covariates:", p, "\n")
    cat("Number of direct relationships:", sum(h_true), "\n")
    cat("Number of linear extensions:", count_linear_extensions(h_true), "\n")
    
    # 10. Generate subsets for sampling total orders
    subsets <- list()
    if(min_sub!=n){
      for (i in 1:N) {
        # Random choice set (subset of items)
        choice_size <- sample(min_sub:n, 1)
        choice_set <- sample(items, choice_size)
        subsets[[i]] <- choice_set
      }
    } else {subsets = rep(list(items),N)} 
    
    # 11. Generate total orders from subsets
    observed_orders <- list()
    choice_sets <- list()
    
    for (i in seq_len(N)) {
      choice_set       <- subsets[[i]]
      choice_sets[[i]] <- paste0("Item_", choice_set)      # keep as strings
      
      if (noise_option == "queue_jump") {
        order <- generate_total_order_queue_jump(
          subset     = choice_set,
          items_all  = items,
          h_global   = h_true,
          prob_noise = prob_noise_true
          
        )
      }
      observed_orders[[i]] <- paste0("Item_", order_global)
    }
    
    # 12. Prepare output data in the format expected by the inference module
    output_data <- list(
      observed_orders = observed_orders,  # List of observed total orders
      choice_sets = choice_sets,  # List of choice sets
      items = list(names = paste0("Item_", items)),  # Item names
      parameters = list(
        n = n,
        N = N,
        K = K,
        rho_true = rho_true,
        prob_noise_true = ifelse(is.null(config$generation$prob_noise_true), 0.1, config$generation$prob_noise_true)
      ),
      true_partial_order = h_true,  # True partial order matrix
      beta_true = beta_true,  # True covariate effects
      X = X,  # Covariate matrix
      U_true = U,  # True latent positions
      alpha_true = alpha,  # True covariate effects
      eta_true = eta  # True transformed latent positions
    )
    
    return(output_data)
    
  }, error = function(e) {
    cat("Error in generate_data:", e$message, "\n")
    stop(e)
  })
}

#' Main function for data generation
#' 
#' @param config_path Path to configuration file
#' @param output_dir Output directory for generated data
#' @return Path to saved data file
#' @export
generate_data_main <- function(config_path = "config/data_generator_config.yaml", 
                               output_dir = "output") {
  tryCatch({
    # Load configuration
    config <- load_config(config_path)
    
    # Create output directory
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Generate data
    data <- generate_data(config)
    
    # Save data
    output_path <- file.path(output_dir, "synthetic_data.json")
    jsonlite::write_json(data, output_path, pretty = TRUE, auto_unbox = TRUE)
    
    cat("\nData saved to", output_path, "\n")
    
    # Print statistics
    cat("\nSampling Statistics:\n")
    observed_orders <- data$observed_orders
    unique_orders <- unique(lapply(observed_orders, function(x) paste(x, collapse = ",")))
    cat("Number of unique total orders:", length(unique_orders), "\n")
    
    # Count occurrences of each ordering
    order_strings <- sapply(observed_orders, function(x) paste(x, collapse = ","))
    order_counts <- table(order_strings)
    
    # Find most common ordering
    most_common <- names(order_counts)[which.max(order_counts)]
    most_common_count <- max(order_counts)
    cat("Most common ordering:", most_common, "(appears", most_common_count, "times)\n")
    
    return(output_path)
    
  }, error = function(e) {
    cat("Error in generate_data_main:", e$message, "\n")
    stop(e)
  })
} 