# Bayesian Partial Order Inference - Inference Engine
# Comprehensive inference function for partial order inference
# Author: Converted for Prof Geoff Nicholls, University of Oxford

# Required libraries
library(yaml)
library(jsonlite)
library(mvtnorm)

#' Load data from JSON file
#' 
#' @param data_path Path to data file
#' @return List containing data
load_data <- function(data_path) {
  if (!file.exists(data_path)) {
    stop(paste("Data file not found at:", data_path))
  }
  
  tryCatch({
    jsonlite::fromJSON(data_path, simplifyVector = FALSE)
  }, error = function(e) {
    stop(paste("Error parsing JSON file:", e$message))
  })
}

#' Run MCMC inference on the data
#' 
#' @param data List containing observed data
#' @param config Configuration list containing MCMC parameters
#' @param use_rj_mcmc Whether to use reversible jump MCMC (default: FALSE)
#' @return List containing MCMC results
#' @export
run_inference <- function(data, config, use_rj_mcmc = FALSE) {
  tryCatch({
    # Extract data
    observed_orders <- data$observed_orders
    choice_sets <- data$choice_sets
    parameters <- data$parameters
    
    # Get MCMC parameters from config
    num_iterations <- config$mcmc$num_iterations
    K <- config$generation$K
    
    # Get update probabilities
    if (use_rj_mcmc) {
      mcmc_pt <- c(
        config$mcmc$update_probabilities$rho,
        config$mcmc$update_probabilities$noise,
        config$mcmc$update_probabilities$U,
        config$mcmc$update_probabilities$beta,
        config$mcmc$update_probabilities$K
      )
    } else {
      mcmc_pt <- c(
        config$mcmc$update_probabilities$rho,
        config$mcmc$update_probabilities$noise,
        config$mcmc$update_probabilities$U,
        config$mcmc$update_probabilities$beta
      )
    }
    
    dr <- config$rho$dr
    drbeta <- config$mcmc$drbeta
    noise_option <- config$noise$noise_option
    sigma_mallow <- config$noise$sigma_mallow
    sigma_beta <- config$mcmc$sigma_beta
    
    # Get prior parameters
    rho_prior <- config$prior$rho_prior
    noise_beta_prior <- config$prior$noise_beta_prior
    mallow_ua <- config$prior$mallow_ua
    
    # Get covariate effects - handle JSON conversion properly
    beta_true <- if (!is.null(data$beta_true)) {
      if (is.list(data$beta_true)) unlist(data$beta_true) else data$beta_true
    } else {
      rep(0, config$covariates$p)
    }
    
    X_data <- if (!is.null(data$X)) {
      if (is.list(data$X)) {
        # Convert list of lists to matrix
        do.call(rbind, lapply(data$X, function(x) if (is.list(x)) unlist(x) else x))
      } else {
        as.matrix(data$X)
      }
    } else {
      matrix(0, nrow = config$covariates$p, ncol = parameters$n)
    }
    
    # Ensure X is in the correct format (p x n)
    X <- as.matrix(X_data)
    if (nrow(X) != config$covariates$p) {
      stop(paste("X matrix has", nrow(X), "rows but expected", config$covariates$p))
    }
    
    # Run MCMC simulation
    if (use_rj_mcmc) {
      # Reversible Jump MCMC
      K_prior <- config$prior$K_prior
      mcmc_results <- mcmc_partial_order_k(
        observed_orders = observed_orders,
        choice_sets = choice_sets,
        num_iterations = num_iterations,
        X = X,
        dr = dr,
        drbeta = drbeta,
        sigma_mallow = sigma_mallow,
        sigma_beta = sigma_beta,
        noise_option = noise_option,
        mcmc_pt = mcmc_pt,
        rho_prior = rho_prior,
        noise_beta_prior = noise_beta_prior,
        mallow_ua = mallow_ua,
        K_prior = K_prior,
        random_seed = config$mcmc$random_seed
      )
    } else {
      # Fixed dimension MCMC
      mcmc_results <- mcmc_partial_order(
        observed_orders = observed_orders,
        choice_sets = choice_sets,
        num_iterations = num_iterations,
        K = K,
        X = X,
        dr = dr,
        drbeta = drbeta,
        sigma_mallow = sigma_mallow,
        sigma_beta = sigma_beta,
        noise_option = noise_option,
        mcmc_pt = mcmc_pt,
        rho_prior = rho_prior,
        noise_beta_prior = noise_beta_prior,
        mallow_ua = mallow_ua,
        random_seed = config$mcmc$random_seed
      )
    }
    
    # Compute the final inferred partial order 'h' from the MCMC trace
    if (!is.null(mcmc_results$h_trace) && length(mcmc_results$h_trace) > 0) {
      burn_in <- config$visualization$burn_in
      burn_in_index <- max(1, burn_in %/% 100)  # Convert to trace index
      
      if (burn_in_index >= length(mcmc_results$h_trace)) {
        cat("Warning: burn_in (", burn_in, ") is larger than trace length (", 
            length(mcmc_results$h_trace), "). Using last available samples.\n")
        burn_in_index <- max(1, length(mcmc_results$h_trace) - 10)
      }
      
      post_burn_in_trace <- mcmc_results$h_trace[(burn_in_index + 1):length(mcmc_results$h_trace)]
      
      if (length(post_burn_in_trace) == 0) {
        stop("No valid data after burn-in period")
      }
      
      # Convert list of matrices to 3D array for easier averaging
      n_items <- nrow(post_burn_in_trace[[1]])
      h_array <- array(0, dim = c(n_items, n_items, length(post_burn_in_trace)))
      
      for (i in seq_along(post_burn_in_trace)) {
        h_array[,,i] <- post_burn_in_trace[[i]]
      }
      
      # Compute the mean over the trace
      h_final <- apply(h_array, c(1,2), mean)
      
      # Apply a threshold (e.g. 0.5) and perform transitive reduction
      threshold <- 0.5
      h_final_inferred <- transitive_reduction((h_final >= threshold) * 1)
      
      # Add the final inferred partial order to the results
      mcmc_results$h <- h_final_inferred
      mcmc_results$h_posterior_probs <- h_final
    } else {
      cat("Warning: h_trace not found in MCMC results; setting 'h' to an empty matrix.\n")
      mcmc_results$h <- matrix(0, nrow = parameters$n, ncol = parameters$n)
      mcmc_results$h_posterior_probs <- matrix(0, nrow = parameters$n, ncol = parameters$n)
    }
    
    # Add final states for other parameters
    if (!is.null(mcmc_results$Z_trace) && length(mcmc_results$Z_trace) > 0) {
      mcmc_results$Z <- mcmc_results$Z_trace[[length(mcmc_results$Z_trace)]]
    } else {
      mcmc_results$Z <- matrix(0, nrow = parameters$n, ncol = K)
    }
    
    mcmc_results$beta <- beta_true
    
    if (!is.null(mcmc_results$rho_trace) && length(mcmc_results$rho_trace) > 0) {
      mcmc_results$rho <- mcmc_results$rho_trace[length(mcmc_results$rho_trace)]
    } else {
      mcmc_results$rho <- 0.0
    }
    
    if (!is.null(mcmc_results$prob_noise_trace) && length(mcmc_results$prob_noise_trace) > 0) {
      mcmc_results$prob_noise <- mcmc_results$prob_noise_trace[length(mcmc_results$prob_noise_trace)]
    } else {
      mcmc_results$prob_noise <- 0.0
    }
    
    if (!is.null(mcmc_results$mallow_theta_trace) && length(mcmc_results$mallow_theta_trace) > 0) {
      mcmc_results$mallow_theta <- mcmc_results$mallow_theta_trace[length(mcmc_results$mallow_theta_trace)]
    } else {
      mcmc_results$mallow_theta <- 1.0
    }
    
    # Package trace information
    trace_keys <- c('Z_trace', 'h_trace', 'rho_trace', 'prob_noise_trace', 'mallow_theta_trace')
    if (use_rj_mcmc) {
      trace_keys <- c(trace_keys, 'K_trace')
    }
    
    trace_info <- list()
    for (key in trace_keys) {
      if (!is.null(mcmc_results[[key]])) {
        trace_info[[key]] <- mcmc_results[[key]]
      } else {
        trace_info[[key]] <- list()
      }
    }
    mcmc_results$trace <- trace_info
    
    return(mcmc_results)
    
  }, error = function(e) {
    cat("Error in run_inference:", e$message, "\n")
    stop(e)
  })
}

#' Save inference results to JSON file
#' 
#' @param results List containing MCMC results
#' @param output_dir Output directory
#' @param data_name Base name for output files
#' @export
save_results <- function(results, output_dir, data_name) {
  tryCatch({
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Convert results to JSON-serializable format
    results_serializable <- results
    
    # Convert matrices and arrays to lists
    for (name in names(results_serializable)) {
      if (is.matrix(results_serializable[[name]]) || is.array(results_serializable[[name]])) {
        results_serializable[[name]] <- as.list(results_serializable[[name]])
      }
    }
    
    # Handle trace information
    if (!is.null(results_serializable$trace)) {
      for (trace_name in names(results_serializable$trace)) {
        trace_data <- results_serializable$trace[[trace_name]]
        if (is.list(trace_data) && length(trace_data) > 0) {
          # Convert matrices in trace to lists
          results_serializable$trace[[trace_name]] <- lapply(trace_data, function(x) {
            if (is.matrix(x)) as.list(x) else x
          })
        }
      }
    }
    
    # Save results to JSON file
    results_path <- file.path(output_dir, paste0(data_name, "_results.json"))
    jsonlite::write_json(results_serializable, results_path, pretty = TRUE, auto_unbox = TRUE)
    cat("\nResults saved to", results_path, "\n")
    
    # Save partial order matrix separately as RDS file
    if (!is.null(results$h)) {
      h_path <- file.path(output_dir, paste0(data_name, "_partial_order.rds"))
      saveRDS(results$h, h_path)
      cat("Partial order matrix saved to", h_path, "\n")
    }
    
  }, error = function(e) {
    cat("Error in save_results:", e$message, "\n")
    stop(e)
  })
}

#' Generate and save plots
#' 
#' @param results List containing MCMC results
#' @param data List containing original data
#' @param config Configuration list
#' @param output_dir Output directory
#' @param data_name Base name for output files
#' @export
generate_plots <- function(results, data, config, output_dir, data_name) {
  tryCatch({
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Get true parameters if they exist
    true_params <- list()
    if (!is.null(data$parameters)) {
      if (!is.null(data$parameters$rho_true)) {
        true_params$rho_true <- data$parameters$rho_true
      }
      if (!is.null(data$parameters$prob_noise_true)) {
        true_params$prob_noise_true <- data$parameters$prob_noise_true
      }
    }
    
    # Get burn-in from config
    burn_in <- config$visualization$burn_in
    
    # Plot MCMC results using the enhanced plotting function
    if (!is.null(results$trace)) {
      plot_mcmc_results(results, true_params, config, burn_in = burn_in %/% 100)
    }
    
    # Get item names
    items <- if (!is.null(data$items$names)) {
      data$items$names
    } else {
      paste0("Item_", 0:(nrow(results$h) - 1))
    }
    
    # Plot inferred partial order
    if (!is.null(results$h)) {
      cat("Plotting inferred partial order...\n")
      plot_partial_order(results$h, items, 
                        title = "Inferred Partial Order (Posterior Mean)",
                        vertex_color = "lightblue", 
                        vertex_frame_color = "darkblue",
                        edge_color = "darkblue")
    }
    
    # Plot true partial order if available
    if (!is.null(data$true_partial_order)) {
      cat("Plotting true partial order...\n")
      true_po <- if (is.list(data$true_partial_order)) {
        do.call(rbind, data$true_partial_order)
      } else {
        as.matrix(data$true_partial_order)
      }
      
      plot_partial_order(true_po, items, 
                        title = "True Partial Order",
                        vertex_color = "lightgray", 
                        vertex_frame_color = "darkgray",
                        edge_color = "darkblue")
      
      # Compare relationships
      if (!is.null(results$h)) {
        missing_relationships <- compare_partial_orders(true_po, results$h, "missing")
        redundant_relationships <- compare_partial_orders(true_po, results$h, "redundant")
        
        # Print relationship comparisons
        if (length(missing_relationships) > 0) {
          cat("\nMissing Relationships (edges present in true PO but absent in inferred PO):\n")
          for (rel in missing_relationships) {
            cat(items[rel[1]], "<", items[rel[2]], "\n")
          }
        } else {
          cat("\nNo missing relationships. The inferred partial order captures all true relationships.\n")
        }
        
        if (length(redundant_relationships) > 0) {
          cat("\nRedundant Relationships (edges present in inferred PO but absent in true PO):\n")
          for (rel in redundant_relationships) {
            cat(items[rel[1]], "<", items[rel[2]], "\n")
          }
        } else {
          cat("\nNo redundant relationships. The inferred partial order is a subset of the true partial order.\n")
        }
        
        # Calculate accuracy metrics
        accuracy <- sum(results$h == true_po) / length(true_po)
        cat("\nOverall accuracy:", round(accuracy * 100, 1), "%\n")
        cat("True edges:", sum(true_po), ", Inferred edges:", sum(results$h), "\n")
      }
    }
    
  }, error = function(e) {
    cat("Error in generate_plots:", e$message, "\n")
    stop(e)
  })
}

#' Compare two partial orders and find differences
#' 
#' @param true_po True partial order matrix
#' @param inferred_po Inferred partial order matrix
#' @param type Type of comparison ("missing" or "redundant")
#' @return List of edge differences
compare_partial_orders <- function(true_po, inferred_po, type = "missing") {
  differences <- list()
  
  for (i in 1:nrow(true_po)) {
    for (j in 1:ncol(true_po)) {
      if (type == "missing" && true_po[i,j] == 1 && inferred_po[i,j] == 0) {
        differences <- append(differences, list(c(i, j)))
      } else if (type == "redundant" && true_po[i,j] == 0 && inferred_po[i,j] == 1) {
        differences <- append(differences, list(c(i, j)))
      }
    }
  }
  
  return(differences)
}

#' Main function to run the inference pipeline
#' 
#' @param data_path Path to data file
#' @param config_path Path to configuration file
#' @param output_dir Output directory
#' @param use_rj_mcmc Whether to use reversible jump MCMC
#' @export
inference_main <- function(data_path, config_path, output_dir = "output", use_rj_mcmc = FALSE) {
  tryCatch({
    cat("Loading configuration from:", config_path, "\n")
    config <- load_config(config_path)
    
    cat("Loading data from:", data_path, "\n")
    data <- load_data(data_path)
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Run inference
    cat("Running", ifelse(use_rj_mcmc, "Reversible Jump", "Fixed Dimension"), "MCMC inference...\n")
    results <- run_inference(data, config, use_rj_mcmc)
    
    # Save results
    data_name <- tools::file_path_sans_ext(basename(data_path))
    save_results(results, output_dir, data_name)
    
    # Generate plots
    generate_plots(results, data, config, output_dir, data_name)
    
    cat("\nInference completed successfully!\n")
    
    return(results)
    
  }, error = function(e) {
    cat("Error in inference_main:", e$message, "\n")
    stop(e)
  })
} 