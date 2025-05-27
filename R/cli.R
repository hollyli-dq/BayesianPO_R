#!/usr/bin/env Rscript
# Command-line interface for Bayesian Partial Order Inference
# Author: Converted for Prof Geoff Nicholls, University of Oxford

# Required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(jsonlite)
})

# Source required functions
source("R/utilities.R")
source("R/mcmc.R")
source("R/mcmc_rj.R")
source("R/analysis.R")
source("R/data_generator.R")
source("R/po_inference.R")

#' Parse command-line arguments
#' 
#' @return List containing parsed arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(c("--data-config"), type = "character", 
                default = "config/data_generator_config.yaml",
                help = "Path to data generation config file [default: %default]"),
    
    make_option(c("--mcmc-config"), type = "character", 
                default = "config/mcmc_config.yaml",
                help = "Path to MCMC config file [default: %default]"),
    
    make_option(c("--generate-data"), action = "store_true", 
                default = FALSE,
                help = "Generate synthetic data only"),
    
    make_option(c("--inference-only"), action = "store_true", 
                default = FALSE,
                help = "Run inference only with existing data"),
    
    make_option(c("--use-rj-mcmc"), action = "store_true", 
                default = FALSE,
                help = "Use reversible jump MCMC for variable dimension inference"),
    
    make_option(c("--output-dir"), type = "character", 
                default = NULL,
                help = "Output directory for results (overrides config setting)"),
    
    make_option(c("--iterations"), type = "integer", 
                default = NULL,
                help = "Number of MCMC iterations (overrides config setting)"),
    
    make_option(c("--burn-in"), type = "integer", 
                default = NULL,
                help = "Burn-in period for MCMC (overrides config setting)"),
    
    make_option(c("--latent-dim"), type = "integer", 
                default = NULL,
                help = "Latent dimension K (overrides config setting)"),
    
    make_option(c("--random-seed"), type = "integer", 
                default = NULL,
                help = "Random seed for reproducibility (overrides config setting)"),
    
    make_option(c("--verbose"), action = "store_true", 
                default = FALSE,
                help = "Enable verbose output")
  )
  
  parser <- OptionParser(
    option_list = option_list,
    description = "Bayesian Partial Order Inference",
    epilogue = paste(
      "Examples:",
      "  Rscript R/cli.R --generate-data",
      "  Rscript R/cli.R --inference-only --use-rj-mcmc",
      "  Rscript R/cli.R --iterations 5000 --burn-in 1000",
      "",
      sep = "\n"
    )
  )
  
  return(parse_args(parser))
}

#' Save generated data to JSON file
#' 
#' @param data Generated data list
#' @param output_dir Output directory
#' @param data_name Base name for data file
#' @return Path to saved data file
save_generated_data <- function(data, output_dir, data_name) {
  tryCatch({
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save data to JSON file
    data_path <- file.path(output_dir, paste0(data_name, ".json"))
    jsonlite::write_json(data, data_path, pretty = TRUE, auto_unbox = TRUE)
    cat("\nGenerated data saved to", data_path, "\n")
    
    return(data_path)
    
  }, error = function(e) {
    cat("Error in save_generated_data:", e$message, "\n")
    stop(e)
  })
}

#' Main function to run the command-line interface
#' 
#' @return Exit code (0 for success, 1 for error)
main <- function() {
  # Parse command-line arguments
  args <- parse_arguments()
  
  tryCatch({
    # Set up paths
    if (!file.exists(args$`data-config`)) {
      stop(paste("Data generator config not found at:", args$`data-config`))
    }
    if (!file.exists(args$`mcmc-config`)) {
      stop(paste("MCMC config not found at:", args$`mcmc-config`))
    }
    
    # Load configurations
    data_gen_config <- load_config(args$`data-config`)
    mcmc_config <- load_config(args$`mcmc-config`)
    
    if (args$verbose) {
      cat("Loaded data generation config from:", args$`data-config`, "\n")
      cat("Loaded MCMC config from:", args$`mcmc-config`, "\n")
    }
    
    # Override config settings with command-line arguments
    if (!is.null(args$`output-dir`)) {
      if (is.null(mcmc_config$data)) mcmc_config$data <- list()
      mcmc_config$data$output_dir <- args$`output-dir`
    }
    
    if (!is.null(args$iterations)) {
      mcmc_config$mcmc$num_iterations <- args$iterations
    }
    
    if (!is.null(args$`burn-in`)) {
      if (is.null(mcmc_config$visualization)) mcmc_config$visualization <- list()
      mcmc_config$visualization$burn_in <- args$`burn-in`
    }
    
    if (!is.null(args$`latent-dim`)) {
      mcmc_config$generation$K <- args$`latent-dim`
    }
    
    if (!is.null(args$`random-seed`)) {
      mcmc_config$mcmc$random_seed <- args$`random-seed`
    }
    
    # Set up output directory (default to results)
    output_dir <- if (!is.null(mcmc_config$data$output_dir)) {
      mcmc_config$data$output_dir
    } else {
      "results"
    }
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    cat("Output directory:", output_dir, "\n")
    
    # Generate data if specified
    if (args$`generate-data` || (!args$`inference-only` && 
                                 (is.null(mcmc_config$data$generate_data) || mcmc_config$data$generate_data))) {
      cat("\nGenerating synthetic data...\n")
      data <- generate_data(data_gen_config)
      
      # Save generated data
      data_name <- "synthetic_data"
      data_path <- save_generated_data(data, output_dir, data_name)
      
      # Update mcmc config with the path to generated data
      mcmc_config$data$path <- file.path(output_dir, paste0(data_name, ".json"))
      
      if (args$`generate-data`) {
        cat("Data generation completed successfully.\n")
        return(0)
      }
    }
    
    # Run inference if specified
    if (!args$`generate-data`) {
      # Load existing data
      data_path <- mcmc_config$data$path
      if (!file.exists(data_path)) {
        stop(paste("Data file not found at:", data_path))
      }
      
      cat("Loading data from:", data_path, "\n")
      data <- jsonlite::fromJSON(data_path, simplifyVector = FALSE)
      data_name <- tools::file_path_sans_ext(basename(data_path))
      
      # Run inference
      mcmc_type <- if (args$`use-rj-mcmc`) "Reversible Jump" else "Fixed Dimension"
      cat("\nRunning", mcmc_type, "MCMC inference...\n")
      
      if (args$verbose) {
        cat("MCMC parameters:\n")
        cat("  Iterations:", mcmc_config$mcmc$num_iterations, "\n")
        cat("  Burn-in:", mcmc_config$visualization$burn_in, "\n")
        cat("  Latent dimensions:", mcmc_config$generation$K, "\n")
        cat("  Random seed:", mcmc_config$mcmc$random_seed, "\n")
      }
      
      results <- run_inference(data, mcmc_config, args$`use-rj-mcmc`)
      
      # Save results
      save_results(results, output_dir, data_name)
      
      # Generate plots
      generate_plots(results, data, mcmc_config, output_dir, data_name)
      
      cat("\nInference completed successfully!\n")
      
      # Print summary statistics
      if (!is.null(results$h)) {
        cat("\nSummary Statistics:\n")
        cat("  Inferred edges:", sum(results$h), "\n")
        
        if (!is.null(data$true_partial_order)) {
          true_po <- if (is.list(data$true_partial_order)) {
            do.call(rbind, data$true_partial_order)
          } else {
            as.matrix(data$true_partial_order)
          }
          
          accuracy <- sum(results$h == true_po) / length(true_po)
          cat("  True edges:", sum(true_po), "\n")
          cat("  Accuracy:", round(accuracy * 100, 1), "%\n")
        }
        
        if (args$`use-rj-mcmc` && !is.null(results$K_trace)) {
          final_K <- results$K_trace[length(results$K_trace)]
          cat("  Final K:", final_K, "\n")
        }
        
        if (!is.null(results$overall_acceptance_rate)) {
          cat("  Acceptance rate:", round(results$overall_acceptance_rate * 100, 1), "%\n")
        }
      }
    }
    
    return(0)
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    if (args$verbose) {
      cat("Traceback:\n")
      traceback()
    }
    return(1)
  })
}

# Run main function if script is executed directly
if (!interactive()) {
  quit(status = main())
} 