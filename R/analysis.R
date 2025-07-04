# Bayesian Partial Order Inference - Analysis and Visualization
# Main analysis functions, visualization, and example usage
# Author: Converted for Prof Geoff Nicholls, University of Oxford

# Required libraries
library(ggplot2)
library(igraph)
library(reshape2)

# ============================================================================
# VSP-STYLE VISUALIZATION FUNCTIONS
# ============================================================================

#' Show DAG (Directed Acyclic Graph) from adjacency matrix
#' 
#' @param m Adjacency matrix
#' @param vertex_color Vertex fill color
#' @param vertex_frame_color Vertex border color
#' @param edge_color Edge color
#' @param ... Additional plotting parameters
#' @export
showDAG <- function(m = NULL, vertex_color = "white", vertex_frame_color = "black", 
                   edge_color = "black", ...) {
  # plots a partial order from an adjacency matrix
  g <- graph_from_adjacency_matrix(m, mode = "directed")
  h <- g  # seems like a bug in layout function puts vertices on top of one another
  plot(g, layout = layout_with_sugiyama(h)$layout, 
       vertex.color = vertex_color,
       vertex.frame.color = vertex_frame_color,
       vertex.label.color = "black",
       vertex.label.family = "serif",
       edge.color = edge_color,
       edge.width = 1.5,
       main.family = "serif",
       ...)
}

#' Transitive reduction function
#' 
#' @param g Adjacency matrix
#' @return Transitive reduction matrix
#' @export
transitive_reduction_vsp <- function(g) {
  if (!(class(g)[1] %in% c("matrix", "graphNEL"))) 
    stop("Input must be an adjacency matrix or graphNEL object")
  if (class(g)[1] == "graphNEL") {
    g <- as(g, "matrix")
  }
  g <- transitive_closure_vsp(g, mat = TRUE)
  g <- g - diag(diag(g))
  type <- (g > 1) * 1 - (g < 0) * 1
  for (y in 1:nrow(g)) {
    for (x in 1:nrow(g)) {
      if (g[x, y] != 0) {
        for (j in 1:nrow(g)) {
          if ((g[y, j] != 0) & sign(type[x, j]) * sign(type[x, y]) * sign(type[y, j]) != -1) {
            g[x, j] <- 0
          }
        }
      }
    }
  }
  g
}

#' Transitive closure function
#' 
#' @param g Adjacency matrix
#' @param mat Return as matrix
#' @param loops Include loops
#' @return Transitive closure matrix
#' @export
transitive_closure_vsp <- function(g, mat = FALSE, loops = TRUE) {
  if (!(class(g)[1] %in% c("graphNEL", "matrix")))         
    stop("Input must be either graphNEL object or adjacency matrix")    
  g <- as(g, "matrix")    
  n <- ncol(g)    
  matExpIterativ <- function(x, pow, y = x, z = x, i = 1) {        
    while (i < pow) {            
      z <- z %*% x            
      y <- y + z            
      i <- i + 1
    }        
    return(y)
  }    
  h <- matExpIterativ(g, n)    
  h <- (h > 0) * 1    
  dimnames(h) <- dimnames(g)    
  if (!loops)         
    diag(h) <- rep(0, n)    
  else diag(h) <- rep(1, n)    
  if (!mat)         
    h <- as(h, "graphNEL")    
  return(h)
}

#' Calculate DAG depth from transitive reduction
#' 
#' @param tr Transitive reduction matrix
#' @return DAG depth
#' @export
dagdepth <- function(tr) {
  if (length(tr) == 1) {return(1)}
  n <- dim(tr)[1]
  cs <- apply(tr, 2, sum)
  if (sum(cs) == 0) {return(1)}
  csi <- (cs == 0)
  bs <- apply(tr, 1, sum)
  bsi <- (bs == 0)
  free <- which(bsi & csi)
  k <- length(free)
  if (k == n) {return(1)}
  if (k > 0) { 
    tr <- tr[-free, -free]
    cs <- apply(tr, 2, sum)
    csi <- (cs == 0)
    bs <- apply(tr, 1, sum)
    bsi <- (bs == 0)
  }
  tops <- which(csi)
  bots <- which(bsi)
  if (length(bots) > length(tops)) {tops <- bots}
  return(1 + dagdepth(tr[-tops, -tops]))
}

#' Show consensus DAG with probability thresholds
#' 
#' @param con_po Consensus partial order matrix
#' @param threshold Threshold for edge inclusion
#' @param label_threshold Threshold for edge labeling
#' @param size Vertex size
#' @export
showDAGcon <- function(con_po, threshold = 0.5, label_threshold = 0.9, size = 10) {
  con_po1 <- (con_po > threshold) + 0.0
  con_po2 <- (con_po > label_threshold) + 0.0
  con_po1 <- transitive_reduction_vsp(con_po1)
  con_po2 <- transitive_reduction_vsp(con_po2)
  g <- graph_from_adjacency_matrix(con_po1, mode = "directed")
  g_label <- graph_from_adjacency_matrix(con_po2)
  el1 <- apply(get.edgelist(g), 1, paste, collapse = "-")
  el2 <- apply(get.edgelist(g_label), 1, paste, collapse = "-")
  E(g)$color <- ifelse(el1 %in% el2, "red", "darkgrey")
  g <- set_edge_attr(g, name = 'weight', value = -1)
  dist <- (-distances(g, v = (V(g)), mode = 'out'))
  layers <- apply(dist, 1, max)
  layers <- max(layers) - layers
  plot(g, layout = layout_with_sugiyama(g)$layout, vertex.size = size, 
       edge.arrow.size = 0.2, main = 'Consensus Order')
}

#' Enhanced partial order plotting
#' 
#' @param h Partial order matrix
#' @param item_names Item names
#' @param title Plot title
#' @param vertex_size Vertex size
#' @param edge_arrow_size Arrow size
#' @param vertex_color Vertex fill color
#' @param vertex_frame_color Vertex border color
#' @param vertex_label_color Vertex label color
#' @param edge_color Edge color
#' @export
plot_partial_order <- function(h,
                               item_names = NULL,
                               title = "Partial Order",
                               # --- visual tuning --------------------------------------------------
                               vertex_size         = 30,
                               edge_arrow_size     = 0.5,
                               vertex_color        = "white",
                               vertex_frame_color  = "black",
                               vertex_label_color  = "black",
                               edge_color          = "black",
                               # new knobs
                               label_cex           = 1.4,   # text size multiplier
                               frame_width         = 2,     # vertex border thickness
                               label_family        = "serif",# or "sans", "mono" …
                               layout_fun          = igraph::layout_with_sugiyama)
{
  stopifnot(is.matrix(h), nrow(h) == ncol(h))
  n <- nrow(h)
  if (is.null(item_names))
    item_names <- paste0("Item_", seq_len(n))
  
  # 1. transitive reduction for clearer picture --------------------------------
  h_red <- transitive_reduction_vsp(h)
  
  # 2. build graph -------------------------------------------------------------
  g <- igraph::graph_from_adjacency_matrix(h_red, mode = "directed")
  igraph::V(g)$name <- item_names
  
  # 3. lay out once (Sugiyama returns a list) ----------------------------------
  lay <- layout_fun(g)
  if (is.list(lay)) lay <- lay$layout   # igraph 1.4+ returns list for Sugiyama
  
  # 4. draw --------------------------------------------------------------------
  plot(g,
       layout               = lay,
       main                 = title,
       main.cex             = 1.2,               # title size
       main.family          = label_family,
       vertex.size          = vertex_size,
       vertex.color         = vertex_color,
       vertex.frame.color   = vertex_frame_color,
       vertex.frame.width   = frame_width,
       vertex.label         = igraph::V(g)$name,
       vertex.label.cex     = label_cex,
       vertex.label.family  = label_family,
       vertex.label.color   = vertex_label_color,
       edge.arrow.size      = edge_arrow_size,
       edge.color           = edge_color,
       edge.width           = 1.5)
  invisible(g)
}


#' Plot MCMC inferred variables with comprehensive diagnostics
#' 
#' @param mcmc_results MCMC results object
#' @param true_param List of true parameter values
#' @param config Configuration list with prior parameters
#' @param burn_in Burn-in period
#' @param output_filename Output filename for saving plot
#' @param output_filepath Output file path
#' @export
plot_mcmc_results <- function(mcmc_results, 
                             true_param = list(),
                             config = list(),
                             burn_in = 0,
                             output_filename = "",
                             output_filepath = ".") {
  
  # Helper function for null coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Extract main MCMC traces and apply burn-in
  traces <- list()
  true_values <- list()
  
  # Check what traces are actually available and extract them
  available_traces <- names(mcmc_results)
  
  # Map actual trace names to display names
  trace_mapping <- list(
    rho = "rho_trace",
    prob_noise = "prob_noise_trace", 
    K = "K_trace"
  )
  
  # Define variables to plot with their properties
  var_configs <- list(
    rho = list(color = '#1f77b4', prior = 'beta', 
               prior_params = list(shape1 = 1.0, shape2 = config$prior$rho_prior %||% 1.0), 
               truncated = TRUE),
    prob_noise = list(color = 'orange', prior = 'beta', 
                     prior_params = list(shape1 = 1.0, shape2 = config$prior$noise_beta_prior %||% 1.0)),
    K = list(color = 'darkcyan', prior = 'poisson', 
             prior_params = list(lambda = config$prior$K_prior %||% 1.0))
  )
  
  # Extract traces and true values for available parameters
  for (var_name in names(var_configs)) {
    trace_key <- trace_mapping[[var_name]]
    if (!is.null(trace_key) && trace_key %in% available_traces && !is.null(mcmc_results[[trace_key]])) {
      trace_data <- mcmc_results[[trace_key]]
      if (length(trace_data) > burn_in) {
        traces[[var_name]] <- trace_data[(burn_in + 1):length(trace_data)]
        true_key <- paste0(var_name, "_true")
        true_values[[var_name]] <- true_param[[true_key]]
      }
    }
  }
  
  # Create plots
  n_vars <- length(traces)
  if (n_vars == 0) {
    warning("No valid traces found for plotting")
    return()
  }
  
  # Set up plotting device
  if (output_filename != "") {
    pdf(file.path(output_filepath, output_filename), width = 12, height = 4 * n_vars)
  }
  
  par(mfrow = c(n_vars, 2), mar = c(4, 4, 3, 2))
  
  # Plot each variable
  for (var_name in names(traces)) {
    trace <- traces[[var_name]]
    var_config <- var_configs[[var_name]]
    true_val <- true_values[[var_name]]
    
    # Trace plot
    iterations <- (burn_in + 1):(burn_in + length(trace))
    
    if (is.matrix(trace) || is.array(trace)) {
      # For matrix/array traces, plot the mean
      mean_trace <- apply(trace, 1, mean)
      plot(iterations, mean_trace, type = "l", col = var_config$color, 
           lwd = 1.2, main = paste("Trace:", var_name),
           xlab = "Iteration", ylab = paste("Mean", var_name),
           family = "serif")
    } else {
      plot(iterations, trace, type = "l", col = var_config$color, 
           lwd = 1.2, main = paste("Trace:", var_name),
           xlab = "Iteration", ylab = var_name,
           family = "serif")
    }
    grid(col = "gray", lty = 3, lwd = 0.5)
    
    # Add burn-in line if burn_in > 0
    if (burn_in > 0) {
      abline(v = burn_in, col = "red", lty = 2, lwd = 1)
    }
    
    # Density plot
    if (is.matrix(trace) || is.array(trace)) {
      plot_data <- apply(trace, 1, mean)
    } else {
      plot_data <- trace
    }
    
    # Handle special case for rho with truncation
    if (var_name == "rho" && var_config$truncated) {
      tol <- 1e-4
      trunc_point <- 1 - tol
      plot_data <- plot_data[plot_data <= trunc_point]
      
      hist(plot_data, breaks = 30, freq = FALSE, col = var_config$color, 
           border = "black", main = paste("Density:", var_name),
           xlab = var_name, ylab = "Density", xlim = c(0, trunc_point),
           family = "serif")
      
      # Add theoretical prior PDF
      if (!is.null(var_config$prior_params)) {
        x_vals <- seq(0, trunc_point, length.out = 1000)
        norm_const <- pbeta(trunc_point, var_config$prior_params$shape1, var_config$prior_params$shape2)
        norm_const <- max(norm_const, 1e-15)
        prior_pdf <- dbeta(x_vals, var_config$prior_params$shape1, var_config$prior_params$shape2) / norm_const
        lines(x_vals, prior_pdf, col = "black", lwd = 2, lty = 1)
      }
      
    } else {
      hist(plot_data, breaks = 30, freq = FALSE, col = var_config$color, 
           border = "black", main = paste("Density:", var_name),
           xlab = var_name, ylab = "Density", family = "serif")
      
      # Add prior distribution if specified
      if (!is.null(var_config$prior) && !is.null(var_config$prior_params)) {
        x_range <- range(plot_data)
        
        if (var_config$prior == "beta") {
          x_vals <- seq(0, 1, length.out = 1000)
          prior_pdf <- dbeta(x_vals, var_config$prior_params$shape1, var_config$prior_params$shape2)
          # Scale to match histogram
          scale_factor <- max(hist(plot_data, breaks = 30, plot = FALSE)$density) / max(prior_pdf) * 0.8
          lines(x_vals, prior_pdf * scale_factor, col = "black", lwd = 2, lty = 2)
        } else if (var_config$prior == "poisson") {
          k_range <- seq(1, max(plot_data) + 2)
          lambda_param <- var_config$prior_params$lambda
          norm_const <- 1.0 - exp(-lambda_param)
          pmf_vals <- dpois(k_range, lambda_param) / norm_const
          # Scale to match histogram
          hist_max <- max(hist(plot_data, breaks = max(plot_data), plot = FALSE)$density)
          scale_factor <- hist_max / max(pmf_vals) * 0.8
          points(k_range, pmf_vals * scale_factor, col = "black", pch = 16, cex = 0.8)
          lines(k_range, pmf_vals * scale_factor, col = "black", lwd = 2, lty = 2)
        }
      }
    }
    
    # Add true value if available
    if (!is.null(true_val)) {
      if (is.matrix(true_val) || is.array(true_val)) {
        true_val <- mean(true_val)
      }
      abline(v = true_val, col = "red", lwd = 2, lty = 2)
    }
    
    # Add sample mean
    sample_mean <- mean(plot_data)
    abline(v = sample_mean, col = "green", lwd = 2, lty = 2)
    
    # Add legend
    legend_items <- c("Sample Mean")
    legend_colors <- c("green")
    legend_lty <- c(2)
    
    if (!is.null(true_val)) {
      legend_items <- c("True Value", legend_items)
      legend_colors <- c("red", legend_colors)
      legend_lty <- c(2, legend_lty)
    }
    
    if (!is.null(var_config$prior)) {
      legend_items <- c(legend_items, "Prior")
      legend_colors <- c(legend_colors, "black")
      legend_lty <- c(legend_lty, 2)
    }
    
    legend("topright", legend = legend_items, col = legend_colors, 
           lty = legend_lty, lwd = 2, cex = 0.8)
  }
  
  if (output_filename != "") {
    dev.off()
    cat("Saved MCMC parameter plots to '", file.path(output_filepath, output_filename), "'\n")
  }
}

#' Plot top partial orders from MCMC results
#' 
#' @param results MCMC results
#' @param burn_in Burn-in period
#' @param top_n Number of top partial orders to show
#' @export
plot_top_partial_orders <- function(results, burn_in = 0, top_n = 3) {
  post_samples <- (burn_in + 1):length(results$h_trace)
  po_posterior <- results$h_trace[post_samples]
  
  # Convert to strings for counting
  str2po <- function(x, n) {
    matrix(as.numeric(strsplit(x, " ")[[1]]), nrow = n)
  }
  
  PO_unlist <- table(sapply(po_posterior, paste, collapse = " "))
  top_pos <- sort(PO_unlist, decreasing = TRUE)[1:min(top_n, length(PO_unlist))]
  
  n <- nrow(po_posterior[[1]])
  par(mfrow = c(1, min(top_n, length(top_pos))))
  
  # Define colors for different ranks
  colors <- c("lightgray", "lightcoral", "lightblue", "lightgreen", "lightyellow")
  frame_colors <- c("darkgray", "darkred", "darkblue", "darkgreen", "orange")
  edge_colors <- c("darkblue", "darkred", "darkblue", "darkgreen", "orange")
  
  for (i in 1:min(top_n, length(top_pos))) {
    po_matrix <- str2po(names(top_pos)[i], n)
    showDAG(transitive_reduction_vsp(po_matrix), 
            main = paste("Rank", i, "VSP"),
            vertex.size = 20,
            edge.arrow.size = 0.3,
            vertex_color = colors[i],
            vertex_frame_color = frame_colors[i],
            edge_color = edge_colors[i])
  }
  
  if (top_n > 1) {
    mtext(paste("Top", min(top_n, length(top_pos)), "VSPs"), 
          side = 3, line = -2, outer = TRUE, cex = 1.5, family = "serif")
  }
  
  par(mfrow = c(1, 1))
}

#' Plot consensus partial order
#' 
#' @param results MCMC results
#' @param burn_in Burn-in period
#' @param threshold Threshold for consensus
#' @export
plot_consensus_order <- function(results, burn_in = 0, threshold = 0.5) {
  post_samples <- (burn_in + 1):length(results$h_trace)
  po_posterior <- results$h_trace[post_samples]
  
  # Calculate consensus partial order
  consensus_po <- Reduce('+', po_posterior) / length(po_posterior)
  
  par(mfrow = c(1, 1))
  showDAGcon(consensus_po, threshold = threshold)
}

# ============================================================================
# MAIN ANALYSIS FUNCTIONS
# ============================================================================

#' Comprehensive MCMC analysis and visualization
#' 
#' @param results MCMC results
#' @param burn_in Burn-in period
#' @param save_plots Whether to save plots
#' @export
analyze_mcmc_results <- function(results, burn_in = 0, save_plots = FALSE) {
  
  cat("=== MCMC Results Analysis ===\n\n")
  
  # Basic statistics
  N <- length(results$rho_trace)
  post_samples <- (burn_in + 1):N
  
  cat("Total iterations:", N, "\n")
  cat("Burn-in period:", burn_in, "\n")
  cat("Posterior samples:", length(post_samples), "\n")
  cat("Overall acceptance rate:", 
      round(results$overall_acceptance_rate * 100, 2), "%\n\n")
  
  # Parameter estimates
  cat("=== Parameter Estimates ===\n")
  cat("Rho (correlation):\n")
  cat("  Mean:", round(mean(results$rho_trace[post_samples]), 3), "\n")
  cat("  95% CI:", round(quantile(results$rho_trace[post_samples], c(0.025, 0.975)), 3), "\n")
  
  if ("prob_noise_trace" %in% names(results)) {
    cat("Noise probability:\n")
    cat("  Mean:", round(mean(results$prob_noise_trace[post_samples]), 3), "\n")
    cat("  95% CI:", round(quantile(results$prob_noise_trace[post_samples], c(0.025, 0.975)), 3), "\n")
  }
  
  # Generate plots
  cat("\n=== Generating Plots ===\n")
  
  # 1. Trace plots
  plot_mcmc_results(results, true_param = list(rho_true = 0.7), config = list(prior = list(rho_prior = 1.0)), burn_in = burn_in)
  
  # 2. Posterior distributions
  plot_posterior_distributions(results, burn_in)
  
  # 3. Top partial orders
  plot_top_partial_orders(results, burn_in, top_n = 3)
  
  # 4. Consensus order
  plot_consensus_order(results, burn_in)
  
  cat("Analysis complete!\n")
}
# EXAMPLE USAGE AND DEMONSTRATION
# ============================================================================

#' Run example analysis
#' 
#' @export
run_example_analysis <- function() {
  cat("Running example Bayesian partial order analysis...\n")
  
  # Generate synthetic data
  synthetic_data <- generate_synthetic_data(
    n_items = 5,
    n_observations = 30,
    K = 2,
    rho_true = 0.7,
    prob_noise_true = 0.1,
    random_seed = 42
  )
  
  # Run MCMC
  mcmc_results <- mcmc_partial_order(
    observed_orders = synthetic_data$observed_orders,
    choice_sets = synthetic_data$choice_sets,
    num_iterations = 1000,
    K = 2,
    X = t(matrix(rnorm(2 * 5), nrow = 2)),
    random_seed = 42
  )
  
  # Analyze results
  analyze_mcmc_results(mcmc_results, burn_in = 200)
  
  return(mcmc_results)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x 
#' Compare MCMC results for a specific parameter between fixed and RJ models
#'
#' Compare MCMC results for a specific parameter between fixed and RJ models
#'
#' @param mcmc_fixed MCMC results from fixed-dimension model (list with traces)
#' @param mcmc_rj MCMC results from reversible jump model (list with traces)
#' @param param Name of parameter to analyze ("rho", "prob_noise", "tau", etc.)
#' @param burn_in Burn-in period to remove from traces (default = 0)
#' @param true_value Optional true parameter value to display (default = NULL)
#' @param trace_name_fixed Name of trace in mcmc_fixed (default = paste0(param, "_trace"))
#' @param trace_name_rj Name of trace in mcmc_rj (default = paste0(param, "_trace"))
#' @param output_filename Output filename for saving plot (default = ""). Leave empty to display plot
#' @param output_filepath Output directory path (default = ".")
#' @param colors Vector of two colors for fixed/RJ results (default = c("blue", "red"))
#' @export
mcmc_result_compare_parameter_analysis <- function(
    mcmc_fixed, 
    mcmc_rj, 
    param, 
    burn_in = 0,
    true_value = NULL,
    trace_name_fixed = NULL,
    trace_name_rj = NULL,
    output_filename = "",
    output_filepath = ".",
    colors = c("blue", "red")
) {
  # Helper for null coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Determine trace names
  trace_name_fixed <- trace_name_fixed %||% paste0(param, "_trace")
  trace_name_rj <- trace_name_rj %||% paste0(param, "_trace")
  
  # Validate traces exist
  if (!trace_name_fixed %in% names(mcmc_fixed)) {
    stop("Fixed trace '", trace_name_fixed, "' not found in mcmc_fixed")
  }
  if (!trace_name_rj %in% names(mcmc_rj)) {
    stop("RJ trace '", trace_name_rj, "' not found in mcmc_rj")
  }
  
  # Extract traces
  trace_fixed_orig <- mcmc_fixed[[trace_name_fixed]]
  trace_rj_orig <- mcmc_rj[[trace_name_rj]]
  
  # Apply burn-in
  extract_post_burnin <- function(trace, burn_in) {
    if (burn_in > 0) {
      if (is.matrix(trace)) {
        trace[(burn_in + 1):nrow(trace), , drop = FALSE]
      } else {
        trace[(burn_in + 1):length(trace)]
      }
    } else {
      trace
    }
  }
  
  trace_fixed <- extract_post_burnin(trace_fixed_orig, burn_in)
  trace_rj <- extract_post_burnin(trace_rj_orig, burn_in)
  
  # Process vector parameters
  process_vector_trace <- function(trace) {
    if (is.matrix(trace) && ncol(trace) > 1) {
      rowMeans(trace)
    } else if (is.matrix(trace)) {
      as.vector(trace)
    } else {
      trace
    }
  }
  
  trace_fixed <- process_vector_trace(trace_fixed)
  trace_rj <- process_vector_trace(trace_rj)
  
  # Handle true value
  if (!is.null(true_value)) {
    if (length(true_value) > 1) {
      true_value <- mean(true_value)
      warning("true_value length >1 - using mean for plotting")
    }
  }
  
  # Setup graphics only if saving to file
  if (nzchar(output_filename)) {
    pdf(file.path(output_filepath, output_filename), width = 10, height = 5)
    on.exit(dev.off(), add = TRUE)
  }
  
  # Set up plot layout with better legend positioning
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  # Use layout for better control
  layout(matrix(c(1, 2), nrow = 1, ncol = 2))
  par(mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 2, 0))
  
  # 1. Trace Plot with improved legend positioning
  iterations_fixed <- (burn_in + 1):(burn_in + length(trace_fixed))
  iterations_rj <- (burn_in + 1):(burn_in + length(trace_rj))
  
  y_range <- range(c(trace_fixed, trace_rj), na.rm = TRUE)
  
  plot(iterations_fixed, trace_fixed, type = "l", col = colors[1],
       ylim = y_range, xlab = "Iteration", ylab = param,
       main = "Trace Comparison", cex.main = 1.2)
  lines(iterations_rj, trace_rj, col = colors[2])
  grid()
  
  if (burn_in > 0) {
    abline(v = burn_in, col = "gray50", lty = 2, lwd = 1.5)
  }
  
  # Place legend in top-left where there's typically more space
  legend("topleft", 
         legend = c(paste("Fixed (n =", length(trace_fixed)), 
                    paste("RJ (n =", length(trace_rj)), 
                    if(burn_in>0) "Burn-in end"),
         col = c(colors, if(burn_in>0) "gray50"),
         lty = c(1, 1, if(burn_in>0) 2),
         lwd = c(2, 2, if(burn_in>0) 1.5),
         cex = 0.9,
         bg = "NA")
  
  # 2. Density Plot with improved legend positioning
  dens_fixed <- density(trace_fixed)
  dens_rj <- density(trace_rj)
  
  x_range <- range(c(dens_fixed$x, dens_rj$x))
  y_range_dens <- range(c(dens_fixed$y, dens_rj$y))
  
  # Extend y-range slightly for legend space
  y_range_dens[2] <- y_range_dens[2] * 1.1
  
  plot(dens_fixed, col = colors[1], lwd = 2, 
       main = "Density Comparison", cex.main = 1.2,
       xlim = x_range, ylim = y_range_dens, 
       xlab = param, ylab = "Density")
  lines(dens_rj, col = colors[2], lwd = 2)
  grid()
  
  # Add means and true value
  abline(v = mean(trace_fixed), col = colors[1], lty = 2, lwd = 1.5)
  abline(v = mean(trace_rj), col = colors[2], lty = 2, lwd = 1.5)
  
  legend_items <- c("Fixed density", "RJ density", 
                    "Fixed mean", "RJ mean")
  legend_colors <- c(colors, colors)
  legend_lty <- c(1, 1, 1.5, 1.5)
  legend_lwd <- c(2, 2, 1.5, 1.5)
  
  if (!is.null(true_value)) {
    abline(v = true_value, col = "black", lty = 1, lwd = 1.5)
    legend_items <- c(legend_items, "True value")
    legend_colors <- c(legend_colors, "black")
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 3)
  }
  
  # Place legend at top where there's more space
  legend("top", 
         legend = legend_items, 
         col = legend_colors, 
         lty = legend_lty,
         lwd = legend_lwd,
         ncol = 2,
         cex = 0.8,
         bg = NA,         # Transparent background (no color)
         box.lty = 0)     # Remove box border
  
  # Add overall title
  title(paste(param, " comparison: Fixed vs Reversible Jump"), 
        outer = TRUE, cex.main = 1.5)
  
  # Return nothing - just show the plot
  invisible()
}






#' Plot MCMC inferred variables for TWO data sources (R + Python)
#'
#' @param mcmc_r  List with R-side MCMC output (names: rho_trace, ...).
#' @param mcmc_py List with Python-side MCMC output (same structure).
#' @param true_param (optional) list of true values, e.g. list(rho_true = 0.7)
#' @param config     (optional) prior hyper-params list (same as before)
#' @param burn_in    integer burn-in to drop from both traces
#' @param out_file   optional PDF filename; if "", plots go to screen
#' @param out_dir    directory for the PDF
#'
#' @return invisible NULL (makes side-by-side trace+density panels)
#'
plot_mcmc_results_dual <- function(mcmc_r,
                                   mcmc_py,
                                   true_param = list(),
                                   config     = list(),
                                   burn_in    = 0,
                                   out_file   = "",
                                   out_dir    = ".") {
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  suppressPackageStartupMessages({ library(viridis) })
  
  # ----- helper to fetch & burn --------------------------------------------
  get_trace <- function(obj, key) {
    tr <- obj[[key]]
    if (is.null(tr)) return(NULL)
    if (length(tr) > burn_in) tr[-seq_len(burn_in)] else NULL
  }
  
  trace_keys <- list(rho = "rho_trace",
                     prob_noise = "prob_noise_trace",
                     K   = "K_trace")
  
  # colour per source
  pal <- c(R  = "steelblue",
           Py = "darkorange")
  
  # collect for each variable
  plot_list <- list()
  for (v in names(trace_keys)) {
    key <- trace_keys[[v]]
    
    tr_r  <- get_trace(mcmc_r,  key)
    tr_py <- get_trace(mcmc_py, key)
    
    if (is.null(tr_r) && is.null(tr_py)) next
    
    # convert matrices to scalar summary (mean) if needed
    summarise <- function(x) if (is.matrix(x) || is.array(x)) apply(x, 1, mean) else x
    tr_r  <- if (!is.null(tr_r))  summarise(tr_r)  else NULL
    tr_py <- if (!is.null(tr_py)) summarise(tr_py) else NULL
    
    df <- data.frame(
      iter = c(seq_along(tr_r),            seq_along(tr_py)),
      val  = c(tr_r,                       tr_py),
      src  = factor(c(rep("R",  length(tr_r)),
                      rep("Py", length(tr_py))))
    )
    
    plot_list[[v]] <- df
  }
  
  if (length(plot_list) == 0) stop("No trace keys found in either object")
  
  # open pdf if requested
  if (nzchar(out_file)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    pdf(file.path(out_dir, out_file),
        width = 8, height = 4 * length(plot_list))
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(length(plot_list), 2), mar = c(4, 4, 3, 2))
  
  for (v in names(plot_list)) {
    df <- plot_list[[v]]
    
    # --- TRACE --------------------------------------------------------------
    plot(df$iter, df$val, type = "n",
         main = paste("Trace:", v), xlab = "Iter", ylab = v,
         family = "serif")
    grid(col = "grey80", lty = 3)
    by(df, df$src,
       function(sub) lines(sub$iter, sub$val, col = pal[as.character(sub$src)[1]], lwd = 1))
    
    # mark burn-in
    if (burn_in > 0) abline(v = burn_in, col = "red", lty = 2)
    
    legend("topright", legend = names(pal), col = pal, lwd = 2, cex = .8)
    
    # --- DENSITY ------------------------------------------------------------
    hist(df$val[df$src == "R"], breaks = 30, freq = FALSE,
         col = adjustcolor(pal["R"], .4), border = pal["R"],
         xlab = v, ylab = "Density",
         main = paste("Density:", v), family = "serif")
    rug(df$val[df$src == "R"], col = pal["R"])
    lines(density(df$val[df$src == "R"]), col = pal["R"], lwd = 1)
    
    if (any(df$src == "Py")) {
      hist(df$val[df$src == "Py"], breaks = 30, freq = FALSE, add = TRUE,
           col = adjustcolor(pal["Py"], .4), border = pal["Py"])
      rug(df$val[df$src == "Py"], col = pal["Py"])
      lines(density(df$val[df$src == "Py"]), col = pal["Py"], lwd = 1)
    }
    
    # true value?
    true_key <- paste0(v, "_true")
    if (!is.null(true_param[[true_key]]))
      abline(v = true_param[[true_key]], col = "black", lwd = 2, lty = 2)
    
    # sample means
    abline(v = mean(df$val[df$src == "R"]),  col = pal["R"],  lty = 3, lwd = 2)
    if (any(df$src == "Py"))
      abline(v = mean(df$val[df$src == "Py"]), col = pal["Py"], lty = 3, lwd = 2)
  }
  
  if (nzchar(out_file)) dev.off()
  invisible(NULL)
}
