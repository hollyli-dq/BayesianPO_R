### Below is the test for data generator 

source("scripts/install_packages.R")
# Load the BayesianPO package
source("R/utilities.R")
source("R/mcmc.R")
source("R/mcmc_rj.R")  # Load reversible jump MCMC
source("R/analysis.R")

library(mvtnorm)
library(ggplot2)
library(igraph)
library(graph)

# Set seed for reproducibility
set.seed(42)

# 1. Set up parameters for partial order generation
n <- 5  # Number of nodes/items in the partial order
N <- 30  # Number of total orders to sample from the partial order

# 2. Load configuration parameters (example config - replace with your actual config)
config <- list(
  mcmc = list(K = 3),
  prior = list(rho_prior = 0.16667, mallow_ua = 1.0,noise_beta_prior =9.0),
  noise = list(noise_option = "queue_jump"),
  covariates = list(p = 2, beta_true = c(0.5, -0.3))
)

K <- config$mcmc$K  # Number of dimensions for latent positions
rho_prior_val <- config$prior$rho_prior  # Prior parameter for correlation
noise_option <- config$noise$noise_option  # Noise model specification
mallow_ua <- config$prior$mallow_ua  # Mallows model parameter

items <- 1:n  # Create list of item indices

# 3. Generate true correlation parameter from prior
rho_true <- rbeta(1, 1, rho_prior_val)  # Sample from Beta(1, rho_prior)
cat(sprintf("True correlation parameter (rho): %.4f\n", rho_true))

# 4. Set up covariates and regression parameters
p <- config$covariates$p  # Number of covariates
beta_true <- config$covariates$beta_true   # True regression coefficients
X <- matrix(rnorm(p * n), nrow = p, ncol = n)  # p x n design matrix

# 5. Compute assessor-specific effects (α = Xᵀβ)
alpha <- as.vector(t(X) %*% beta_true)  # n-dimensional vector
cat("\nThe covariates effects (alpha):\n")
print(alpha)

# 6. Generate latent positions for each item
Sigma <- build_sigma_rho(K, rho_true)  # Correlation matrix
U <- mvtnorm::rmvnorm(n, mean = rep(0, K), sigma = Sigma)  # n x K matrix
cat("\nBase U matrix (latent positions):\n")
print(U)

# 7. Generate latent positions with covariates effects (η = g(U) + α)
eta <- transform_U_to_eta(U, alpha)
cat("\nAdjusted latent positions (eta):\n")
print(eta)

# 8. Generate partial order from adjusted latent positions
h_full <- generate_partial_order(eta)  # Full partial order
h_true <- transitive_reduction(h_full)  # Transitive reduction
cat("\nPartial Order (adjacency matrix):\n")
print(h_true)

# 9. Print descriptive statistics of the generated partial order
cat("\nPartial Order Statistics:\n")
cat(sprintf("Number of items: %d\n", n))
cat(sprintf("Number of covariates: %d\n", p))
cat(sprintf("Number of direct relationships: %d\n", sum(h_true)))
cat(sprintf("Number of linear extensions: %d\n", nle(h_true)))

plot_partial_order(
  h_true, items,
  title         = "True Partial Order",
  vertex_size   = 38,
  edge_arrow_size = 0.6,
  label_cex     = 1.3,  # even bigger text
  frame_width   = 2.3   # bolder outlines
)

generate_subsets <- function(N, n, min_sub= 2, max_size = NULL) {
  if (is.null(max_size)) max_size <- n
  subsets <- list()
  
  for (i in 1:N) {
    size <- sample(min_sub:max_size, 1)
    subset <- sample(1:n, size, replace = FALSE)
    subsets[[i]] <- sort(subset)  # Return sorted for consistency
  }
  return(subsets)
}

subsets <- generate_subsets(N, n)
cat("\nGenerated subsets (first 5 shown):\n")
print(head(subsets, 5))

h_tc <- transitive_closure(h_true)

# 13. Generate noise probability from prior
# -----------------------------------------------------------------------------
noise_beta_prior <- config$prior$noise_beta_prior
prob_noise <- rbeta(1, 1, noise_beta_prior)
prob_noise_true <- prob_noise
cat(sprintf("\nGenerated noise probability: %.4f\n", prob_noise_true))


total_orders <- list()

for (choice_set in subsets) {
  if (noise_option == "queue_jump") {
    y_generated <- generate_total_order_queue_jump(
      subset = choice_set,
      items_all = items,
      h_global = h_tc,
      prob_noise = prob_noise
    )
  } else if (noise_option == "mallows_noise") {
    # Mallows model not implemented yet - placeholder
    stop("Mallows noise model not implemented in this translation")
  }
  total_orders <- c(total_orders, list(y_generated))
}

# 15. Analyze the order distribution
# -----------------------------------------------------------------------------
# Convert orders to string representation for easier comparison
order_strings <- sapply(total_orders, function(x) paste(x, collapse = " > "))

# Count unique orders
unique_orders <- unique(order_strings)
num_unique <- length(unique_orders)

# Find most common ordering
order_counts <- table(order_strings)
most_common <- names(which.max(order_counts))

cat("\nSampling Statistics:\n")
cat(sprintf("Number of total orders generated: %d\n", length(total_orders)))
cat(sprintf("Number of unique total orders: %d\n", num_unique))
cat(sprintf("Most common ordering: %s\n", most_common))

# 16. Visualize some example orders (optional)
# -----------------------------------------------------------------------------
if (length(total_orders) > 0) {
  cat("\nExample generated orders (first 3):\n")
  for (i in 1:min(3, length(total_orders))) {
    cat(sprintf("Order %d (Set size %d): %s\n", 
                i, 
                length(total_orders[[i]]),
                paste(total_orders[[i]], collapse = " > ")))
  }
}


# Save results to JSON
# -----------------------------------------------------------------------------
library(jsonlite)

# Create data structure matching Python format
data <- list(
  observed_orders = total_orders,
  choice_sets = subsets,
  items = items,
  parameters = list(
    n = n,
    N = N,
    K = K,
    rho_true = round(rho_true, 4),
    prob_noise_true = round(prob_noise_true, 4)
  ),
  true_partial_order = h_tc,
  beta_true = beta_true,
  X = t(X),  # Transpose to match Python's orientation (n x p)
  U_true = U,
  alpha_true = alpha,
  eta_true = eta
)

# Custom conversion function for JSON serialization
convert_for_json <- function(obj) {
  if (is.matrix(obj)) {
    # Convert matrix to list of lists (row-wise)
    lapply(1:nrow(obj), function(i) as.list(obj[i, ]))
  } else if (is.list(obj)) {
    # Recursively process list elements
    lapply(obj, convert_for_json)
  } else if (is.vector(obj) && length(obj) > 1) {
    # Convert vectors to lists
    as.list(obj)
  } else {
    # Return scalars as-is
    obj
  }
}

# Apply conversion to all elements
json_data <- convert_for_json(data)

# Save to JSON file
write_json(
  json_data, 
  "partial_order_data.json",
  pretty = TRUE,
  auto_unbox = TRUE,  # Remove unnecessary brackets for scalars
  digits = 6          # Control floating point precision
)

cat("\nData saved to partial_order_data.json\n")
