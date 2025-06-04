library(mvtnorm)
source("examples/simple_test.R")
source("R/utilities.R")
source("R/mcmc.R")
source("R/mcmc_rj.R")  # Load reversible jump MCMC
source("R/analysis.R")


set.seed(42)
K      <- 3
rho    <- 0.4
n      <- 5
min_sub <-3 
Sigma   <- build_sigma_rho(K, rho)
Z       <- rmvnorm(n, sigma = Sigma)
alpha   <- rnorm(n)
eta     <- transform_U_to_eta(Z, alpha)
h       <- generate_partial_order(eta)

# ---- write CSV fixtures ----------------------------------------------------
if (!dir.exists("fixtures")) dir.create("fixtures")
write.table(Sigma,  "fixtures/Sigma.csv",  sep = ",", row.names = FALSE, col.names = FALSE)
write.table(Z,      "fixtures/Z.csv",      sep = ",", row.names = FALSE, col.names = FALSE)
write.table(alpha,  "fixtures/alpha.csv",  sep = ",", row.names = FALSE, col.names = FALSE)
write.table(eta,    "fixtures/eta.csv",    sep = ",", row.names = FALSE, col.names = FALSE)
write.table(h,      "fixtures/h.csv",      sep = ",", row.names = FALSE, col.names = FALSE)


# derived deterministic results ---------------------------------------------
tc <- transitive_closure(h)
tr <- transitive_reduction(h)
le_count <- count_linear_extensions(h)

write.table(tc, "fixtures/tc.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(tr, "fixtures/tr.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write("\n"  , file = "fixtures/le_count.txt") # ensure file exists even if write below fails
write(le_count, file = "fixtures/le_count.txt")

#BELOW IS FOR THE MCMC 

source("R/mcmc.R")        # file that contains log_likelihood_*()
K           <- 2
rho_true    <- 0.3
n_items     <- 4
Sigma_true  <- build_sigma_rho(K, rho_true)
Z_true      <- rmvnorm(n_items, sigma = Sigma_true)
alpha_true  <- rnorm(n_items)
eta_true    <- transform_U_to_eta(Z_true, alpha_true)
h_true      <- generate_partial_order(eta_true)

# make two observed orders on the full item set
items <- paste0("I", 1:n_items)
observed_orders <- list(items, rev(items))
choice_sets     <- list(items, items)

# mapping
item_to_index   <- setNames(seq_along(items), items)

# convert observed orders to index form for Python
obs_idx <- lapply(observed_orders, \(ord) unname(item_to_index[ord]))

# save fixtures ---------------------------------------------------
dir.create("llfixture", showWarnings = FALSE)
write.table(h_true,      "ll_fixture/h.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)
# ---------- save list-of-lists as rectangular CSV -----------------
max_len <- max(vapply(obs_idx, length, 0))
pad_to  <- function(x, n) { length(x) <- n; x }        # pad with NA
obs_mat <- do.call(rbind, lapply(obs_idx, pad_to, max_len))
write.table(obs_mat,
            "llfixture/obs_idx.csv",
            sep = ",", row.names = FALSE, col.names = FALSE, na = "")

max_len_cs <- max(vapply(choice_sets, length, 0))
cs_mat <- do.call(rbind, lapply(choice_sets, pad_to, max_len_cs))
write.table(cs_mat,
            "llfixture/choice_sets.csv",
            sep = ",", row.names = FALSE, col.names = FALSE, na = "")

write.table(item_to_index,
            "llfixture/item_to_index.csv",
            sep = ",", row.names = TRUE, col.names = FALSE)

# A single likelihood value from R (queue-jump)
ll_R <- log_likelihood_queue_jump(
  h_true, obs_idx, choice_sets, item_to_index,
  prob_noise = 0.15)
cat("R   log-likelihood:", ll_R, "\n")





n                <- 5          # number of items
N                <- 10         # number of total orders
K                <- 3           # latent dimensions
p                <- 2           # number of covariates
min_sub<- 3 
rho_prior        <- 0.16667     # prior mean for rho  (Beta(1, (1/ρ₀)–1))
prob_noise_true  <- 0.10        # queue-jump noise used to *simulate* data
beta_true        <- c(0.5, -0.3)# true regression coefficients
rng_seed         <- 42          # reproducibility
noise_option<-  "queue_jump"
set.seed(rng_seed)

# draw rho_true from Beta(1, b) with mean = rho_prior
rho_true <- rbeta(1, 1, (1 / rho_prior) - 1)

# design matrix  X  of size p × n  (independent standard-normal entries)
X <- matrix(rnorm(p * n), nrow = p, ncol = n)

synthetic_data <- generate_synthetic_data(
  n_items         = n,
  n_observations  = N,
  min_sub.        = min_sub,
  K               = K,
  rho_true        = rho_true,
  prob_noise_true = prob_noise_true,
  beta_true       = beta_true,
  X               = t(X),          # generator expects n × p
  random_seed     = rng_seed
)


# Extract components
items <- synthetic_data$items
observed_orders <- synthetic_data$observed_orders
choice_sets <- synthetic_data$choice_sets
h_true <- synthetic_data$h_true

plot_partial_order(
  h_true, items,
  title         = "True Partial Order",
  vertex_size   = 38,
  edge_arrow_size = 0.6,
  label_cex     = 1.3,  # even bigger text
  frame_width   = 2.3   # bolder outlines
)
