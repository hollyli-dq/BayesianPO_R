library(mvtnorm)

source("R/utilities.R")

set.seed(42)
K      <- 3
rho    <- 0.4
n      <- 5

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
