# Test Script for Bayesian Partial Order Inference in R
# This script tests the main functionality of the converted R package

context("Basic Functionality Tests")

test_that("Package loads correctly", {
  expect_true(require(BayesianPO, quietly = TRUE))
})

test_that("Synthetic data generation works", {
  test_data <- generate_synthetic_data(
    n_items = 4,
    n_observations = 20,
    K = 2,
    rho_true = 0.7,
    prob_noise_true = 0.1,
    random_seed = 42
  )
  
  expect_equal(length(test_data$items), 4)
  expect_equal(length(test_data$observed_orders), 20)
  expect_equal(test_data$rho_true, 0.7)
  expect_equal(test_data$prob_noise_true, 0.1)
  expect_true(is.list(test_data$observed_orders))
  expect_true(is.list(test_data$choice_sets))
  expect_true(is.matrix(test_data$X))
})

test_that("MCMC inference runs without errors", {
  test_data <- generate_synthetic_data(
    n_items = 4,
    n_observations = 20,
    K = 2,
    rho_true = 0.7,
    prob_noise_true = 0.1,
    random_seed = 42
  )
  
  expect_silent({
    results <- mcmc_partial_order(
      observed_orders = test_data$observed_orders,
      choice_sets = test_data$choice_sets,
      num_iterations = 100,  # Very short run for testing
      K = 2,
      X = test_data$X,
      noise_option = "queue_jump",
      random_seed = 42
    )
  })
  
  expect_true(is.list(results))
  expect_true("final_rho" %in% names(results))
  expect_true("final_prob_noise" %in% names(results))
  expect_true("overall_acceptance_rate" %in% names(results))
})

test_that("Complete analysis pipeline works", {
  test_data <- generate_synthetic_data(
    n_items = 4,
    n_observations = 20,
    K = 2,
    rho_true = 0.7,
    prob_noise_true = 0.1,
    random_seed = 42
  )
  
  expect_silent({
    analysis_results <- run_po_analysis(
      observed_orders = test_data$observed_orders,
      choice_sets = test_data$choice_sets,
      X = test_data$X,
      num_iterations = 200,
      K = 2,
      burn_in = 50,
      random_seed = 42
    )
  })
  
  expect_true(is.list(analysis_results))
  expect_true("summary_stats" %in% names(analysis_results))
  expect_true("h_final" %in% names(analysis_results))
  expect_true("mcmc_results" %in% names(analysis_results))
})

test_that("Utility functions work correctly", {
  # Test correlation matrix
  sigma <- build_sigma_rho(3, 0.5)
  expect_equal(dim(sigma), c(3, 3))
  expect_equal(diag(sigma), c(1, 1, 1))
  expect_equal(sigma[1, 2], 0.5)
  
  # Test transitive closure/reduction
  test_matrix <- matrix(c(0,1,0,0,0,1,0,0,0), nrow=3)
  tc <- transitive_closure(test_matrix)
  tr <- transitive_reduction(tc)
  expect_true(is.matrix(tc))
  expect_true(is.matrix(tr))
  
  # Test transformation
  U_test <- matrix(rnorm(6), nrow=2, ncol=3)
  alpha_test <- c(0.1, 0.2)
  eta_test <- transform_U_to_eta(U_test, alpha_test)
  expect_equal(dim(eta_test), c(2, 3))
})

test_that("Different noise models work", {
  test_data <- generate_synthetic_data(
    n_items = 4,
    n_observations = 10,
    K = 2,
    rho_true = 0.7,
    prob_noise_true = 0.1,
    random_seed = 42
  )
  
  # Test Mallows model
  expect_silent({
    mallows_results <- mcmc_partial_order(
      observed_orders = test_data$observed_orders,
      choice_sets = test_data$choice_sets,
      num_iterations = 50,
      K = 2,
      X = test_data$X,
      noise_option = "mallows_noise",
      random_seed = 42
    )
  })
  
  expect_true(is.list(mallows_results))
})

test_that("Data format validation works", {
  test_data <- generate_synthetic_data(
    n_items = 4,
    n_observations = 10,
    K = 2,
    random_seed = 42
  )
  
  expect_true(is.list(test_data$observed_orders))
  expect_true(is.list(test_data$choice_sets))
  expect_true(is.matrix(test_data$X))
  expect_equal(nrow(test_data$X), length(test_data$items))
})

test_that("Example analysis runs", {
  # This is a longer test, so we'll make it conditional
  skip_on_cran()
  
  expect_silent({
    example_results <- run_example_analysis()
  })
  
  expect_true(is.list(example_results))
  expect_true("data" %in% names(example_results))
  expect_true("results" %in% names(example_results))
}) 