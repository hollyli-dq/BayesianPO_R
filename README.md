# Bayesian Partial Order Inference

A comprehensive R implementation for Bayesian inference of partial orders from preference data using MCMC methods.

## 🚀 Quick Start

``` r
# 1. Install required packages
source("scripts/install_packages.R")
```

## 📁 Project Structure

```         
BayesianPO/
├── README.md                    # This file
├── DESCRIPTION                  # R package description
├── NAMESPACE                    # R package namespace
│
├── R/                          # Core R functions
│   ├── data_generator.R        # Data generation functions
│   ├── mcmc.R                  # Fixed dimension MCMC
│   ├── mcmc_rj.R              # Reversible jump MCMC
│   ├── utilities.R            # Utility functions
│   └── analysis.R             # Analysis and plotting
│
├── config/                     # Configuration files
│   ├── data_generator_config.yaml
│   └── mcmc_config.yaml
│
├── scripts/                    # Utility scripts
│   └── install_packages.R     # Package installation script
│
├── docs/                       # Documentation
│   ├── MCMC_Simulation_Tutorial.Rmd # Detailed tutorial
├── data/                       # Data directory (for your datasets, including the sushi)
├── results/                    # Output directory for results
```

## 🎯 Main Features

-   **Data Generation**: Synthetic partial order data with noise models
-   **Fixed Dimension MCMC**: Standard MCMC for known dimensionality
-   **Reversible Jump MCMC**: Automatic dimension selection
-   **Comprehensive Analysis**: Plotting and comparison tools

## 📖 Documentation

| Document | Description |
|----|----|
| [`docs/MCMC_Simulation_Tutorial.Rmd`](docs/MCMC_Simulation_Tutorial.Rmd) | Detailed tutorial |

## 🔧 Installation

### Prerequisites

-   R version 4.0.0 or higher
-   Required packages: `yaml`, `jsonlite`, `mvtnorm`, `optparse`, `igraph`

### Automated Installation

``` r
source("scripts/install_packages.R")
```

## 💻 Usage

### R Functions

``` r
# Load functions
source("R/data_generator.R")
source("R/po_inference.R")

# Generate data
config <- load_config("config/data_generator_config.yaml")
data <- generate_data(config)

# Run inference
mcmc_config <- load_config("config/mcmc_config.yaml")
results <- run_inference(data, mcmc_config, use_rj_mcmc = TRUE)

# Save and analyze results
save_results(results, "results", "my_analysis")
generate_plots(results, data, mcmc_config, "results", "my_analysis")
```

## 🔬 Methods

### Noise Models

-   **Queue Jump**: Models preference reversals as adjacent swaps

### MCMC Algorithms

-   **Fixed Dimension**: Standard MCMC with known latent dimensions
-   **Reversible Jump**: Automatic model selection with birth/death moves

### Key Features

-   Proper posterior averaging with thresholding
-   Transitive reduction for minimal partial orders
-   Comprehensive diagnostics and plotting
-   JSON-based data interchange

## ⚙️ Configuration

### Data Generation Config (`config/data_generator_config.yaml`)

``` yaml
generation:
  n: 5                    # Number of items
  N: 50                   # Number of observations
  prob_noise_true: 0.1    # Noise probability

covariates:
  p: 2                    # Number of covariates
  beta_true: [0.5, -0.3]  # True coefficients
```

### MCMC Config (`config/mcmc_config.yaml`)

``` yaml
mcmc:
  num_iterations: 1000    # MCMC iterations
  random_seed: 42         # Random seed
  
update_probabilities:
  rho: 0.2               # Update probabilities
  noise: 0.2
  U: 0.4
  beta: 0.2
  K: 0.0                 # Set to 0.1+ for RJ-MCMC
```

## 🚀 Getting Started

1.  **Clone/Download** the repository
2.  **Install packages**: `source("scripts/install_packages.R")`
3.  **Explore documentation**: Check `docs/` folder
4.  **Customize configs**: Edit files in `config/`

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📧 Contact

For questions or issues, please open a GitHub issue or contact the maintainers.

## 
