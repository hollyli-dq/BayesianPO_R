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
