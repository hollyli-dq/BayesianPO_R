# Bayesian Partial Order Inference

A comprehensive R implementation for Bayesian inference of partial orders from preference data using MCMC methods.

## 🚀 Quick Start

``` r
# 1. Install required packages
source("scripts/install_packages.R")

# 2. Test the installation
source("examples/simple_test.R")
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
│   ├── po_inference.R          # Inference functions
│   ├── mcmc.R                  # Fixed dimension MCMC
│   ├── mcmc_rj.R              # Reversible jump MCMC
│   ├── utilities.R            # Utility functions
│   └── analysis.R             # Analysis and plotting
│   └── sushi.R             # Analysis Sushi data 
│
├── config/                     # Configuration files
│   ├── data_generator_config.yaml
│   └── mcmc_config.yaml
│
├── examples/                   # Example scripts and tests
│   ├── simple_test.R          # Basic functionality test
│   ├── simple_inference_test.R # Inference-only test
│   └── test_comprehensive_functions.R # Complete test suite
│
├── scripts/                    # Utility scripts
│   └── install_packages.R     # Package installation script
│
├── docs/                       # Documentation
│   ├── README_COMPREHENSIVE_FUNCTIONS.md # Function documentation
│   ├── QUICK_SETUP.md         # Quick setup guide
│   ├── MCMC_Simulation_Tutorial.Rmd # Detailed tutorial
│   └── sushi_study.rmd         # The sushi study 
├── data/                       # Data directory (for your datasets, including the sushi)
├── results/                    # Output directory for results
├── tests/                      # Unit tests (R package structure)
├── man/                        # Manual pages (R package structure)
└── inst/                       # Installed files (R package structure)
```

## 🎯 Main Features

-   **Data Generation**: Synthetic partial order data with noise models
-   **Fixed Dimension MCMC**: Standard MCMC for known dimensionality
-   **Reversible Jump MCMC**: Automatic dimension selection
-   **Command Line Interface**: Easy-to-use CLI for batch processing
-   **Comprehensive Analysis**: Plotting and comparison tools

## 📖 Documentation

| Document | Description |
|--------------------------------|----------------------------------------|
| [`docs/QUICK_SETUP.md`](docs/QUICK_SETUP.md) | 3-step setup guide |
| [`docs/README_COMPREHENSIVE_FUNCTIONS.md`](docs/README_COMPREHENSIVE_FUNCTIONS.md) | Complete function documentation |
| [`docs/MCMC_Simulation_Tutorial.Rmd`](docs/MCMC_Simulation_Tutorial.Rmd) | Detailed tutorial |

## 🔧 Installation

### Prerequisites

-   R version 4.0.0 or higher
-   Required packages: `yaml`, `jsonlite`, `mvtnorm`, `optparse`, `igraph`

### Automated Installation

``` r
source("scripts/install_packages.R")
```

### Manual Installation

``` r
install.packages(c("yaml", "jsonlite", "mvtnorm", "optparse", "igraph"))
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

## 🧪 Examples

| Example | Description | Usage |
|----------------------|--------------------------------|------------------|
| [`examples/simple_test.R`](examples/simple_test.R) | Basic data generation test | `source("examples/simple_test.R")` |
| [`examples/simple_cli_test.R`](examples/simple_cli_test.R) | Full pipeline with both MCMC methods | `source("examples/simple_cli_test.R")` |
| [`examples/simple_inference_test.R`](examples/simple_inference_test.R) | Inference-only example | `source("examples/simple_inference_test.R")` |
| [`examples/test_comprehensive_functions.R`](examples/test_comprehensive_functions.R) | Complete test suite | `source("examples/test_comprehensive_functions.R")` |

## 🔬 Methods

### Noise Models

-   **Queue Jump**: Models preference reversals as adjacent swaps
-   **Mallows**: Distance-based noise model for rankings

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
3.  **Test installation**: `source("examples/simple_test.R")`
4.  **Run full pipeline**: `source("examples/simple_cli_test.R")`
5.  **Explore documentation**: Check `docs/` folder
6.  **Customize configs**: Edit files in `config/`
7.  **Run your analysis**: Use CLI or R functions

## 🔍 Troubleshooting

### Getting Help

-   Check [`docs/QUICK_SETUP.md`](docs/QUICK_SETUP.md) for step-by-step guide
-   Review [`docs/README_COMPREHENSIVE_FUNCTIONS.md`](docs/README_COMPREHENSIVE_FUNCTIONS.md) for detailed documentation
-   Run `Rscript R/cli.R --help` for CLI options

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📧 Contact

For questions or issues, please open a GitHub issue or contact the maintainers.

## 
