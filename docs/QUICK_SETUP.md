# Quick Setup Guide - Bayesian Partial Order Inference

## 🚀 Get Started in 3 Steps

### Step 1: Install Required Packages
```r
# Run the automated installation script
source("install_packages.R")
```

This will:
- ✅ Check for missing packages
- ✅ Install from multiple CRAN mirrors
- ✅ Test package loading
- ✅ Display system information

### Step 2: Test the Installation
```r
# Quick test of data generation
source("simple_test.R")
```

Expected output:
```
Testing data generation...
True correlation parameter (rho): 0.XXXX
Data generation successful!
Generated 20 observed orders
True partial order has X edges
Data saved to output/simple_test_data.json
```

### Step 3: Run Full Pipeline
```r
# Test both MCMC methods
source("simple_cli_test.R")
```

Expected output:
```
=== Testing Comprehensive Functions ===

1. Data Generation
------------------
Data generation completed!

2. Fixed Dimension MCMC
-----------------------
Fixed MCMC completed!
Inferred edges: X
Acceptance rate: XX.X %

3. Reversible Jump MCMC
-----------------------
RJ-MCMC completed!
Inferred edges: X
Acceptance rate: XX.X %
Final K: X

=== All Tests Completed Successfully ===
```

## 🎯 Command Line Usage

### Basic Commands
```bash
# Show help
Rscript R/cli.R --help

# Generate data only
Rscript R/cli.R --generate-data

# Run full pipeline with custom settings
Rscript R/cli.R --iterations 1000 --use-rj-mcmc --verbose
```

### Example Workflow
```bash
# 1. Generate synthetic data
Rscript R/cli.R --generate-data --output-dir my_results

# 2. Run fixed dimension MCMC
Rscript R/cli.R --inference-only --iterations 2000 --burn-in 500

# 3. Run reversible jump MCMC
Rscript R/cli.R --inference-only --use-rj-mcmc --iterations 2000
```

## 📁 File Structure

After setup, you'll have:
```
BayesianPO/
├── R/
│   ├── cli.R                    # Command line interface
│   ├── data_generator.R         # Data generation functions
│   ├── po_inference.R          # Inference functions
│   ├── mcmc.R                  # Fixed dimension MCMC
│   ├── mcmc_rj.R               # Reversible jump MCMC
│   ├── utilities.R             # Utility functions
│   └── analysis.R              # Analysis and plotting
├── config/
│   ├── data_generator_config.yaml
│   └── mcmc_config.yaml
├── output/                     # Generated results
├── install_packages.R          # Package installer
└── test files...
```

## ⚙️ Configuration

### Quick Config Changes

**Change dataset size:**
```yaml
# config/data_generator_config.yaml
generation:
  n: 8        # Number of items (was 5)
  N: 100      # Number of observations (was 50)
```

**Change MCMC settings:**
```yaml
# config/mcmc_config.yaml
mcmc:
  num_iterations: 2000    # More iterations (was 1000)
  
update_probabilities:
  K: 0.1                  # Enable RJ-MCMC (was 0.0)
```

## 🔧 Troubleshooting

### Package Installation Issues
```r
# If installation fails, try:
install.packages("package_name", repos = "https://cloud.r-project.org/")

# For igraph issues on Linux/Mac:
# Install system dependencies first
# Ubuntu: sudo apt-get install libxml2-dev libcurl4-openssl-dev
# Mac: brew install libxml2
```

### Memory Issues
```r
# Increase memory limit (Windows)
memory.limit(size = 8000)

# Or reduce dataset size in configs
```

### YAML Parsing Errors
- Ensure config files end with a newline
- Check for proper indentation (spaces, not tabs)
- Verify quotes around string values

## 📊 Expected Performance

| Operation | Time | Memory | Output |
|-----------|------|--------|--------|
| Data Generation | < 1s | ~10MB | JSON file |
| Fixed MCMC (1000 iter) | ~20s | ~50MB | 30-40% acceptance |
| RJ-MCMC (1000 iter) | ~30s | ~75MB | 50-60% acceptance |

## 🎓 Next Steps

1. **Modify configs** for your specific use case
2. **Run longer chains** for better convergence
3. **Analyze results** using the plotting functions
4. **Compare methods** (fixed vs. reversible jump)
5. **Scale up** to larger datasets

## 📚 Full Documentation

- `README_COMPREHENSIVE_FUNCTIONS.md` - Complete function documentation
- `COMPREHENSIVE_FUNCTIONS_SUMMARY.md` - Implementation summary
- Individual R files have detailed function documentation

## ✅ Success Indicators

You're ready to go when:
- ✅ All packages install without errors
- ✅ `simple_test.R` generates data successfully  
- ✅ `simple_cli_test.R` runs both MCMC methods
- ✅ Output files appear in `output/` directory
- ✅ Acceptance rates are reasonable (>20%)

Happy inferencing! 🎉 