# MCMC Scaling Analysis Script

This directory contains a script to evaluate the computational performance and timing characteristics of the BayesianPO MCMC algorithms across different problem sizes (number of nodes/items in the partial order).

## Overview

The scaling analysis evaluates:
- **Data Generation Time**: How long it takes to generate synthetic partial order data
- **MCMC Runtime**: Computational time for both Fixed-K and Reversible Jump MCMC
- **Acceptance Rates**: MCMC performance indicators
- **Scaling Behavior**: How computational time changes with problem size

## Scripts

### 1. `scaling_analysis.R` - Timing Analysis
**Purpose**: Computational scaling analysis across node sizes 5, 6, 7, 8, 9, 10, 15, 20, 25

**Features**:
- Runs both Fixed-K MCMC and Reversible Jump MCMC
- 50,000 MCMC iterations per algorithm
- Timing measurements for data generation and MCMC execution
- Acceptance rate monitoring
- Generates timing analysis plots
- Creates comprehensive CSV results

**Runtime**: ~30 minutes to several hours depending on your machine

**Usage**:
```bash
cd scripts
Rscript scaling_analysis.R
```

**Outputs**:
- `scaling_analysis_results.csv` - Detailed timing results table
- `scaling_timing_analysis.pdf` - Timing performance plots
- `scaling_analysis_intermediate.csv` - Intermediate results (auto-saved)

### 2. `view_results.R` - Results Viewer
**Purpose**: Display results in a formatted, easy-to-read table

**Features**:
- Automatically detects which results file to load
- Formats timing information nicely
- Provides performance analysis and complexity estimates
- Works with both full analysis and quick test results

**Usage**:
```bash
cd scripts
Rscript view_results.R
```

## Results Table Columns

### Scaling Analysis (`scaling_analysis_results.csv`)
| Column | Description |
|--------|-------------|
| `n_nodes` | Number of items/nodes in partial order |
| `data_gen_time` | Data generation time (seconds) |
| `fixed_mcmc_time` | Fixed-K MCMC runtime (seconds) |
| `rj_mcmc_time` | Reversible Jump MCMC runtime (seconds) |
| `total_time` | Total runtime (data generation + both MCMC runs) |
| `fixed_acceptance_rate` | Fixed MCMC acceptance rate |
| `rj_acceptance_rate` | RJ-MCMC acceptance rate |
| `rj_mean_K` | RJ-MCMC estimated number of dimensions |

## Configuration

### Key Parameters (editable in scripts)

**Problem Size**:
- `node_sizes`: Vector of node counts to test
- `N`: Number of observed orderings per problem
- `K`: True number of latent dimensions

**MCMC Settings**:
- `iters`: Number of MCMC iterations
- `burn_in`: Burn-in period
- `mcmc_pt`: Update probabilities for different move types

**Data Generation**:
- `rho_true`: True correlation parameter
- `prob_noise_true`: Queue-jump noise level
- `beta_true`: True regression coefficients

### Customization Examples

**Test larger problems**:
```r
node_sizes <- c(10, 15, 20, 25, 30, 40, 50)
```

**Longer MCMC for better accuracy**:
```r
mcmc_config$iters <- 100000
burn_in <- 10000
```

**Different noise levels**:
```r
prob_noise_true <- 0.05  # Less noise
# or
prob_noise_true <- 0.20  # More noise
```

## Interpreting Results

### Timing Analysis
- **Linear scaling**: Time increases proportionally with nodes O(n)
- **Quadratic scaling**: Time increases with nodes² O(n²)
- **Cubic scaling**: Time increases with nodes³ O(n³)
- **Exponential scaling**: Time increases exponentially (problematic)

### Performance Metrics
- **Acceptance Rate**: MCMC proposal acceptance rate (optimal: 20-50%)
- **K Estimation**: RJ-MCMC's estimated number of latent dimensions
- **Total Runtime**: Combined time for data generation and both MCMC algorithms

### Expected Patterns
- Data generation time should scale roughly O(n²) due to partial order operations
- MCMC time typically scales O(n²) to O(n³) depending on implementation
- Acceptance rates should remain relatively stable across problem sizes
- RJ-MCMC should consistently estimate K ≈ 3 (the true value)

## Troubleshooting

### Memory Issues
If you encounter memory problems with large node sizes:
1. Reduce the number of MCMC iterations
2. Reduce the number of observed orderings (`N`)
3. Test smaller node sizes first

### Long Runtime
For the full analysis:
1. Start with the quick test to verify everything works
2. Consider running overnight for the full analysis
3. Monitor intermediate results files

### Accuracy Issues
If accuracy is unexpectedly low:
1. Check that burn-in period is sufficient
2. Verify MCMC acceptance rates (should be ~20-50%)
3. Increase number of iterations
4. Check for numerical issues in the log files

## Expected Runtime Estimates

| Script | Node Sizes | Approximate Runtime |
|--------|------------|-------------------|
| Quick Test | 5,6,7,8 | 5-15 minutes |
| Full Analysis | 5-10 | 30-60 minutes |
| Full Analysis | 5-25 | 2-6 hours |

*Times vary significantly based on hardware and exact configuration*

## Files Generated

- `*.csv`: Results tables (importable into Excel, Python, etc.)
- `*.pdf`: Plots and visualizations
- `*_intermediate.csv`: Partial results (for long-running analyses)

## Next Steps

After running the analysis:
1. Use `view_results.R` to get a formatted summary
2. Import CSV files into your preferred analysis software
3. Examine the PDF plots for visual insights
4. Consider parameter tuning based on acceptance rates and accuracy
5. Scale up to larger problems if performance is acceptable 