# Package Installation Script for Bayesian Partial Order Inference
# This script installs all required R packages for the comprehensive functions

cat("=== Bayesian Partial Order Inference - Package Installation ===\n\n")

# Define required packages
required_packages <- c(
  "yaml",        # Configuration file parsing
  "jsonlite",    # JSON data handling and serialization
  "mvtnorm",     # Multivariate normal distributions
  "optparse",    # Command-line argument parsing
  "igraph"       # Graph visualization and analysis
)

# Optional packages (usually pre-installed with R)
optional_packages <- c(
  "tools",       # File path utilities
  "stats",       # Statistical functions
  "utils",       # Utility functions
  "grDevices",   # Graphics devices
  "graphics"     # Base graphics
)

cat("Required packages:\n")
for(pkg in required_packages) {
  cat(" -", pkg, "\n")
}

cat("\nChecking package installation status...\n")

# Check which packages are not installed
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) > 0) {
  cat("\nMissing packages found:", paste(missing_packages, collapse = ", "), "\n")
  cat("Installing missing packages...\n")
  
  # Try multiple repositories in case one fails
  repos <- c(
    "https://cran.rstudio.com/",
    "https://cloud.r-project.org/",
    "https://cran.r-project.org/"
  )
  
  for(pkg in missing_packages) {
    installed <- FALSE
    for(repo in repos) {
      tryCatch({
        install.packages(pkg, repos = repo, dependencies = TRUE)
        cat("✓ Successfully installed:", pkg, "\n")
        installed <- TRUE
        break
      }, error = function(e) {
        cat("✗ Failed to install", pkg, "from", repo, "\n")
      })
    }
    
    if(!installed) {
      cat("ERROR: Could not install package:", pkg, "\n")
      cat("Please try installing manually with: install.packages('", pkg, "')\n", sep = "")
    }
  }
} else {
  cat("✓ All required packages are already installed!\n")
}

cat("\nTesting package loading...\n")

# Test loading all packages
all_loaded <- TRUE
for(pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat("✓", pkg, "loaded successfully\n")
  }, error = function(e) {
    cat("✗", pkg, "failed to load:", e$message, "\n")
    all_loaded <- FALSE
  })
}

# Check optional packages
cat("\nChecking optional packages...\n")
for(pkg in optional_packages) {
  if(pkg %in% installed.packages()[,"Package"]) {
    tryCatch({
      library(pkg, character.only = TRUE, quietly = TRUE)
      cat("✓", pkg, "available\n")
    }, error = function(e) {
      cat("⚠", pkg, "installed but failed to load\n")
    })
  } else {
    cat("⚠", pkg, "not installed (usually comes with base R)\n")
  }
}

# Final status
cat("\n=== Installation Summary ===\n")
if(all_loaded) {
  cat("✅ SUCCESS: All required packages are installed and working!\n")
  cat("\nYou can now run the comprehensive functions:\n")
  cat("- source('examples/simple_test.R')           # Test data generation\n")
  cat("- source('examples/simple_cli_test.R')       # Test full pipeline\n")
  cat("- Rscript R/cli.R --help                    # Command line interface\n")
} else {
  cat("❌ ISSUES: Some packages failed to load properly.\n")
  cat("\nTroubleshooting steps:\n")
  cat("1. Check your R version: R.version.string\n")
  cat("2. Update R if version < 4.0.0\n")
  cat("3. Try installing packages manually\n")
  cat("4. Check for system dependencies (especially for igraph)\n")
}

# Display system information
cat("\n=== System Information ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")

# Package versions
cat("\n=== Package Versions ===\n")
for(pkg in required_packages) {
  if(pkg %in% installed.packages()[,"Package"]) {
    version <- packageVersion(pkg)
    cat(pkg, ":", as.character(version), "\n")
  }
}

cat("\n=== Installation Complete ===\n") 