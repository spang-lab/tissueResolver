# This script runs tests and updates documentation if all tests pass. To run
# this script got to the root directory and call Rscript tests/validate_pkg.R
# from the terminal.

# load dependencies
suppressPackageStartupMessages({
  library(devtools)
  library(testthat)
})

# Run tests
test(stop_on_failure = TRUE)
document()
