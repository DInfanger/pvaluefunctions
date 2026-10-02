# Creates `regression.rds`, the reference output for test-regression.R.
#
# The fixture was created with the code of pvaluefunctions 1.6.3 (git commit
# fd4d419), i.e. BEFORE the refactoring of the development version, so that the
# tests check that the refactoring did not change any results. It should not be
# re-created with the current code unless a numerical change is intended.
#
# Usage (from the package root, with the 1.6.3 sources attached by
# pkgload::load_all() or installed):
#
#   Rscript tests/testthat/fixtures/make-regression.R

source("tests/testthat/fixtures/regression-configs.R")

grDevices::pdf(NULL)

result <- lapply(regression_configs(), function(args) {
  # Version 1.6.3 only accepts the name of the function
  if (is.function(args$trans)) {
    assign("rse_fun", args$trans, envir = globalenv())
    args$trans <- "rse_fun"
  }
  args$plot <- FALSE

  res <- suppressMessages(suppressWarnings(do.call(pvaluefunctions::conf_dist, args)))
  res[c("res_frame", "conf_frame", "counternull_frame", "point_est", "aucc_frame")]
})

saveRDS(result, "tests/testthat/fixtures/regression.rds", compress = "xz")
