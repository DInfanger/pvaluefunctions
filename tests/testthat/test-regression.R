# The reference output in fixtures/regression.rds was created with
# pvaluefunctions 1.6.3, i.e. before the refactoring (see
# fixtures/make-regression.R). These tests make sure that the results did not
# change.

reference <- readRDS(test_path("fixtures", "regression.rds"))
configs <- regression_configs()

# -log2(p) is dominated by floating point noise for the extreme tails of the
# exact Pearson distribution (p below about 1e-6)
mask_extreme_s_val <- function(frame) {
  frame$s_val[is.na(frame$p_two) | frame$p_two < 1e-6] <- NA
  frame
}

check_configuration <- function(name) {
  args <- configs[[name]]
  args$plot <- FALSE

  res <- suppressMessages(suppressWarnings(do.call(conf_dist, args)))
  ref <- reference[[name]]

  # Version 1.6.3 read the median of the t- and normal distributions off a
  # coarse grid. It is now exactly the estimate (see test-location-scale.R).
  location_scale_types <- c(
    "ttest", "linreg", "gammareg", "general_t",
    "logreg", "poisreg", "coxreg", "general_z"
  )
  if (args$type %in% location_scale_types) {
    ref$point_est$est_median <- res$point_est$est_median
  }

  # Missing point estimates (propdiff) were logical in version 1.6.3 and are
  # now numeric
  est_cols <- c("est_mean", "est_median", "est_mode")
  ref$point_est[est_cols] <- lapply(ref$point_est[est_cols], as.numeric)

  expect_equal(
    mask_extreme_s_val(res$res_frame),
    mask_extreme_s_val(ref$res_frame),
    tolerance = 1e-6,
    label = paste(name, "res_frame")
  )
  expect_equal(
    res$conf_frame,
    ref$conf_frame,
    tolerance = 1e-6,
    label = paste(name, "conf_frame")
  )
  expect_equal(
    res$counternull_frame,
    ref$counternull_frame,
    tolerance = 1e-6,
    label = paste(name, "counternull_frame")
  )
  expect_equal(
    res$point_est,
    ref$point_est,
    tolerance = 1e-6,
    label = paste(name, "point_est")
  )
  expect_equal(
    res$aucc_frame,
    ref$aucc_frame,
    tolerance = 1e-6,
    label = paste(name, "aucc_frame")
  )
}

test_that("results are unchanged compared with version 1.6.3", {
  # The exact distribution of Pearson's correlation is tested separately
  fast <- setdiff(names(configs), c("pearson", "pearson_one_sided"))

  for (name in fast) {
    check_configuration(name)
  }
})

test_that("exact Pearson results are unchanged compared with version 1.6.3", {
  skip_on_cran()

  check_configuration("pearson")
  check_configuration("pearson_one_sided")
})

test_that("the reference fixture covers all configurations", {
  expect_setequal(names(reference), names(configs))
})
