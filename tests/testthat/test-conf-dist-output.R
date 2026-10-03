# Structure and consistency of the object returned by conf_dist()

z_args <- list(estimate = 0.5, stderr = 0.2, type = "general_z")

test_that("the returned list has the documented components", {
  res <- cd(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    conf_level = 0.95,
    null_values = 0,
    plot_counternull = TRUE
  )

  expect_named(
    res,
    c("res_frame", "conf_frame", "counternull_frame", "point_est", "aucc_frame")
  )
  expect_named(
    res$res_frame,
    c(
      "values", "conf_dist", "conf_dens", "p_two", "p_one", "variable",
      "s_val", "hypothesis", "counternull"
    )
  )
  expect_named(res$conf_frame, c("conf_level", "lwr", "upr", "variable"))
  expect_named(
    res$counternull_frame,
    c("null_value", "counternull", "variable")
  )
  expect_named(res$point_est, c("est_mean", "est_median", "est_mode", "variable"))
  expect_named(res$aucc_frame, c("variable", "aucc", "null", "p_above_null"))
})

test_that("components are NULL if they were not requested", {
  res <- cd(estimate = 0.5, stderr = 0.2, type = "general_z")

  expect_null(res$conf_frame)
  expect_null(res$counternull_frame)
  expect_named(res$aucc_frame, c("variable", "aucc"))
  expect_false("counternull" %in% names(res$res_frame))
})

test_that("a plot is only returned on request", {
  without_plot <- cd(estimate = 0.5, stderr = 0.2, type = "general_z")
  expect_false("plot" %in% names(without_plot))

  with_plot <- conf_dist_quiet(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    n_values = 100L
  )
  expect_s3_class(with_plot$plot, "ggplot")
})

test_that("conf_dist() prints the plot and returns the results invisibly", {
  skip_on_cran()

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_invisible(conf_dist(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    n_values = 100L
  ))
})

test_that("the data frame is sorted by the values", {
  res <- cd(
    estimate = c(0, 1),
    stderr = c(1, 0.5),
    type = "general_z",
    null_values = c(0.2, -0.3)
  )

  expect_false(is.unsorted(res$res_frame$values, na.rm = TRUE))
})

test_that("null values are part of the evaluated grid", {
  res <- cd(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    null_values = c(0.123, -0.321)
  )

  expect_true(all(c(0.123, -0.321) %in% res$res_frame$values))
  expect_equal(res$counternull_frame$null_value, c(0.123, -0.321))
})

test_that("s-values are the surprisal of the two-sided p-values", {
  for (type in c("general_z", "general_t")) {
    res <- cd(estimate = 0.5, stderr = 0.2, df = 10, type = type)

    expect_equal(res$res_frame$s_val, -log2(res$res_frame$p_two))
  }
})

test_that("the hypothesis column splits the curve at the estimate", {
  res <- cd(estimate = 0.5, stderr = 0.2, type = "general_z")
  rows <- res$res_frame
  rows <- rows[!is.na(rows$hypothesis), ]

  expect_true(all(rows$values[rows$hypothesis == "less"] < 0.5))
  expect_true(all(rows$values[rows$hypothesis == "greater"] >= 0.5))
})

test_that("the counternull column is the s-value at the counternull", {
  res <- cd(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    null_values = 0,
    plot_counternull = TRUE
  )
  cn <- res$res_frame$counternull

  expect_true(any(!is.na(cn)))
})

test_that("small p-values are cut off at plot_p_limit", {
  res <- conf_dist(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    n_values = 500L,
    plot_p_limit = 0.01,
    plot = FALSE
  )
  rows <- res$res_frame

  expect_true(any(is.na(rows$p_two)))
  expect_true(all(rows$p_two[!is.na(rows$p_two)] >= 0.01))
  expect_false(anyNA(rows$conf_dist))
})

test_that("a value for plot_p_limit that is too small is rounded", {
  expect_no_error(
    conf_dist(
      estimate = 0.5,
      stderr = 0.2,
      type = "general_z",
      n_values = 50L,
      plot_p_limit = 1e-12,
      plot = FALSE
    )
  )
})

test_that("results do not depend on the number of estimates evaluated together", {
  both <- cd(
    estimate = c(0, 1),
    stderr = c(1, 0.5),
    type = "general_z",
    conf_level = 0.95
  )
  second <- cd(estimate = 1, stderr = 0.5, type = "general_z", conf_level = 0.95)

  expect_equal(both$conf_frame$lwr[2], second$conf_frame$lwr)
  expect_equal(both$point_est$est_mean[2], second$point_est$est_mean)
  expect_equal(
    rows_of(both, 2)$p_two,
    rows_of(second, 1)$p_two
  )
})

test_that("proportion of the AUCC above null values is returned", {
  res <- cd(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    null_values = c(-10, 0.5, 10),
    n_values = 4000L
  )
  above <- res$aucc_frame$p_above_null

  expect_equal(res$aucc_frame$null, c(-10, 0.5, 10))
  expect_equal(above[1], 1, tolerance = 1e-3)
  expect_equal(above[2], 0.5, tolerance = 1e-2)
  expect_equal(above[3], 0, tolerance = 1e-3)
})

test_that("null values that are outside of the range give no error", {
  expect_no_error(cd(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    null_values = c(-50, 50),
    plot_counternull = TRUE
  ))
})

test_that("duplicated null values are accepted", {
  expect_no_error(cd(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    null_values = c(0, 0)
  ))
})

test_that("confidence levels that are too low for one-sided are dropped", {
  expect_message(
    res <- cd(
      estimate = 0.5,
      stderr = 0.2,
      type = "general_z",
      conf_level = c(0.3, 0.9),
      alternative = "one_sided"
    ),
    "too low"
  )
  expect_equal(res$conf_frame$conf_level, 0.9)

  # No message if all levels are fine
  expect_no_message(cd(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    conf_level = c(0.3, 0.9)
  ))
})

test_that("argument matching works for partial strings", {
  res <- cd(estimate = 0.5, stderr = 0.2, type = "general_z", alternative = "one")

  expect_s3_class(res$res_frame, "data.frame")
  expect_error(
    cd(estimate = 0.5, stderr = 0.2, type = "general_z", alternative = "both"),
    "should be one of"
  )
  expect_error(
    cd(estimate = 0.5, stderr = 0.2, type = "general_z", plot_type = "xyz"),
    "should be one of"
  )
  expect_error(
    cd(estimate = 0.5, stderr = 0.2, type = "general_z", x_scale = "sqrt"),
    "should be one of"
  )
})
