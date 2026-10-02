# Input validation. Errors are checked with a short regular expression for the
# key phrase of the message; a few representative messages are snapshotted.

z_args <- list(estimate = 0.5, stderr = 0.2, type = "general_z")

# Calls cd() with the base arguments modified by `...`
cd_z <- function(...) {
  args <- utils::modifyList(z_args, list(...))
  do.call(cd, args)
}

test_that("valid edge cases are accepted", {
  expect_no_error(cd_z(conf_level = numeric(0)))
  expect_no_error(cd_z(df = Inf, type = "linreg"))
  expect_no_error(cd_z(est_names = "x"))
  expect_no_error(cd_z(plot_p_limit = 0))
})

test_that("type and estimate are required and checked", {
  expect_error(cd(stderr = 1, type = "general_z"), "provide an estimate")
  expect_error(cd(estimate = 1, stderr = 1), "type of the estimate")
  expect_error(cd_z(type = "nonsense"), "must be one of")
  expect_error(cd_z(type = c("general_z", "ttest")), "must be one of")
  expect_error(cd_z(type = 1), "must be one of")
  expect_error(cd_z(estimate = "a"), "numeric vector")
  expect_error(cd_z(estimate = NA_real_), "missing or infinite")
  expect_error(cd_z(estimate = Inf), "missing or infinite")
})

test_that("type is not case sensitive", {
  expect_no_error(cd_z(type = "GENERAL_Z"))
})

test_that("logical flags must be TRUE or FALSE", {
  expect_error(cd_z(together = NA), "must be either")
  expect_error(cd_z(together = "yes"), "must be either")
  expect_error(cd_z(inverted = c(TRUE, FALSE)), "must be either")
  expect_error(cd_z(log_yaxis = 1), "must be either")
  expect_error(
    cd(estimate = 0.5, stderr = 0.2, type = "general_z", plot_legend = NULL),
    "must be either"
  )
  expect_error(cd_z(same_color = "no"), "must be either")
  expect_error(cd_z(plot_counternull = NA), "must be either")
})

test_that("numeric options are range checked", {
  expect_error(cd_z(n_values = 1), "single number")
  expect_error(cd_z(n_values = NA), "single number")
  expect_error(cd_z(n_values = c(10, 20)), "single number")
  expect_error(cd_z(plot_p_limit = 1), "single number")
  expect_error(cd_z(plot_p_limit = -0.1), "single number")
  expect_error(cd_z(cut_logyaxis = 0), "single number")
  expect_error(cd_z(cut_logyaxis = 1.5), "single number")
  expect_error(cd_z(nrow = 0), "whole number")
  expect_error(cd_z(nrow = 1.5), "whole number")
  expect_error(cd_z(ncol = "a"), "whole number")
})

test_that("p-value limits are consistent with the options", {
  expect_error(
    cd_z(plot_p_limit = 0, log_yaxis = TRUE),
    "Cannot plot 0 on logarithmic axis"
  )
  expect_error(
    cd_z(plot_p_limit = 0.5, alternative = "one_sided"),
    "below 0.5"
  )
})

test_that("confidence levels, null values, labels and limits are checked", {
  expect_error(cd_z(conf_level = 1), "between 0 and 1")
  expect_error(cd_z(conf_level = c(0.9, 1.2)), "between 0 and 1")
  expect_error(cd_z(conf_level = NA_real_), "missing or infinite")
  expect_error(cd_z(conf_level = "a"), "numeric vector")
  expect_error(cd_z(null_values = "a"), "numeric vector")
  expect_error(cd_z(null_values = NA_real_), "missing or infinite")
  expect_error(cd_z(xlab = c("a", "b")), "x-axis label")
  expect_error(cd_z(xlim = 1), "two limits")
  expect_error(cd_z(xlim = c(0, 1, 2)), "two limits")
  expect_error(cd_z(xlim = c(0, NA)), "x-axis limits")
  expect_error(cd_z(xlim = c(0, Inf)), "x-axis limits")
})

test_that("required arguments by type are enforced", {
  expect_error(cd(estimate = 1, type = "ttest", df = 10), "t-statistic")
  expect_error(cd(estimate = 1, type = "ttest", tstat = 2), "t-statistic")
  expect_error(cd(estimate = 1, type = "linreg", df = 10), "standard error")
  expect_error(cd(estimate = 1, type = "linreg", stderr = 1), "degrees of freedom")
  expect_error(cd(estimate = 1, type = "general_t", stderr = 1), "degrees of freedom")
  expect_error(cd(estimate = 1, type = "logreg"), "standard error")
  expect_error(cd(estimate = 1, type = "poisreg"), "standard error")
  expect_error(cd(estimate = 1, type = "coxreg"), "standard error")
  expect_error(cd(estimate = 0.5, type = "pearson"), "sample size")
  expect_error(cd(estimate = 0.5, type = "prop"), "sample size")
  expect_error(cd(estimate = 4, type = "var"), "Sample size")
})

test_that("stderr, df and n must be positive", {
  expect_error(cd_z(stderr = 0), "larger than 0")
  expect_error(cd_z(stderr = -1), "larger than 0")
  expect_error(cd_z(stderr = NA_real_), "missing or infinite")
  expect_error(cd_z(stderr = "a"), "numeric vector")
  expect_error(cd_z(type = "linreg", df = 0), "larger than 0")
  expect_error(cd_z(type = "linreg", df = NA_real_), "missing")
  expect_error(cd(estimate = 4, n = 0, type = "var"), "larger than 0")
  expect_error(cd(estimate = 0.3, n = -5, type = "prop"), "larger than 0")
})

test_that("lengths of the arguments must match the estimates", {
  expect_error(
    cd_z(estimate = c(1, 2), stderr = 1),
    "same length as estimates"
  )
  expect_error(
    cd(estimate = c(1, 2), type = "linreg", stderr = c(1, 1), df = 10),
    "same length as estimates"
  )
  # The third argument used to be ignored
  expect_error(
    cd(estimate = c(1, 2), type = "linreg", stderr = c(1, 1), df = c(10, 10, 10)),
    "same length as estimates"
  )
  expect_error(
    cd(estimate = c(1, 2), type = "ttest", tstat = c(2, 3), df = c(10, 10, 10)),
    "same length as estimates"
  )
  expect_error(
    cd(estimate = c(1, 2), type = "ttest", tstat = c(2, 3, 4), df = c(10, 10)),
    "same length as estimates"
  )
  expect_error(
    cd(estimate = c(0.3, 0.4), type = "prop", n = 50),
    "same length as estimates"
  )
})

test_that("estimate names must match and be unique", {
  expect_error(
    cd_z(estimate = c(1, 2), stderr = c(1, 1), est_names = "a"),
    "estimate names"
  )
  expect_error(
    cd_z(estimate = c(1, 2), stderr = c(1, 1), est_names = c("a", "a")),
    "must be unique"
  )
})

test_that("nrow * ncol must be large enough for separate panels", {
  expect_error(
    cd_z(
      estimate = c(1, 2, 3),
      stderr = c(1, 1, 1),
      together = FALSE,
      nrow = 1,
      ncol = 2
    ),
    "nrow \\* ncol"
  )
  expect_no_error(
    cd_z(
      estimate = c(1, 2, 3),
      stderr = c(1, 1, 1),
      together = TRUE,
      nrow = 1,
      ncol = 2
    )
  )
})

test_that("the t-statistic must be consistent with the estimate", {
  ttest <- function(estimate, tstat) {
    cd(estimate = estimate, tstat = tstat, df = 20, type = "ttest")
  }

  expect_no_error(ttest(0.8, 2.1))
  expect_no_error(ttest(-0.8, -2.1))
  expect_error(ttest(0.8, -2.1), "same sign")
  expect_error(ttest(0.8, 0), "same sign")
  expect_error(ttest(0, 2), "same sign")
})

test_that("variance estimates are checked", {
  expect_error(cd(estimate = 0, n = 10, type = "var"), "larger than 0")
  expect_error(cd(estimate = -1, n = 10, type = "var"), "larger than 0")
  expect_error(cd(estimate = 4, n = 1, type = "var"), "larger than 1")
})

test_that("proportions are checked", {
  expect_error(cd(estimate = 1.2, n = 10, type = "prop"), "between 0 and 1")
  expect_error(cd(estimate = -0.1, n = 10, type = "prop"), "between 0 and 1")
  expect_error(
    cd(estimate = 0.3, n = 10, type = "prop", null_values = 1),
    "excluding"
  )
  expect_error(
    cd(estimate = 0.3, n = 10, type = "prop", null_values = 0),
    "excluding"
  )
})

test_that("differences of proportions are checked", {
  expect_error(
    cd(estimate = 0.3, n = 10, type = "propdiff"),
    "exactly two"
  )
  expect_error(
    cd(estimate = c(0.3, 0.4), n = 10, type = "propdiff"),
    "exactly two"
  )
  expect_error(
    cd(estimate = c(0.3, 0.4), n = c(10, 10), type = "propdiff", est_names = c("a", "b")),
    "only one estimate name"
  )
  expect_error(
    cd(estimate = c(0.3, 0.4), n = c(10, 10), type = "propdiff", plot_type = "pdf"),
    "only P-value functions"
  )
  expect_warning(
    cd(estimate = c(0.333, 0.4), n = c(10, 10), type = "propdiff"),
    "not integer"
  )
})

test_that("correlations are checked", {
  for (type in c("pearson", "spearman", "kendall")) {
    expect_error(cd(estimate = 1, n = 20, type = type), "strictly between")
    expect_error(cd(estimate = -1, n = 20, type = type), "strictly between")
    expect_error(
      cd(estimate = 0.3, n = 20, type = type, null_values = 1.5),
      "between -1 and 1"
    )
  }

  expect_error(cd(estimate = 0.3, n = 3, type = "pearson"), "at least 4")
  expect_error(cd(estimate = 0.3, n = 3, type = "spearman"), "at least 4")
  expect_error(cd(estimate = 0.3, n = 4, type = "kendall"), "at least 5")
})

test_that("approximations for rank correlations warn", {
  expect_warning(
    cd(estimate = 0.95, n = 30, type = "spearman"),
    "Spearman"
  )
  expect_warning(
    cd(estimate = 0.3, n = 8, type = "spearman"),
    "Spearman"
  )
  expect_warning(
    cd(estimate = 0.85, n = 30, type = "kendall"),
    "Kendall"
  )
})

test_that("conf_dist() ignores confidence levels that are too low", {
  expect_message(
    res <- cd_z(conf_level = c(0.2, 0.9), alternative = "one_sided"),
    "too low"
  )
  expect_equal(res$conf_frame$conf_level, 0.9)
})

test_that("representative error messages are stable", {
  expect_snapshot(cd_z(stderr = -1), error = TRUE)
  expect_snapshot(cd_z(together = NA), error = TRUE)
  expect_snapshot(
    cd_z(estimate = c(1, 2), stderr = c(1, 1), est_names = c("a", "a")),
    error = TRUE
  )
})
