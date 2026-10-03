# Variances, proportions and differences of proportions.

# Reference implementation of Newcombe's (1998) method 10: the intervals of the
# two proportions are Wilson score intervals with continuity correction, taken
# from stats::prop.test().
newcombe_cc <- function(x1, n1, x2, n2, level) {
  ci1 <- stats::prop.test(x1, n1, conf.level = level, correct = TRUE)$conf.int
  ci2 <- stats::prop.test(x2, n2, conf.level = level, correct = TRUE)$conf.int
  p1 <- x1 / n1
  p2 <- x2 / n2

  c(
    max(-1, p1 - p2 - sqrt((p1 - ci1[1])^2 + (ci2[2] - p2)^2)),
    min(1, p1 - p2 + sqrt((ci1[2] - p1)^2 + (p2 - ci2[1])^2))
  )
}

#-------------------------------------------------------------------------------
# Variance
#-------------------------------------------------------------------------------

test_that("variance confidence limits are chi-squared based", {
  levels <- c(0.95, 0.8)
  res <- cd(estimate = 4, n = 20, type = "var", conf_level = levels)
  df <- 19

  expect_equal(res$conf_frame$lwr, 4 * df / qchisq((1 + levels) / 2, df))
  expect_equal(res$conf_frame$upr, 4 * df / qchisq((1 - levels) / 2, df))
})

test_that("variance distribution matches the closed form", {
  res <- cd(estimate = 4, n = 20, type = "var", null_values = 3)
  rows <- res$res_frame
  df <- 19

  expect_equal(rows$conf_dist, 1 - pchisq(df * 4 / rows$values, df))
  expect_equal(rows$p_two, 1 - 2 * abs(rows$conf_dist - 0.5))
  expect_true(all(rows$values > 0))
})

test_that("variance counternull has the tail probability of the null value", {
  res <- cd(estimate = 4, n = 20, type = "var", null_values = c(3, 6))
  df <- 19
  cdf <- function(x) 1 - pchisq(df * 4 / x, df)

  expect_equal(
    1 - cdf(res$counternull_frame$counternull),
    cdf(res$counternull_frame$null_value)
  )
})

test_that("variance point estimates match the closed form", {
  n <- c(20, 30)
  est <- c(4, 6)
  df <- n - 1
  res <- cd(estimate = est, n = n, type = "var")

  expect_equal(res$point_est$est_mean, est * df / (df - 2))
  expect_equal(res$point_est$est_mode, est * df / (df + 2))
  expect_equal(res$point_est$est_median, est * df / qchisq(0.5, df))
})

test_that("the variance mean is missing for n <= 3", {
  # The mean of the scaled inverse chi-squared distribution needs df > 2
  res <- cd(estimate = c(4, 4, 4), n = c(2, 3, 4), type = "var")

  expect_equal(res$point_est$est_mean, c(NA, NA, 4 * 3))
})

#-------------------------------------------------------------------------------
# Proportions
#-------------------------------------------------------------------------------

test_that("proportion limits are Wilson score limits", {
  levels <- c(0.95, 0.8)
  res <- cd(estimate = 15 / 50, n = 50, type = "prop", conf_level = levels)

  for (i in seq_along(levels)) {
    wilson <- stats::prop.test(15, 50, conf.level = levels[i], correct = FALSE)
    expect_equal(c(res$conf_frame$lwr[i], res$conf_frame$upr[i]), wilson$conf.int[1:2])
  }
})

test_that("one-sided proportion limits use the matching two-sided level", {
  one_sided <- cd(
    estimate = 0.3,
    n = 50,
    type = "prop",
    conf_level = 0.95,
    alternative = "one_sided"
  )
  wilson <- stats::prop.test(15, 50, conf.level = 0.9, correct = FALSE)

  expect_equal(c(one_sided$conf_frame$lwr, one_sided$conf_frame$upr), wilson$conf.int[1:2])
})

test_that("proportion distribution is centred on the estimate", {
  res <- cd(estimate = 0.3, n = 50, type = "prop", null_values = 0.5)
  rows <- res$res_frame

  expect_equal(max(rows$p_two), 1)
  expect_equal(rows$values[which.max(rows$p_two)], 0.3)
  expect_true(all(rows$values > 0 & rows$values < 1))
  expect_true(all(res$counternull_frame$counternull < 0.3))
  expect_equal(nrow(res$aucc_frame), 1)
})

test_that("wilson_cicc() equals prop.test() with continuity correction", {
  for (n in c(10, 37, 120)) {
    # stats::prop.test() skips the continuity correction when the estimate is
    # exactly 0.5, so avoid that case
    for (successes in c(1, floor(n / 3), floor(n / 2) + 1, n - 2)) {
      for (level in c(0.8, 0.95)) {
        expected <- stats::prop.test(
          successes,
          n,
          conf.level = level,
          correct = TRUE
        )$conf.int[1:2]

        expect_equal(
          wilson_cicc(successes / n, n, level),
          expected,
          tolerance = 1e-8,
          label = paste("successes", successes, "n", n, "level", level)
        )
      }
    }
  }
})

test_that("wilson_ci() equals prop.test() without continuity correction", {
  expected <- stats::prop.test(7, 31, conf.level = 0.9, correct = FALSE)$conf.int[1:2]

  expect_equal(wilson_ci(7 / 31, 31, 0.9), expected)
})

#-------------------------------------------------------------------------------
# Difference of proportions
#-------------------------------------------------------------------------------

test_that("wilson_cicc_diff() follows Newcombe's method 10", {
  for (level in c(0.8, 0.95, 0.99)) {
    expect_equal(
      wilson_cicc_diff(c(56 / 70, 48 / 80), c(70, 80), level),
      newcombe_cc(56, 70, 48, 80, level),
      tolerance = 1e-8
    )
    expect_equal(
      wilson_cicc_diff(c(5 / 20, 12 / 30), c(20, 30), level),
      newcombe_cc(5, 20, 12, 30, level),
      tolerance = 1e-8
    )
  }
})

test_that("difference of proportions has the expected confidence limits", {
  res <- cd(
    estimate = c(68 / 100, 98 / 150),
    n = c(100, 150),
    type = "propdiff",
    conf_level = c(0.95, 0.8),
    n_values = 100L
  )

  expected <- rbind(
    newcombe_cc(68, 100, 98, 150, 0.95),
    newcombe_cc(68, 100, 98, 150, 0.8)
  )

  expect_equal(res$conf_frame$lwr, expected[, 1], tolerance = 1e-8)
  expect_equal(res$conf_frame$upr, expected[, 2], tolerance = 1e-8)
})

test_that("one-sided limits for a difference use the matching two-sided level", {
  one_sided <- cd(
    estimate = c(0.4, 0.3),
    n = c(60, 60),
    type = "propdiff",
    conf_level = 0.95,
    alternative = "one_sided",
    n_values = 100L
  )
  expected <- newcombe_cc(24, 60, 18, 60, 0.9)

  expect_equal(c(one_sided$conf_frame$lwr, one_sided$conf_frame$upr), expected, tolerance = 1e-8)
})

test_that("difference of proportions returns a p-value function only", {
  res <- cd(
    estimate = c(0.4, 0.3),
    n = c(60, 60),
    type = "propdiff",
    null_values = 0,
    n_values = 100L
  )
  rows <- res$res_frame

  # No confidence distribution or density for this type
  expect_true(all(is.na(rows$conf_dist)))
  expect_true(all(is.na(rows$conf_dens)))
  expect_false(all(is.na(rows$p_two)))
  expect_true(all(rows$p_two[!is.na(rows$p_two)] > 0))
  expect_equal(nrow(res$counternull_frame), 1)
})

test_that("the counternull of a difference of proportions is found", {
  estimate <- c(0.4, 0.3)
  n <- c(60, 60)
  res <- cd(
    estimate = estimate,
    n = n,
    type = "propdiff",
    null_values = c(-0.1, 0.05, 0.3, 1),
    n_values = 100L
  )
  counternull <- res$counternull_frame

  # The estimated difference is 0.1: the counternull lies on the other side
  expect_gt(counternull$counternull[1], 0.1)
  expect_gt(counternull$counternull[2], 0.1)
  expect_lt(counternull$counternull[3], 0.1)
  # A null value outside of the range of the confidence intervals has none
  expect_true(is.na(counternull$counternull[4]))

  # The null value and its counternull are the two limits of one interval
  for (j in 1:3) {
    null_value <- counternull$null_value[j]
    side <- if (null_value < 0.1) 1 else 2
    level <- uniroot(
      function(level) {
        wilson_cicc_diff(estimate, n, level)[side] - null_value
      },
      interval = c(1e-10, 1 - 1e-10)
    )$root

    expect_equal(
      wilson_cicc_diff(estimate, n, level)[3 - side],
      counternull$counternull[j],
      tolerance = 1e-2
    )
  }
})

test_that("null values in the gap of a difference of proportions have no counternull", {
  # The continuity correction leaves a gap around the estimated difference
  # (-0.1) in which uniroot() used to fail
  estimate <- c(0.3, 0.4)
  n <- c(10, 10)
  gap <- wilson_cicc_diff(estimate, n, 1e-15)
  null_values <- c(gap[1] + 0.01, estimate[1] - estimate[2], gap[2] - 0.01)

  res <- cd(estimate = estimate, n = n, type = "propdiff", null_values = null_values)

  expect_equal(res$counternull_frame$counternull, rep(NA_real_, 3))
})

test_that("the confidence density of a proportion matches the closed form", {
  # Wilson's score: z(x) = sqrt(n) * (x - p) / sqrt(x * (1 - x))
  p <- 0.3
  n <- 40
  res <- cd(estimate = p, n = n, type = "prop")
  rows <- res$res_frame
  x <- rows$values
  z <- sqrt(n) * (x - p) / sqrt(x * (1 - x))
  dz <- sqrt(n) * (x + p - 2 * p * x) / (2 * (x * (1 - x))^1.5)

  expect_equal(rows$conf_dist, pnorm(z))
  expect_equal(rows$conf_dens, dnorm(z) * dz)
  # The density integrates to 1
  expect_equal(trapz_sorted(x, rows$conf_dens), 1, tolerance = 1e-3)
})

test_that("the mean of a proportion is found for narrow distributions", {
  # integrate() over (0, 1) used to miss the peak and returned 0
  res <- cd(estimate = 0.3, n = 1e6, type = "prop")

  expect_equal(res$point_est$est_mean, 0.3, tolerance = 1e-5)
})

test_that("proportions of 0 and 1 give valid results", {
  # The z-score was 0 / 0 at the estimate, which gave NaN p-values
  for (successes in c(0, 20)) {
    res <- cd(
      estimate = successes / 20,
      n = 20,
      type = "prop",
      conf_level = 0.95,
      null_values = 0.5
    )
    rows <- res$res_frame

    expect_false(anyNA(rows$conf_dist))
    expect_false(anyNA(rows$p_two))
    expect_false(anyNA(rows$s_val))
    expect_true(all(rows$p_two >= 0 & rows$p_two <= 1))
    expect_equal(max(rows$p_two), 1)

    wilson <- stats::prop.test(successes, 20, correct = FALSE)$conf.int[1:2]
    expect_equal(c(res$conf_frame$lwr, res$conf_frame$upr), wilson)
  }

  # The distributions for 0 and 1 are mirror images
  mean_0 <- cd(estimate = 0, n = 20, type = "prop")$point_est$est_mean
  mean_1 <- cd(estimate = 1, n = 20, type = "prop")$point_est$est_mean
  expect_equal(mean_0, 1 - mean_1)
  expect_gt(mean_0, 0)
})

test_that("wilson_cicc() handles 0 and n successes without warnings", {
  # The square root of a negative number was taken (and silently dropped)
  for (successes in c(0, 20)) {
    expected <- stats::prop.test(successes, 20, correct = TRUE)$conf.int[1:2]

    expect_no_warning(limits <- wilson_cicc(successes / 20, 20, 0.95))
    expect_equal(limits, expected, tolerance = 1e-8)
  }
})

test_that("differences of proportions of 0 and 1 give no warnings", {
  for (estimate in list(c(0, 0), c(0, 1), c(1, 1))) {
    expect_no_warning(
      cd(
        estimate = estimate,
        n = c(20, 20),
        type = "propdiff",
        conf_level = 0.95,
        null_values = 0
      )
    )
  }
})
