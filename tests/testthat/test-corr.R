# Correlation coefficients: Fisher's z based approximations for Spearman and
# Kendall, the exact distribution for Pearson.

fisher_limits <- function(r, se, level) {
  q <- qnorm((1 + level) / 2)

  list(lwr = tanh(atanh(r) - q * se), upr = tanh(atanh(r) + q * se))
}

test_that("Spearman and Kendall use Fisher's z with their approximate se", {
  levels <- c(0.95, 0.8)

  spearman <- cd(
    estimate = 0.5,
    n = 40,
    type = "spearman",
    conf_level = levels
  )
  expected <- fisher_limits(0.5, sqrt((1 + 0.5^2 / 2) / (40 - 3)), levels)

  expect_equal(spearman$conf_frame$lwr, expected$lwr)
  expect_equal(spearman$conf_frame$upr, expected$upr)

  kendall <- cd(estimate = 0.4, n = 40, type = "kendall", conf_level = levels)
  expected <- fisher_limits(0.4, sqrt(0.437 / (40 - 4)), levels)

  expect_equal(kendall$conf_frame$lwr, expected$lwr)
  expect_equal(kendall$conf_frame$upr, expected$upr)
})

test_that("the Fisher z distribution matches the closed form", {
  res <- cd(
    estimate = 0.5,
    n = 40,
    type = "spearman",
    null_values = c(-0.2, 0.2)
  )
  rows <- res$res_frame
  se <- sqrt((1 + 0.5^2 / 2) / (40 - 3))
  z <- (atanh(0.5) - atanh(rows$values)) / se

  expect_equal(rows$conf_dist, 1 - pnorm(z))
  expect_equal(rows$p_two, 2 * (1 - pnorm(abs(z))))
  expect_true(all(abs(rows$values) < 1))
})

test_that("rank correlations have a counternull on Fisher's z scale", {
  res <- cd(
    estimate = 0.5,
    n = 40,
    type = "kendall",
    null_values = c(0, 0.2)
  )

  expect_equal(
    res$counternull_frame$counternull,
    tanh(2 * atanh(0.5) - atanh(c(0, 0.2)))
  )
})

test_that("correlation confidence limits are symmetric in the sign", {
  pos <- cd(estimate = 0.4, n = 30, type = "spearman", conf_level = 0.9)
  neg <- cd(estimate = -0.4, n = 30, type = "spearman", conf_level = 0.9)

  expect_equal(neg$conf_frame$lwr, -pos$conf_frame$upr)
  expect_equal(neg$conf_frame$upr, -pos$conf_frame$lwr)
})

test_that("exact Pearson distribution is a proper distribution", {
  skip_on_cran()

  res <- cd(estimate = 0.3, n = 30, type = "pearson", n_values = 400L)
  rows <- res$res_frame[order(res$res_frame$values), ]

  # Density integrates to one and agrees with the distribution function
  expect_equal(sum(diff(rows$values) * (rows$conf_dens[-1] + rows$conf_dens[-length(rows$conf_dens)]) / 2), 1, tolerance = 1e-2)
  expect_true(all(diff(rows$conf_dist) >= -1e-8))

  # Median and mode of the confidence distribution are close to the estimate
  expect_true(abs(res$point_est$est_mode - 0.3) < 0.02)
  expect_true(abs(res$point_est$est_median - 0.3) < 0.02)
  expect_true(abs(res$point_est$est_mean - 0.3) < 0.02)
})

test_that("exact Pearson p-values are not negative in the extreme tails", {
  skip_on_cran()

  # The integration error made the cdf slightly larger than 1 close to r = 1,
  # which gave negative p-values and NaN s-values
  expect_no_warning(
    res <- cd(estimate = 0.6, n = 15, type = "pearson", n_values = 500L)
  )
  rows <- res$res_frame

  expect_true(all(rows$conf_dist >= 0 & rows$conf_dist <= 1))
  expect_true(all(rows$p_two >= 0))
  expect_false(anyNA(rows$s_val))
})

test_that("exact Pearson confidence limits have the stated coverage", {
  skip_on_cran()

  levels <- c(0.95, 0.8)
  res <- cd(
    estimate = 0.3,
    n = 30,
    type = "pearson",
    n_values = 20L,
    conf_level = levels
  )
  limits <- res$conf_frame

  # The distribution function at the limits gives back the tail probabilities
  check <- cd(
    estimate = 0.3,
    n = 30,
    type = "pearson",
    n_values = 20L,
    null_values = c(limits$lwr, limits$upr)
  )$res_frame
  cdf_at <- function(x) {
    check$conf_dist[match(x, check$values)]
  }

  expect_equal(cdf_at(limits$lwr), (1 - levels) / 2, tolerance = 1e-3)
  expect_equal(cdf_at(limits$upr), 1 - (1 - levels) / 2, tolerance = 1e-3)
})

test_that("exact Pearson limits agree with Fisher's z for large n", {
  skip_on_cran()

  n <- 500
  res <- cd(
    estimate = 0.3,
    n = n,
    type = "pearson",
    n_values = 20L,
    conf_level = 0.95
  )
  expected <- fisher_limits(0.3, 1 / sqrt(n - 3), 0.95)

  expect_equal(res$conf_frame$lwr, expected$lwr, tolerance = 5e-3)
  expect_equal(res$conf_frame$upr, expected$upr, tolerance = 5e-3)
})

test_that("exact Pearson results are symmetric in the sign", {
  skip_on_cran()

  args <- list(n = 25, type = "pearson", n_values = 20L, conf_level = 0.9)
  pos <- do.call(cd, c(args, estimate = 0.4))
  neg <- do.call(cd, c(args, estimate = -0.4))

  expect_equal(neg$conf_frame$lwr, -pos$conf_frame$upr, tolerance = 1e-4)
  expect_equal(neg$conf_frame$upr, -pos$conf_frame$lwr, tolerance = 1e-4)
  expect_equal(neg$point_est$est_mean, -pos$point_est$est_mean, tolerance = 1e-6)
})

test_that("exact Pearson counternull has the same tail probability as the null", {
  skip_on_cran()

  args <- list(estimate = 0.3, n = 30, type = "pearson", n_values = 20L)
  first <- do.call(cd, c(args, null_values = 0))
  counternull <- first$counternull_frame$counternull

  expect_gt(counternull, 0.3)

  both <- do.call(cd, c(args, list(null_values = c(0, counternull))))$res_frame
  cdf_null <- both$conf_dist[match(0, both$values)]
  cdf_counternull <- both$conf_dist[match(counternull, both$values)]

  expect_equal(1 - cdf_counternull, cdf_null, tolerance = 1e-3)
})

test_that("exact Pearson works for large sample sizes (gamma overflow)", {
  skip_on_cran()

  for (n in c(200, 1000)) {
    res <- cd(
      estimate = 0.3,
      n = n,
      type = "pearson",
      n_values = 20L,
      conf_level = 0.95,
      null_values = 0
    )
    expected <- fisher_limits(0.3, 1 / sqrt(n - 3), 0.95)

    expect_true(all(is.finite(unlist(res$point_est))))
    expect_true(all(is.finite(unlist(res$conf_frame))))
    expect_true(is.finite(res$counternull_frame$counternull))
    expect_equal(res$conf_frame$lwr, expected$lwr, tolerance = 5e-3)
    expect_equal(res$conf_frame$upr, expected$upr, tolerance = 5e-3)
    # The mode must not be pushed to the boundary of the parameter space
    expect_true(abs(res$point_est$est_mode - 0.3) < 0.05)
    expect_gt(res$counternull_frame$counternull, 0.3)
  }
})

test_that("exact Pearson works for several estimates", {
  skip_on_cran()

  res <- cd(
    estimate = c(0.2, 0.6),
    n = c(25, 60),
    type = "pearson",
    n_values = 20L,
    conf_level = 0.95
  )

  expect_equal(nrow(res$point_est), 2)
  expect_equal(nrow(res$conf_frame), 2)
  expect_true(all(is.finite(unlist(res$point_est))))
  expect_gt(res$conf_frame$lwr[2], res$conf_frame$lwr[1])
  # The estimates are different, so are their sample sizes
  single <- cd(
    estimate = 0.6,
    n = 60,
    type = "pearson",
    n_values = 20L,
    conf_level = 0.95
  )
  expect_equal(res$conf_frame$lwr[2], single$conf_frame$lwr, tolerance = 1e-6)
})

test_that("exact Pearson works for very narrow distributions", {
  skip_on_cran()

  # integrate() used to miss the narrow peak of the density, so that the
  # confidence limits could not be found
  res <- cd(
    estimate = 0.99,
    n = 500,
    type = "pearson",
    n_values = 300L,
    conf_level = 0.95,
    null_values = 0.98
  )
  rows <- res$res_frame

  # Fisher's z is accurate for this sample size
  fisher <- tanh(atanh(0.99) + c(-1, 1) * qnorm(0.975) / sqrt(500 - 3))
  expect_equal(c(res$conf_frame$lwr, res$conf_frame$upr), fisher, tolerance = 1e-3)
  expect_equal(max(rows$conf_dist), 1, tolerance = 1e-4)
  expect_equal(res$point_est$est_mean, 0.99, tolerance = 1e-3)
  expect_gt(res$counternull_frame$counternull, 0.99)
})

test_that("the mean of the Fisher z distribution is found for narrow distributions", {
  # integrate() over (-1, 1) used to miss the peak and returned 0
  res <- cd(estimate = 0.8, n = 1e5, type = "spearman")

  expect_equal(res$point_est$est_mean, 0.8, tolerance = 1e-4)
})

# The exact Pearson tests above are skipped on CRAN; this one is fast and
# keeps the exact distribution covered there
test_that("exact Pearson results are plausible (also run on CRAN)", {
  res <- cd(
    estimate = 0.3,
    n = 30,
    type = "pearson",
    n_values = 40L,
    conf_level = 0.95,
    null_values = 0
  )
  fisher <- fisher_limits(0.3, 1 / sqrt(27), 0.95)

  expect_equal(res$conf_frame$lwr, fisher$lwr, tolerance = 0.05)
  expect_equal(res$conf_frame$upr, fisher$upr, tolerance = 0.05)
  expect_true(all(res$res_frame$p_two >= 0 & res$res_frame$p_two <= 1))
  expect_gt(res$counternull_frame$counternull, 0.3)
})
