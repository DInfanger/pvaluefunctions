# t and normal based confidence distributions. All expected values are
# computed in closed form with stats functions.

est <- 0.5
se <- 0.2

test_that("normal confidence distribution matches the closed form", {
  res <- cd(
    estimate = est,
    stderr = se,
    type = "general_z",
    null_values = c(0, 0.3)
  )
  rows <- res$res_frame
  z <- (rows$values - est) / se

  expect_equal(rows$conf_dist, pnorm(z))
  expect_equal(rows$conf_dens, dnorm(z) / se)
  expect_equal(rows$p_two, 2 * (1 - pnorm(abs(z))))
  expect_equal(rows$p_one, pmin(pnorm(z), 1 - pnorm(z)))
  expect_equal(rows$s_val, -log2(rows$p_two))
})

test_that("t confidence distribution matches the closed form", {
  res <- cd(
    estimate = est,
    stderr = se,
    df = 12,
    type = "general_t",
    null_values = 0
  )
  rows <- res$res_frame
  t_val <- (rows$values - est) / se

  expect_equal(rows$conf_dist, pt(t_val, df = 12))
  expect_equal(rows$conf_dens, dt(t_val, df = 12) / se)
  expect_equal(rows$p_two, 2 * (1 - pt(abs(t_val), df = 12)))
})

test_that("linreg and gammareg are t based like general_t", {
  args <- list(estimate = est, stderr = se, df = 12, null_values = 0)
  general_t <- do.call(cd, c(args, type = "general_t"))

  expect_equal(do.call(cd, c(args, type = "linreg")), general_t)
  expect_equal(do.call(cd, c(args, type = "gammareg")), general_t)
})

test_that("logreg, poisreg and coxreg are normal based like general_z", {
  args <- list(estimate = est, stderr = se, null_values = 0)
  general_z <- do.call(cd, c(args, type = "general_z"))

  expect_equal(do.call(cd, c(args, type = "logreg")), general_z)
  expect_equal(do.call(cd, c(args, type = "poisreg")), general_z)
  expect_equal(do.call(cd, c(args, type = "coxreg")), general_z)
})

test_that("the t-test is a t distribution with the implied standard error", {
  tstat <- 2.1
  res <- cd(estimate = 0.8, tstat = tstat, df = 25, type = "ttest")
  general_t <- cd(estimate = 0.8, stderr = 0.8 / tstat, df = 25, type = "general_t")

  expect_equal(res, general_t)
})

test_that("p-value at the estimate is 1 and the density peaks there", {
  res <- cd(estimate = est, stderr = se, type = "general_z")
  rows <- res$res_frame

  expect_equal(max(rows$p_two), 1)
  expect_equal(rows$values[which.max(rows$p_two)], est)
  expect_equal(rows$conf_dist[rows$values == est], 0.5)
})

test_that("two-sided confidence limits are estimate +- quantile * se", {
  levels <- c(0.95, 0.8, 0.5)
  res <- cd(estimate = est, stderr = se, type = "general_z", conf_level = levels)

  expect_equal(res$conf_frame$conf_level, levels)
  expect_equal(res$conf_frame$lwr, est - qnorm((1 + levels) / 2) * se)
  expect_equal(res$conf_frame$upr, est + qnorm((1 + levels) / 2) * se)

  res_t <- cd(
    estimate = est,
    stderr = se,
    df = 9,
    type = "general_t",
    conf_level = levels
  )

  expect_equal(res_t$conf_frame$lwr, est - qt((1 + levels) / 2, df = 9) * se)
  expect_equal(res_t$conf_frame$upr, est + qt((1 + levels) / 2, df = 9) * se)
})

test_that("one-sided confidence limits use the one-sided quantiles", {
  levels <- c(0.95, 0.8)
  res <- cd(
    estimate = est,
    stderr = se,
    type = "general_z",
    conf_level = levels,
    alternative = "one_sided"
  )

  expect_equal(res$conf_frame$lwr, est + qnorm(1 - levels) * se)
  expect_equal(res$conf_frame$upr, est + qnorm(levels) * se)

  # A one-sided 95% interval has the limits of a two-sided 90% interval
  two_sided <- cd(estimate = est, stderr = se, type = "general_z", conf_level = 0.9)
  expect_equal(res$conf_frame$lwr[1], two_sided$conf_frame$lwr)
  expect_equal(res$conf_frame$upr[1], two_sided$conf_frame$upr)
})

test_that("the counternull mirrors the null value around the estimate", {
  nulls <- c(0, 0.3, 0.9)

  for (type in c("general_z", "general_t")) {
    res <- cd(
      estimate = est,
      stderr = se,
      df = 14,
      type = type,
      null_values = nulls
    )

    expect_equal(res$counternull_frame$null_value, nulls)
    expect_equal(res$counternull_frame$counternull, 2 * est - nulls)
  }
})

test_that("mean, median and mode of the distribution equal the estimate", {
  # Median and mode are read off a grid, so use a fine one
  res <- cd(estimate = est, stderr = se, type = "general_t", df = 20, n_values = 2000L)

  expect_equal(res$point_est$est_mean, est)
  expect_equal(res$point_est$est_median, est, tolerance = 1e-2)
  expect_equal(res$point_est$est_mode, est, tolerance = 1e-2)
})

test_that("the area under the confidence curve matches the closed form", {
  # The area under a two-sided normal p-value function is 2 * E|Z| * se
  res <- cd(estimate = est, stderr = se, type = "general_z", n_values = 4000L)

  expect_equal(res$aucc_frame$aucc, 2 * sqrt(2 / pi) * se, tolerance = 1e-3)
})

test_that("a transformation is applied to values, limits and counternulls", {
  args <- list(
    estimate = est,
    stderr = se,
    type = "logreg",
    conf_level = 0.95,
    null_values = 0
  )
  identity_res <- do.call(cd, args)
  exp_res <- do.call(cd, c(args, trans = "exp"))

  expect_equal(exp_res$res_frame$values, exp(identity_res$res_frame$values))
  expect_equal(exp_res$conf_frame$lwr, exp(identity_res$conf_frame$lwr))
  expect_equal(exp_res$conf_frame$upr, exp(identity_res$conf_frame$upr))
  expect_equal(
    exp_res$counternull_frame$counternull,
    exp(identity_res$counternull_frame$counternull)
  )
  expect_equal(exp_res$res_frame$conf_dist, identity_res$res_frame$conf_dist)
})

test_that("trans can be a name, a different case or a function", {
  args <- list(estimate = est, stderr = se, type = "logreg", conf_level = 0.95)
  reference <- do.call(cd, c(args, trans = "exp"))

  expect_equal(do.call(cd, c(args, trans = "EXP")), reference)
  expect_equal(do.call(cd, c(args, trans = exp)), reference)
  expect_equal(do.call(cd, c(args, trans = function(y) exp(y))), reference)

  # Names are looked up in the environment that calls conf_dist(), so call it
  # directly instead of through the cd() wrapper
  scale_ten <- function(value) value * 10
  by_name <- conf_dist(
    estimate = est,
    stderr = se,
    type = "logreg",
    trans = "scale_ten",
    n_values = 200L,
    plot_p_limit = 0,
    plot = FALSE
  )
  expect_equal(
    by_name$res_frame$values,
    rows_of(do.call(cd, args))$values * 10
  )
})

test_that("an unknown transformation gives an error", {
  expect_error(
    cd(estimate = est, stderr = se, type = "logreg", trans = "no_such_function"),
    "was not found"
  )
})

test_that("several estimates are returned in one data frame", {
  res <- cd(
    estimate = c(0, 1, 2),
    stderr = c(1, 0.5, 2),
    type = "general_z",
    est_names = c("a", "b", "c"),
    conf_level = 0.95
  )

  expect_equal(levels(res$res_frame$variable), c("a", "b", "c"))
  expect_equal(levels(res$conf_frame$variable), c("a", "b", "c"))
  expect_equal(nrow(res$point_est), 3)
  expect_equal(res$conf_frame$lwr, c(0, 1, 2) - qnorm(0.975) * c(1, 0.5, 2))
  expect_equal(nrow(res$aucc_frame), 3)

  # Narrower distributions have smaller areas
  expect_equal(
    res$aucc_frame$aucc,
    2 * sqrt(2 / pi) * c(1, 0.5, 2),
    tolerance = 1e-2
  )
})

test_that("default estimate names are the indices", {
  res <- cd(estimate = c(0, 1), stderr = c(1, 1), type = "general_z")

  expect_equal(levels(res$res_frame$variable), c("1", "2"))
})
