test_that("trapz_sorted() sorts and drops missing values", {
  x <- c(2, 0, NA, 1, 3)
  y <- c(2, 0, 5, 1, NA)

  # Remaining points (0, 0), (1, 1), (2, 2): the area under y = x on [0, 2]
  expect_equal(trapz_sorted(x, y), 2)
})

test_that("compute_aucc() integrates a triangular p-value function", {
  x <- seq(-1, 1, length.out = 201)
  res_frame <- data.frame(
    values = x,
    p_two = 1 - abs(x),
    variable = 1
  )

  no_null <- compute_aucc(res_frame, n_estimates = 1, null_values = NULL)

  expect_named(no_null, c("variable", "aucc"))
  expect_equal(no_null$aucc, 1)

  with_null <- compute_aucc(res_frame, n_estimates = 1, null_values = c(-1, 0, 0.5))

  expect_named(with_null, c("variable", "aucc", "null", "p_above_null"))
  expect_equal(with_null$aucc, c(1, 1, 1))
  expect_equal(with_null$null, c(-1, 0, 0.5))
  # Area above -1, 0 and 0.5: all, half, and 1/8 of the total area. Only grid
  # points strictly above the null value are used, so the strip between the
  # null value and the next grid point (width 0.01) is missing.
  expect_equal(with_null$p_above_null, c(1, 0.5, 0.125), tolerance = 0.02)
})

test_that("compute_aucc() handles several estimates", {
  x <- seq(-1, 1, length.out = 201)
  res_frame <- data.frame(
    values = c(x, 2 * x),
    p_two = c(1 - abs(x), 1 - abs(x)),
    variable = rep(1:2, each = length(x))
  )
  res <- compute_aucc(res_frame, n_estimates = 2, null_values = 0)

  expect_equal(res$variable, 1:2)
  expect_equal(res$aucc, c(1, 2))
  expect_equal(res$p_above_null, c(0.5, 0.5), tolerance = 0.03)
})
