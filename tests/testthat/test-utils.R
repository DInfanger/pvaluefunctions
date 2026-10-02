test_that("same_length() compares all arguments", {
  expect_true(same_length(1:3, 4:6, 7:9))
  expect_true(same_length(1:3))
  expect_false(same_length(1:3, 1:3, 1:2))
  expect_false(same_length(1:2, 1:3, 1:3))
  expect_false(same_length(1:3, 1:2))
})

test_that("is_identity_trans() recognises names and functions", {
  expect_true(is_identity_trans("identity"))
  expect_true(is_identity_trans("IDENTITY"))
  expect_true(is_identity_trans(identity))
  expect_false(is_identity_trans("exp"))
  expect_false(is_identity_trans(exp))
  expect_false(is_identity_trans(function(x) x))
  expect_false(is_identity_trans(NULL))
})

test_that("resolve_trans() handles the built-in transformations", {
  res <- resolve_trans("exp")
  expect_equal(res$name, "exp")
  expect_identical(res$fun, base::exp)

  res <- resolve_trans("Identity")
  expect_equal(res$name, "identity")
  expect_identical(res$fun, base::identity)

  expect_equal(resolve_trans(exp)$name, "exp")
  expect_equal(resolve_trans(identity)$name, "identity")
})

test_that("resolve_trans() accepts functions with any argument name", {
  res <- resolve_trans(function(y) y + 1)

  expect_equal(res$name, "custom")
  expect_equal(res$fun(1), 2)
})

test_that("resolve_trans() finds functions by name in the calling environment", {
  local_fun <- function(a) a * 10
  env <- environment()

  expect_identical(resolve_trans("local_fun", env = env)$fun, local_fun)
})

test_that("resolve_trans() keeps the case of the name but falls back to lower case", {
  MixedCase <- function(x) x + 1
  lower_only <- function(x) x + 2
  env <- environment()

  expect_equal(resolve_trans("MixedCase", env = env)$fun(0), 1)
  expect_equal(resolve_trans("LOWER_ONLY", env = env)$fun(0), 2)
})

test_that("resolve_trans() rejects invalid input", {
  expect_error(
    resolve_trans("no_such_function_anywhere"),
    "was not found"
  )
  expect_error(resolve_trans(1), "must be a function or the name of a function")
  expect_error(resolve_trans(c("exp", "identity")), "must be a function")
  expect_error(resolve_trans(NA_character_), "must be a function")
})
