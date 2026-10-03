# Confidence distributions for proportions and differences of proportions.

wilson_ci <- function(estimate, n, conf_level) {
  z <- qnorm((conf_level + 1) / 2)

  p1 <- estimate + (1 / 2) * z^2 / n
  p2 <- z * sqrt((estimate * (1 - estimate) + (1 / 4) * z^2 / n) / n)
  p3 <- 1 + z^2 / n

  c((p1 - p2) / p3, (p1 + p2) / p3)
}

wilson_cicc <- function(estimate, n, conf_level) {
  z <- qnorm((conf_level + 1) / 2)
  x <- round(estimate * n) # To get number of successes/failures
  estimate_compl <- (1 - estimate) # Complement of estimate

  # The limits are 0 for x = 0 and 1 for x = n by definition (Newcombe 1998).
  # The formulas would take the square root of a negative number there.
  if (x == 0) {
    lower <- 0
  } else {
    lower <- max(
      0,
      (2 *
        x +
        z^2 -
        1 -
        z * sqrt(z^2 - 2 - 1 / n + 4 * estimate * (n * estimate_compl + 1))) /
        (2 * (n + z^2))
    )
  }

  if (x == n) {
    upper <- 1
  } else {
    upper <- min(
      1,
      (2 *
        x +
        z^2 +
        1 +
        z * sqrt(z^2 + 2 - 1 / n + 4 * estimate * (n * estimate_compl - 1))) /
        (2 * (n + z^2))
    )
  }

  c(lower, upper)
}

wilson_cicc_diff <- function(estimate, n, conf_level) {
  est_diff <- (estimate[1] - estimate[2])

  res1 <- wilson_cicc(estimate = estimate[1], n = n[1], conf_level = conf_level)
  res2 <- wilson_cicc(estimate = estimate[2], n = n[2], conf_level = conf_level)

  l1 <- res1[1]
  u1 <- res1[2]
  l2 <- res2[1]
  u2 <- res2[2]

  lim1 <- max(-1, est_diff - sqrt((estimate[1] - l1)^2 + (u2 - estimate[2])^2))
  lim2 <- min(1, est_diff + sqrt((u1 - estimate[1])^2 + (estimate[2] - l2)^2))

  sort(c(lim1, lim2), decreasing = FALSE)
}

cdist_prop1 <- function(
  estimate = NULL,
  n = NULL,
  n_values = NULL,
  conf_level = NULL,
  null_values = NULL,
  alternative = NULL
) {
  # Auxiliary functions: the z-score of Wilson's interval as a function of the
  # true proportion x and its derivative with respect to x

  # For an estimate of 0 or 1, the grid contains x = p at the boundary, where
  # the formulas give 0 / 0. The z-score tends to 0 there and the density
  # diverges.
  cdf_fun <- function(x, n, p) {
    ifelse(x == p, 0, sqrt(n) * (x - p) / sqrt(x * (1 - x)))
  }

  deriv_fun <- function(x, n, p) {
    ifelse(
      x * (1 - x) == 0,
      Inf,
      sqrt(n) * (x + p - 2 * p * x) / (2 * (x * (1 - x))^1.5)
    )
  }

  # Inverse of the z-score: the true proportion x with z-score z (a limit of
  # Wilson's interval). Rounding can push it just outside of [0, 1] for
  # estimates of 0 or 1.
  wilson_inverse <- function(z, n, p) {
    x <- (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n)
    pmin(pmax(x, 0), 1)
  }

  eps <- 1e-10

  res_list <- list()

  conf_list <- list()

  counternull_list <- list()

  for (i in seq_along(estimate)) {
    limits <- wilson_ci(estimate = estimate[i], n = n[i], conf_level = 1 - eps)

    x_calc <- c(
      estimate[i],
      null_values,
      seq(limits[1], limits[2], length.out = n_values)
    )

    z_calc <- cdf_fun(x = x_calc, n = n[i], p = estimate[i])
    cdf_calc <- pnorm(z_calc)

    res_list[[length(res_list) + 1]] <- cdist_res_matrix(
      x = x_calc,
      cdf = cdf_calc,
      dens = dnorm(z_calc) * deriv_fun(x = x_calc, n = n[i], p = estimate[i]),
      i = i
    )

    # Confidence intervals

    if (!is.null(conf_level)) {
      conf_list[[length(conf_list) + 1]] <- cdist_conf_matrix(
        conf_level = conf_level,
        limits = wilson_ci(
          estimate = estimate[i],
          n = n[i],
          conf_level = two_sided_level(conf_level, alternative)
        ),
        i = i
      )
    }

    # Counternulls

    if (!is.null(null_values)) {
      # The counternull has the z-score of the null value with opposite sign
      counternull_list[[length(counternull_list) + 1]] <-
        cdist_counternull_matrix(
          null_values = null_values,
          counternull = wilson_inverse(
            -cdf_fun(x = null_values, n = n[i], p = estimate[i]),
            n = n[i],
            p = estimate[i]
          ),
          i = i
        )
    }
  }

  frames <- assemble_cdist_frames(
    res_list = res_list,
    conf_list = conf_list,
    counternull_list = counternull_list,
    conf_level = conf_level,
    null_values = null_values
  )

  # Point estimators

  point_est_frame <- empty_point_est_frame(length(estimate))

  # The mean is E[x(Z)] with Z ~ N(0, 1). Integrating on the z scale cannot
  # miss the peak of a narrow density (large n) and also covers the point mass
  # of 1/2 at the boundary for estimates of 0 or 1. The median is the estimate
  # (z = 0).
  point_est_frame$est_median <- estimate

  for (i in seq_along(estimate)) {
    point_est_frame$est_mean[i] <- integrate(
      function(z) wilson_inverse(z, n = n[i], p = estimate[i]) * dnorm(z),
      lower = -10,
      upper = 10,
      rel.tol = 1e-10
    )$value
    point_est_frame$est_mode[i] <- point_est_mode(frames$res_frame, i)
  }

  c(frames, list(point_est = point_est_frame))
}

cdist_propdiff <- function(
  estimate = NULL,
  n = NULL,
  n_values = NULL,
  conf_level = NULL,
  null_values = NULL,
  alternative = NULL
) {
  eps <- 1e-15
  est_diff <- estimate[1] - estimate[2]

  # Two-sided p-values only: each confidence level gives the two values at
  # which the p-value function equals 1 - level
  conf_levels <- seq(eps, 1 - eps, length.out = ceiling(n_values / 2))
  x_calc <- vapply(
    conf_levels,
    wilson_cicc_diff,
    estimate = estimate,
    n = n,
    FUN.VALUE = double(2L)
  )

  # The continuity correction leaves a gap around the estimate even for a
  # confidence level of (almost) 0. It is filled with values without p-values.
  n_gap <- 100L
  gap <- wilson_cicc_diff(estimate, n, conf_level = eps)
  values <- c(x_calc[1, ], x_calc[2, ], seq(gap[1], gap[2], length.out = n_gap))
  p_two <- c(1 - conf_levels, 1 - conf_levels, rep(NA, n_gap))

  # No confidence distribution and density for this one
  res_list <- list(cbind(values, NA, NA, p_two, p_two / 2, 1))

  # Confidence intervals

  conf_list <- list()

  if (!is.null(conf_level)) {
    limits <- vapply(
      two_sided_level(conf_level, alternative),
      wilson_cicc_diff,
      estimate = estimate,
      n = n,
      FUN.VALUE = double(2L)
    )

    conf_list[[1]] <- cdist_conf_matrix(
      conf_level = conf_level,
      limits = c(limits[1, ], limits[2, ]),
      i = 1
    )
  }

  # Counternulls: find the confidence level at which one limit equals the
  # null value; the counternull is then the other limit. Null values outside
  # of the computed range or inside the gap have no counternull.

  counternull_of <- function(null_value) {
    if (
      null_value > max(values) ||
        null_value < min(values) ||
        (null_value >= gap[1] && null_value <= gap[2])
    ) {
      return(NA_real_)
    }

    side <- if (null_value > est_diff) 2L else 1L

    null_conf <- uniroot(
      function(conf) {
        wilson_cicc_diff(estimate, n, conf_level = conf)[side] - null_value
      },
      lower = eps,
      upper = 1 - eps,
      tol = 1e-12
    )$root

    wilson_cicc_diff(estimate, n, conf_level = null_conf)[3L - side]
  }

  counternull_list <- list()

  if (!is.null(null_values)) {
    counternull_list[[1]] <- cdist_counternull_matrix(
      null_values = null_values,
      counternull = vapply(null_values, counternull_of, double(1L)),
      i = 1
    )
  }

  frames <- assemble_cdist_frames(
    res_list = res_list,
    conf_list = conf_list,
    counternull_list = counternull_list,
    conf_level = conf_level,
    null_values = null_values
  )

  c(frames, list(point_est = empty_point_est_frame(1)))
}
