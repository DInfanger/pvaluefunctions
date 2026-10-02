# Confidence distributions for proportions and differences of proportions.

wilson_ci <- function(
  estimate,
  n,
  conf_level,
  alternative
) {
  z <- qnorm((conf_level + 1) / 2)

  p1 <- estimate + (1 / 2) * z^2 / n
  p2 <- z * sqrt((estimate * (1 - estimate) + (1 / 4) * z^2 / n) / n)
  p3 <- 1 + z^2 / n

  c((p1 - p2) / p3, (p1 + p2) / p3)
}

wilson_cicc <- function(
  estimate,
  n,
  conf_level,
  alternative
) {
  z <- qnorm((conf_level + 1) / 2)
  x <- round(estimate * n) # To get number of successes/failures
  estimate_compl <- (1 - estimate) # Complement of estimate

  lower <- max(
    0,
    (2 *
      x +
      z^2 -
      1 -
      z * sqrt(z^2 - 2 - 1 / n + 4 * estimate * (n * estimate_compl + 1))) /
      (2 * (n + z^2)),
    na.rm = TRUE
  )
  upper <- min(
    1,
    (2 *
      x +
      z^2 +
      1 +
      z * sqrt(z^2 + 2 - 1 / n + 4 * estimate * (n * estimate_compl - 1))) /
      (2 * (n + z^2)),
    na.rm = TRUE
  )

  c(lower, upper)
}

wilson_cicc_diff <- function(
  estimate,
  n,
  conf_level,
  alternative
) {
  est_diff <- (estimate[1] - estimate[2])

  res1 <- wilson_cicc(
    estimate = estimate[1],
    n = n[1],
    conf_level = conf_level,
    alternative = alternative
  )

  res2 <- wilson_cicc(
    estimate = estimate[2],
    n = n[2],
    conf_level = conf_level,
    alternative = alternative
  )

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
  # Auxilliary functions

  cdf_fun <- function(x, n, p) {
    x <- as.complex(x)
    n <- as.complex(n)
    p <- as.complex(p)

    -Re((1i * sqrt(n) * (p - x)) / (sqrt(x - 1) * sqrt(x)))
  }

  deriv_fun <- function(x, n, p) {
    x <- as.complex(x)
    n <- as.complex(n)
    p <- as.complex(p)

    Re(
      (1i * sqrt(n) * (-x + p * (-1 + 2 * x))) /
        (2 * (-1 + x)^(3 / 2) * x^(3 / 2))
    )
  }

  eps <- 1e-10

  res_list <- list()

  conf_list <- list()

  counternull_list <- list()

  for (i in seq_along(estimate)) {
    limits <- wilson_ci(
      estimate = estimate[i],
      n = n[i],
      conf_level = (1 - eps),
      alternative = alternative
    )

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
      # A one-sided level corresponds to a two-sided interval of level 2 * l - 1
      conf_tmp <- switch(
        alternative,
        two_sided = conf_level,
        one_sided = 2 * conf_level - 1
      )

      conf_list[[length(conf_list) + 1]] <- cdist_conf_matrix(
        conf_level = conf_level,
        limits = wilson_ci(
          estimate = estimate[i],
          n = n[i],
          conf_level = conf_tmp,
          alternative = alternative
        ),
        i = i
      )
    }

    # Counternulls

    if (!is.null(null_values)) {
      counter_tmp <- 1 -
        (pnorm(cdf_fun(x = null_values, n = n[i], p = estimate[i])))

      counternull_list[[length(counternull_list) + 1]] <-
        cdist_counternull_matrix(
          null_values = null_values,
          counternull = vapply(
            counter_tmp,
            function(x, cdf, q) {
              x[which.min(abs(cdf - q))]
            },
            x = x_calc,
            cdf = cdf_calc,
            FUN.VALUE = double(1L)
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

  res_frame <- frames$res_frame
  conf_frame <- frames$conf_frame
  counternull_frame <- frames$counternull_frame

  # Point estimators

  point_est_frame <- empty_point_est_frame(length(estimate))

  mean_fun <- function(x, estimate, n) {
    x *
      dnorm(cdf_fun(x = x, n = n, p = estimate)) *
      deriv_fun(x = x, n = n, p = estimate)
  }

  for (i in seq_along(estimate)) {
    point_est_frame$est_mean[i] <- integrate(
      mean_fun,
      lower = 0,
      upper = 1,
      n = n[i],
      estimate = estimate[i],
      rel.tol = 1e-10
    )$value # Mean
    point_est_frame[i, c("est_median", "est_mode")] <-
      point_est_median_mode(res_frame, i)
  }

  list(
    res_frame = res_frame,
    conf_frame = conf_frame,
    counternull_frame = counternull_frame,
    point_est = point_est_frame
  )
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

  res_list <- list()

  conf_list <- list()

  counternull_list <- list()

  conf_levels <- seq(eps, 1 - eps, length.out = ceiling(n_values / 2))
  x_calc <- vapply(
    conf_levels,
    wilson_cicc_diff,
    estimate = estimate,
    n = n,
    alternative = alternative,
    FUN.VALUE = double(2L)
  )

  val_min <- wilson_cicc_diff(estimate, n, conf_level = eps)
  val_between <- seq(min(val_min), max(val_min), length.out = 100)

  res_mat_tmp <- matrix(NA, nrow = length(x_calc) + 100, ncol = 6)

  res_mat_tmp[, 1] <- c(x_calc[1, ], x_calc[2, ], val_between)
  is.na(res_mat_tmp[, 2]) <- TRUE # No confidence distribution for this one
  is.na(res_mat_tmp[, 3]) <- TRUE # No confidence density for this one
  res_mat_tmp[, 4] <- c(1 - conf_levels, 1 - conf_levels, rep(NA, 100))
  res_mat_tmp[, 5] <- c(
    (1 - conf_levels) / 2,
    (1 - conf_levels) / 2,
    rep(NA, 100)
  )
  res_mat_tmp[, 6] <- rep(1, times = (length(x_calc) + 100))

  res_list[[length(res_list) + 1]] <- res_mat_tmp

  # Confidence intervals

  if (!is.null(conf_level)) {
    conf_tmp <- switch(
      alternative,
      two_sided = conf_level,
      one_sided = 2 * conf_level - 1
    )

    conf_mat_tmp <- matrix(NA, ncol = 4, nrow = length(conf_level))

    limits_tmp <- matrix(NA, ncol = 2, nrow = length(conf_tmp))

    for (j in seq_along(conf_tmp)) {
      limits_tmp[j, ] <- wilson_cicc_diff(
        estimate = estimate,
        n = n,
        conf_level = conf_tmp[j],
        alternative = alternative
      )
    }

    conf_mat_tmp[, 1] <- conf_level
    conf_mat_tmp[, 2] <- limits_tmp[, 1]
    conf_mat_tmp[, 3] <- limits_tmp[, 2]
    conf_mat_tmp[, 4] <- rep(1, length(conf_level))

    conf_list[[length(conf_list) + 1]] <- conf_mat_tmp
  }

  # Counternulls

  if (!is.null(null_values)) {
    cnull_tmp <- rep(NA, length(null_values))

    tmp_fun_up <- function(conf_level, estimate, n, null_values) {
      wilson_cicc_diff(estimate = estimate, n = n, conf_level = conf_level)[2] -
        null_values
    }

    tmp_fun_low <- function(conf_level, estimate, n, null_values) {
      wilson_cicc_diff(estimate = estimate, n = n, conf_level = conf_level)[1] -
        null_values
    }

    for (j in seq_along(null_values)) {
      if (
        (null_values[j] > max(res_mat_tmp[, 1])) ||
          (null_values[j] < min(res_mat_tmp[, 1]))
      ) {
        is.na(cnull_tmp[j]) <- TRUE # Set values in the "gap" to missing for plotting
      } else {
        if (null_values[j] > -diff(estimate)) {
          null_conf <- uniroot(
            tmp_fun_up,
            lower = 1e-15,
            upper = 1 - 1e-15,
            null_values = null_values[j],
            estimate = estimate,
            n = n
          )$root

          cnull_tmp[j] <- wilson_cicc_diff(
            estimate = estimate,
            n = n,
            conf_level = null_conf
          )[1]
        } else if (null_values[j] < -diff(estimate)) {
          null_conf <- uniroot(
            tmp_fun_low,
            lower = 1e-15,
            upper = 1 - 1e-15,
            null_values = null_values[j],
            estimate = estimate,
            n = n
          )$root

          cnull_tmp[j] <- wilson_cicc_diff(
            estimate = estimate,
            n = n,
            conf_level = null_conf
          )[2]
        }
      }
    }

    counternull_mat_tmp <- matrix(NA, ncol = 3, nrow = length(null_values))
    counternull_mat_tmp[, 1] <- null_values
    counternull_mat_tmp[, 2] <- cnull_tmp
    counternull_mat_tmp[, 3] <- rep(1, length(null_values))
    counternull_list[[length(counternull_list) + 1]] <- counternull_mat_tmp
  }

  frames <- assemble_cdist_frames(
    res_list = res_list,
    conf_list = conf_list,
    counternull_list = counternull_list,
    conf_level = conf_level,
    null_values = null_values
  )

  res_frame <- frames$res_frame
  conf_frame <- frames$conf_frame
  counternull_frame <- frames$counternull_frame

  # Point estimators

  point_est_frame <- empty_point_est_frame(1)

  list(
    res_frame = res_frame,
    conf_frame = conf_frame,
    counternull_frame = counternull_frame,
    point_est = point_est_frame
  )
}
