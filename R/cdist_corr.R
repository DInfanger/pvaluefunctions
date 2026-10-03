# Confidence distributions for correlation coefficients.

cdist_corr <- function(
  estimate = NULL,
  stderr = NULL,
  n = NULL,
  n_values = NULL,
  conf_level = NULL,
  null_values = NULL,
  alternative = NULL
) {
  res_list <- list()

  conf_list <- list()

  counternull_list <- list()

  for (i in seq_along(estimate)) {
    limits <- c(-1, 1)

    x_calc <- c(
      estimate[i],
      null_values,
      seq(limits[1], limits[2], length.out = (n_values + 2))
    )
    x_calc <- x_calc[!(x_calc %in% c(-1, 1))]

    z_calc <- (1 / stderr[i]) * (atanh(estimate[i]) - atanh(x_calc))

    res_list[[length(res_list) + 1]] <- cdist_res_matrix(
      x = x_calc,
      cdf = 1 - pnorm(z_calc),
      dens = -dnorm(z_calc) * (-1 / (stderr[i] * (1 - x_calc^2))),
      i = i
    )

    # Confidence intervals

    if (!is.null(conf_level)) {
      quants_tmp <- conf_limit_probs(conf_level, alternative)

      conf_list[[length(conf_list) + 1]] <- cdist_conf_matrix(
        conf_level = conf_level,
        limits = tanh(qnorm(quants_tmp) * stderr[i] + atanh(estimate[i])),
        i = i
      )
    }

    # Counternulls

    if (!is.null(null_values)) {
      counternull_list[[length(counternull_list) + 1]] <-
        cdist_counternull_matrix(
          null_values = null_values,
          counternull = tanh(2 * atanh(estimate[i]) - atanh(null_values)),
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

  mean_fun <- function(x, stderr, estimate) {
    x *
      ((-dnorm((1 / stderr) * (atanh(estimate) - atanh(x))) *
        (-1 / (stderr * (1 - x^2)))))
  }

  for (i in seq_along(estimate)) {
    point_est_frame$est_mean[i] <- integrate(
      mean_fun,
      lower = -1,
      upper = 1,
      stderr = stderr[i],
      estimate = estimate[i],
      rel.tol = 1e-10
    )$value # Mean
    point_est_frame[i, c("est_median", "est_mode")] <-
      point_est_median_mode(frames$res_frame, i)
  }

  c(frames, list(point_est = point_est_frame))
}

cdist_corr_exact <- function(
  estimate = NULL,
  n = NULL,
  n_values = NULL,
  conf_level = NULL,
  null_values = NULL,
  alternative = NULL
) {
  # Define function for the exact density (see Taraldsen 2020)
  # The prefactor is computed on the log scale because gamma() overflows to Inf
  # for n > 172, which would result in NaN
  conf_dens_corr <- function(rho, r, n) {
    nu <- (n - 1)

    # For nu == 2 the exponent is 0, so the term is 1 even for |rho| == 1
    log_rho_term <- if (nu == 2) {
      0
    } else {
      ((nu - 2) / 2) * log1p(-rho^2)
    }

    log_prefactor <- log(nu) +
      log(nu - 1) +
      lgamma(nu - 1) -
      log(sqrt(2 * pi)) -
      lgamma(nu + (1 / 2)) +
      ((nu - 1) / 2) * log1p(-r^2) +
      log_rho_term +
      ((1 - 2 * nu) / 2) * log1p(-r * rho)

    exp(log_prefactor) *
      gsl::hyperg_2F1(
        -1 / 2,
        3 / 2,
        nu + (1 / 2),
        (1 + r * rho) / 2,
        strict = FALSE
      )
  }
  # Tail probabilities P(rho < z) (the cdf) or P(rho > z) by numerical
  # integration of the pdf. The integration error can push the result
  # slightly outside of [0, 1], which gives negative p-values.
  tail_prob <- function(z, r, n, upper = FALSE) {
    probs <- vapply(
      z,
      function(z) {
        integrate(
          conf_dens_corr,
          lower = if (upper) z else -1,
          upper = if (upper) 1 else z,
          r = r,
          n = n,
          subdivisions = 1000L
        )$value
      },
      double(1L)
    )
    pmin(pmax(probs, 0), 1)
  }

  # Function to find confidence intervals based on the cdf
  find_ci <- function(conf_level, r, n, alternative) {
    quants_tmp <- conf_limit_probs(conf_level, alternative)

    limit_at <- function(prob) {
      uniroot(
        function(z) tail_prob(z, r = r, n = n) - prob,
        interval = c(-1, 1),
        maxiter = 2000
      )$root
    }
    c(limit_at(quants_tmp[1]), limit_at(quants_tmp[2]))
  }

  # Function to find counternull values
  # The counternull is the value for which P(rho > x) = P(rho < null value).
  # The tail probabilities are integrated directly instead of using
  # 1 - cdf because this loses all precision for extreme tails (i.e. large n).
  find_counternulls <- function(null_values, r, n) {
    lower_tail_null <- tail_prob(null_values, r = r, n = n)

    if (lower_tail_null <= 1 / 2) {
      zero_fun <- function(x, target, r, n) {
        tail_prob(x, r = r, n = n, upper = TRUE) - target
      }
      target <- lower_tail_null
    } else {
      zero_fun <- function(x, target, r, n) {
        tail_prob(x, r = r, n = n) - target
      }
      target <- tail_prob(null_values, r = r, n = n, upper = TRUE)
    }

    uniroot(
      zero_fun,
      r = r,
      n = n,
      target = target,
      interval = c(-1, 1)
    )$root
  }

  res_list <- list()

  conf_list <- list()

  counternull_list <- list()

  for (i in seq_along(estimate)) {
    limits <- c(-1, 1)

    x_calc <- c(
      estimate[i],
      null_values,
      seq(limits[1], limits[2], length.out = (n_values + 2))
    )
    x_calc <- x_calc[!(x_calc %in% c(-1, 1))]

    res_list[[length(res_list) + 1]] <- cdist_res_matrix(
      x = x_calc,
      cdf = tail_prob(x_calc, r = estimate[i], n = n[i]),
      dens = conf_dens_corr(x_calc, r = estimate[i], n = n[i]),
      i = i
    )

    # Confidence intervals

    if (!is.null(conf_level)) {
      limits_tmp <- vapply(
        conf_level,
        FUN = find_ci,
        r = estimate[i],
        n = n[i],
        alternative = alternative,
        FUN.VALUE = double(2L)
      )

      conf_list[[length(conf_list) + 1]] <- cdist_conf_matrix(
        conf_level = conf_level,
        limits = c(limits_tmp[1, ], limits_tmp[2, ]),
        i = i
      )
    }

    # Counternulls

    if (!is.null(null_values)) {
      counternulls_tmp <- vapply(
        null_values,
        find_counternulls,
        r = estimate[i],
        n = n[i],
        FUN.VALUE = double(1L)
      )

      counternull_list[[length(counternull_list) + 1]] <-
        cdist_counternull_matrix(
          null_values = null_values,
          counternull = counternulls_tmp,
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

  # The mean is obtained by numerical integration
  mean_fun <- function(r, n) {
    int_fun <- function(rho, r, n) {
      rho * conf_dens_corr(rho = rho, r = r, n = n)
    }
    integrate(
      int_fun,
      lower = -1,
      upper = 1,
      r = r,
      n = n,
      rel.tol = 1e-10
    )$value
  }

  # The median is found by numerical root-finding using the cdf
  median_fun <- function(r, n) {
    zero_fun <- function(x, r, n) {
      tail_prob(x, r = r, n = n) - (1 / 2)
    }
    uniroot(zero_fun, r = r, n = n, interval = c(-1, 1))$root
  }

  # Find the root of the derivative of the pdf numerically
  mode_fun <- function(r, n) {
    # `drop_positive_factor` removes the strictly positive factor that
    # multiplies the derivative. This does not change its roots but avoids
    # underflow to exactly 0 at the ends of the search interval for large n.
    pdf_deriv <- function(rho, r, n, drop_positive_factor = FALSE) {
      nu <- (n - 1)

      # The powers and the gamma ratios are combined on the log scale
      # because gamma() overflows to Inf for large n
      log_powers <- (1 / 2 * (nu - 1)) *
        log1p(-r^2) +
        (-1 / 2 - nu) * log1p(-r * rho) +
        (1 / 2 * (nu - 4)) * log1p(-rho^2)
      gamma_ratio_1 <- exp(lgamma(nu - 1) - lgamma(1 / 2 + nu))
      gamma_ratio_2 <- exp(lgamma(nu - 1) - lgamma(3 / 2 + nu))
      positive_factor <- if (drop_positive_factor) 1 else exp(log_powers)

      -1 /
        (8 * sqrt(2 * pi)) *
        positive_factor *
        (nu - 1) *
        nu *
        (4 *
          (2 * (nu - 2) * rho + r * (1 - 2 * nu + 3 * rho^2)) *
          (gsl::hyperg_2F1(
            -1 / 2,
            3 / 2,
            1 / 2 + nu,
            1 / 2 * (1 + r * rho),
            strict = FALSE
          )) *
          gamma_ratio_1 +
          3 *
            r *
            (r * rho - 1) *
            (rho^2 - 1) *
            (gsl::hyperg_2F1(
              1 / 2,
              5 / 2,
              3 / 2 + nu,
              1 / 2 * (1 + rho * r),
              strict = FALSE
            ) *
              gamma_ratio_2))
    }
    bounds <- c(-0.99999, 0.99999)
    deriv_at_bounds <- pdf_deriv(bounds, r = r, n = n)
    drop_factor <- !all(is.finite(deriv_at_bounds)) || any(deriv_at_bounds == 0)

    uniroot(
      pdf_deriv,
      r = r,
      n = n,
      drop_positive_factor = drop_factor,
      lower = bounds[1],
      upper = bounds[2]
    )$root
  }

  for (i in seq_along(estimate)) {
    point_est_frame$est_mean[i] <- mean_fun(r = estimate[i], n = n[i])
    point_est_frame$est_median[i] <- median_fun(r = estimate[i], n = n[i])
    point_est_frame$est_mode[i] <- mode_fun(r = estimate[i], n = n[i])
  }

  c(frames, list(point_est = point_est_frame))
}
