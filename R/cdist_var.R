# Confidence distributions for variances.

cdist_var <- function(
  estimate = NULL,
  n = NULL,
  n_values = NULL,
  conf_level = NULL,
  null_values = NULL,
  alternative = NULL
) {
  eps <- 1e-10

  df <- (n - 1)

  res_list <- list()

  conf_list <- list()

  counternull_list <- list()

  for (i in seq_along(estimate)) {
    limits <- c(
      estimate[i] * df[i] / qchisq(eps, df = df[i], lower.tail = FALSE),
      estimate[i] * df[i] / qchisq(1 - eps, df = df[i], lower.tail = FALSE)
    )

    x_calc <- c(
      estimate[i],
      null_values,
      seq(limits[1], limits[2], length.out = n_values)
    )

    chisq_calc <- df[i] * estimate[i] / x_calc

    res_list[[length(res_list) + 1]] <- cdist_res_matrix(
      x = x_calc,
      cdf = 1 - pchisq(chisq_calc, df = df[i]),
      dens = -dchisq(chisq_calc, df = df[i]) *
        (-((df[i] * estimate[i]) / x_calc^2)),
      i = i
    )

    # Confidence intervals

    if (!is.null(conf_level)) {
      quants_tmp <- conf_limit_probs(conf_level, alternative)

      conf_list[[length(conf_list) + 1]] <- cdist_conf_matrix(
        conf_level = conf_level,
        limits = (estimate[i] * df[i]) /
          qchisq(quants_tmp, df = df[i], lower.tail = FALSE),
        i = i
      )
    }

    # Counternulls

    if (!is.null(null_values)) {
      counter_tmp <- 1 -
        (1 - pchisq(df[i] * estimate[i] / null_values, df = df[i]))

      counternull_list[[length(counternull_list) + 1]] <-
        cdist_counternull_matrix(
          null_values = null_values,
          counternull = estimate[i] *
            df[i] /
            qchisq(counter_tmp, df = df[i], lower.tail = FALSE),
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

  # The mean of the scaled inverse chi-squared distribution only exists for
  # df > 2
  point_est_frame$est_mean <- ifelse(
    df > 2,
    estimate * df / (df - 2),
    NA_real_
  )
  point_est_frame$est_median <- estimate *
    df /
    (2 * stats::qgamma(0.5, df / 2, lower.tail = FALSE))
  point_est_frame$est_mode <- exp((log(estimate) + log(df)) - (log(2 + df)))

  c(frames, list(point_est = point_est_frame))
}
