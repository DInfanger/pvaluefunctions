# Confidence distributions for location-scale families: estimates that follow a
# t-distribution (linear regression, t-tests, ...) or a normal distribution
# (logistic regression, Cox regression, ...) after standardisation.

# Internal helper: quantile, distribution and density function of the
# standardised distribution for estimate i. `df` is NULL for the normal
# distribution.
location_scale_funs <- function(df, i) {
  if (is.null(df)) {
    list(q = qnorm, p = pnorm, d = dnorm)
  } else {
    list(
      q = function(x) qt(x, df = df[i]),
      p = function(x) pt(x, df = df[i]),
      d = function(x) dt(x, df = df[i])
    )
  }
}

cdist_location_scale <- function(
  estimate,
  stderr,
  df,
  n_values,
  conf_level,
  null_values,
  alternative
) {
  eps <- 1e-10

  res_list <- list()
  conf_list <- list()
  counternull_list <- list()

  for (i in seq_along(estimate)) {
    funs <- location_scale_funs(df, i)

    limits <- c(
      funs$q(eps) * stderr[i] + estimate[i],
      funs$q(1 - eps) * stderr[i] + estimate[i]
    )

    x_calc <- c(
      estimate[i],
      null_values,
      seq(limits[1], limits[2], length.out = n_values)
    )

    z_calc <- (x_calc - estimate[i]) / stderr[i]

    res_list[[length(res_list) + 1]] <- cdist_res_matrix(
      x = x_calc,
      cdf = funs$p(z_calc),
      dens = funs$d(z_calc) * (1 / stderr[i]),
      i = i
    )

    # Confidence intervals

    if (!is.null(conf_level)) {
      quants_tmp <- conf_limit_probs(conf_level, alternative)

      conf_list[[length(conf_list) + 1]] <- cdist_conf_matrix(
        conf_level = conf_level,
        limits = funs$q(quants_tmp) * stderr[i] + estimate[i],
        i = i
      )
    }

    # Counternulls

    if (!is.null(null_values)) {
      counternull_list[[length(counternull_list) + 1]] <-
        cdist_counternull_matrix(
          null_values = null_values,
          counternull = 2 * estimate[i] - null_values,
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

  mean_fun <- function(x, estimate, stderr, dens_fun) {
    x * dens_fun((x - estimate) / stderr) * (1 / stderr)
  }

  for (i in seq_along(estimate)) {
    point_est_frame$est_mean[i] <- integrate(
      mean_fun,
      lower = -Inf,
      upper = Inf,
      estimate = estimate[i],
      stderr = stderr[i],
      dens_fun = location_scale_funs(df, i)$d,
      rel.tol = 1e-10
    )$value # Mean
    point_est_frame[i, c("est_median", "est_mode")] <-
      point_est_median_mode(frames$res_frame, i)
  }

  list(
    res_frame = frames$res_frame,
    conf_frame = frames$conf_frame,
    counternull_frame = frames$counternull_frame,
    point_est = point_est_frame
  )
}
