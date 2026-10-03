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

  # Point estimators: the distributions are symmetric around the estimate, so
  # median and mode equal the estimate. The mean of the t-distribution only
  # exists for df > 1 (df is NULL for the normal distribution).
  has_mean <- if (is.null(df)) rep(TRUE, length(estimate)) else df > 1

  point_est_frame <- data.frame(
    est_mean = ifelse(has_mean, estimate, NA_real_),
    est_median = estimate,
    est_mode = estimate,
    variable = seq_along(estimate)
  )

  c(frames, list(point_est = point_est_frame))
}
