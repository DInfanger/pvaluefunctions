# Internal helpers shared by the cdist_*() functions.

# Internal helper: empty point-estimates data frame for n estimates.
empty_point_est_frame <- function(n) {
  data.frame(
    est_mean = rep(NA_real_, n),
    est_median = rep(NA_real_, n),
    est_mode = rep(NA_real_, n),
    variable = seq_len(n)
  )
}

# Internal helper: assemble the shared result, confidence-interval and
# counternull data frames from the per-estimate matrices accumulated in the
# cdist_* functions.
assemble_cdist_frames <- function(
  res_list,
  conf_list,
  counternull_list,
  conf_level,
  null_values
) {
  res_frame <- as.data.frame(do.call(rbind, res_list))
  names(res_frame) <- c(
    "values",
    "conf_dist",
    "conf_dens",
    "p_two",
    "p_one",
    "variable"
  )

  if (!is.null(conf_level)) {
    conf_frame <- as.data.frame(do.call(rbind, conf_list))
    names(conf_frame) <- c("conf_level", "lwr", "upr", "variable")
  } else {
    conf_frame <- NULL
  }

  if (!is.null(null_values)) {
    counternull_frame <- as.data.frame(do.call(rbind, counternull_list))
    names(counternull_frame) <- c("null_value", "counternull", "variable")
  } else {
    counternull_frame <- NULL
  }

  list(
    res_frame = res_frame,
    conf_frame = conf_frame,
    counternull_frame = counternull_frame
  )
}

# Internal helper: median and mode point estimates for estimate i.
point_est_median_mode <- function(res_frame, i) {
  values <- res_frame$values[res_frame$variable == i]
  conf_dist <- res_frame$conf_dist[res_frame$variable == i]
  conf_dens <- res_frame$conf_dens[res_frame$variable == i]

  c(
    values[which.min(abs(conf_dist[-1] - 0.5)) + 1],
    values[which.max(conf_dens)]
  )
}

# Internal helper: probabilities of the lower and upper confidence limits, i.e.
# first all lower and then all upper probabilities.
conf_limit_probs <- function(conf_level, alternative) {
  switch(
    alternative,
    two_sided = c(1 - (conf_level + 1) / 2, (conf_level + 1) / 2),
    one_sided = c((1 - conf_level), conf_level)
  )
}

# Internal helper: result matrix for estimate i with the columns values,
# confidence distribution, confidence density, two-sided p-value, one-sided
# p-value and the index of the estimate.
cdist_res_matrix <- function(x, cdf, dens, i) {
  matrix(
    c(
      x,
      cdf,
      dens,
      1 - 2 * abs(cdf - (1 / 2)),
      (1 / 2) - abs(cdf - (1 / 2)),
      rep(i, length(x))
    ),
    ncol = 6
  )
}

# Internal helper: confidence interval matrix for estimate i. `limits` contains
# all lower limits followed by all upper limits (see conf_limit_probs()).
cdist_conf_matrix <- function(conf_level, limits, i) {
  limits <- matrix(limits, ncol = 2)

  matrix(
    c(conf_level, limits[, 1], limits[, 2], rep(i, length(conf_level))),
    ncol = 4
  )
}

# Internal helper: counternull matrix for estimate i.
cdist_counternull_matrix <- function(null_values, counternull, i) {
  matrix(c(null_values, counternull, rep(i, length(null_values))), ncol = 3)
}
