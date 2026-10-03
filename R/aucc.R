# Area under the confidence curve (AUCC), see Berrar (2017) Mach Learn
# 106:911-949.

# Internal helper: trapezoidal integration of y over x after removing missing
# values and sorting by x.
trapz_sorted <- function(x, y) {
  not_missing <- which(!is.na(x) & !is.na(y))
  x <- x[not_missing]
  y <- y[not_missing]
  ordered <- order(x, decreasing = FALSE)

  x <- x[ordered]
  y <- y[ordered]
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}

# Internal helper: AUCC of the two-sided p-value function of every estimate.
# If null values are given, the proportion of the AUCC above each null value
# is calculated as well.
compute_aucc <- function(res_frame, n_estimates, null_values) {
  if (is.null(null_values)) {
    aucc_frame <- data.frame(
      variable = seq_len(n_estimates),
      aucc = NA
    )
  } else {
    aucc_frame <- data.frame(
      variable = rep(seq_len(n_estimates), each = length(null_values)),
      aucc = NA,
      null = rep(null_values, times = n_estimates),
      p_above_null = NA
    )
  }

  for (i in seq_len(n_estimates)) {
    x <- res_frame$values[res_frame$variable == i]
    y <- res_frame$p_two[res_frame$variable == i]

    aucc_frame$aucc[aucc_frame$variable == i] <- trapz_sorted(x, y)

    for (j in seq_along(null_values)) {
      above_null <- which(x > null_values[j])
      aucc_above_null <- trapz_sorted(x[above_null], y[above_null])

      is_row <- aucc_frame$variable == i & aucc_frame$null == null_values[j]
      aucc_frame$p_above_null[is_row] <- aucc_above_null /
        aucc_frame$aucc[is_row]
    }
  }

  aucc_frame
}
