# Shared helpers for the tests.

# conf_dist() prints the plot if `plot = TRUE`. A null device keeps the tests
# from creating Rplots.pdf or opening a window.
conf_dist_quiet <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  conf_dist(...)
}

# Data only: no plot and fewer grid points than the default to keep the tests
# fast. With the default `plot_p_limit`, values with very small p-values are set
# to missing in the returned data frame, which would get in the way of
# comparisons with closed-form results.
cd <- function(..., n_values = 200L, plot_p_limit = 0) {
  conf_dist(..., n_values = n_values, plot_p_limit = plot_p_limit, plot = FALSE)
}

# Result for one estimate from the returned res_frame
rows_of <- function(res, i = 1L) {
  res$res_frame[res$res_frame$variable == i, ]
}
