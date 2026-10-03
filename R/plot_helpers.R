# Custom y-axis transformations (logarithmic below a cutoff, linear above).

# Define new mixed scale: log for p <= interval_low, else linear

# With rev = TRUE the transformation is mirrored for an inverted y-axis.
magnify_trans <- function(
  interval_low = 0.05,
  interval_high = 1,
  reducer = 0.05,
  reducer2 = 8,
  rev = FALSE
) {
  sgn <- if (rev) -1 else 1
  offset <- if (rev) interval_high else 0

  # Apply `below` to the elements < interval_low and `above` to the elements
  # > interval_high; all other elements (and NAs) are left unchanged.
  piecewise <- function(x, below, above) {
    low <- which(x < interval_low)
    high <- which(x > interval_high)
    x[low] <- below(x[low])
    x[high] <- above(x[high])
    x
  }

  trans <- function(x) {
    sgn *
      piecewise(
        x,
        function(x) log10(x / reducer) / reducer2 + interval_low,
        function(x) {
          log10((x - interval_high) / reducer + interval_high) /
            reducer2 +
            offset
        }
      )
  }

  inv <- function(x) {
    piecewise(
      sgn * x,
      function(x) 10^(-(interval_low - x) * reducer2) * reducer,
      function(x) {
        interval_high +
          10^((x - offset) * reducer2) * reducer -
          interval_high * reducer
      }
    )
  }

  scales::new_transform(
    name = if (rev) "customlog_rev" else "customlog",
    transform = trans,
    inverse = inv,
    domain = c(1e-16, Inf)
  )
}

# Internal helper: data frame with the positions and labels of the confidence
# levels that are displayed in the plot.
make_text_frame <- function(conf_level, alternative, trans) {
  theor_values <- switch(
    alternative,
    two_sided = rep(
      ifelse(trans %in% "exp", 0, -Inf),
      each = length(conf_level)
    ),
    one_sided = rep(Inf, each = length(conf_level))
  )

  p_value <- switch(
    alternative,
    two_sided = round((1 - conf_level), 10),
    one_sided = round((2 - 2 * conf_level), 10)
  )

  data.frame(
    label = (1 - conf_level),
    theor_values = theor_values,
    p_value = p_value
  )
}
