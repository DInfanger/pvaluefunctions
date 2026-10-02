# Custom y-axis transformations (logarithmic below a cutoff, linear above).

# Define new mixed scale: log for p <= interval_low, else linear

magnify_trans_log <- function(
  interval_low = 0.05,
  interval_high = 1,
  reducer = 0.05,
  reducer2 = 8
) {
  trans <- Vectorize(function(
    x,
    i_low = interval_low,
    i_high = interval_high,
    r = reducer,
    r2 = reducer2
  ) {
    if (is.na(x) || (x >= i_low && x <= i_high)) {
      x
    } else if (x < i_low && !is.na(x)) {
      (log10(x / r) / r2 + i_low)
    } else {
      log10((x - i_high) / r + i_high) / r2
    }
  })

  inv <- Vectorize(function(
    x,
    i_low = interval_low,
    i_high = interval_high,
    r = reducer,
    r2 = reducer2
  ) {
    if (is.na(x) || (x >= i_low && x <= i_high)) {
      x
    } else if (x < i_low && !is.na(x)) {
      10^(-(i_low - x) * r2) * r
    } else {
      i_high + 10^(x * r2) * r - i_high * r
    }
  })

  scales::new_transform(
    name = "customlog",
    transform = trans,
    inverse = inv,
    domain = c(1e-16, Inf)
  )
}

magnify_trans_log_rev <- function(
  interval_low = 0.05,
  interval_high = 1,
  reducer = 0.05,
  reducer2 = 8
) {
  trans <- Vectorize(function(
    x,
    i_low = interval_low,
    i_high = interval_high,
    r = reducer,
    r2 = reducer2
  ) {
    -if (is.na(x) || (x >= i_low && x <= i_high)) {
      x
    } else if (x < i_low && !is.na(x)) {
      (log10(x / r) / r2 + i_low)
    } else {
      log10((x - i_high) / r + i_high) / r2 + i_high
    }
  })

  inv <- Vectorize(function(
    x,
    i_low = interval_low,
    i_high = interval_high,
    r = reducer,
    r2 = reducer2
  ) {
    if (is.na(x) || (-x >= i_low && -x <= i_high)) {
      -x
    } else if (-x < i_low && !is.na(x)) {
      (10^(-(i_low + x) * r2) * r)
    } else {
      i_high + 10^(-r2 * (i_high + x)) * r - i_high * r
    }
  })

  scales::new_transform(
    name = "customlog_rev",
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
