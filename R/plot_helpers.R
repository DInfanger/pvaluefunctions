# Custom y-axis transformations (logarithmic below a cutoff, linear above).

# Define new mixed scale: log for p <= interval_low, else linear

magnify_trans_log <- function(
  interval_low = 0.05,
  interval_high = 1,
  reducer = 0.05,
  reducer2 = 8
) {
  trans <- Vectorize(function(x) {
    if (is.na(x) || (x >= interval_low && x <= interval_high)) {
      x
    } else if (x < interval_low) {
      (log10(x / reducer) / reducer2 + interval_low)
    } else {
      log10((x - interval_high) / reducer + interval_high) / reducer2
    }
  })

  inv <- Vectorize(function(x) {
    if (is.na(x) || (x >= interval_low && x <= interval_high)) {
      x
    } else if (x < interval_low) {
      10^(-(interval_low - x) * reducer2) * reducer
    } else {
      interval_high + 10^(x * reducer2) * reducer - interval_high * reducer
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
  trans <- Vectorize(function(x) {
    -if (is.na(x) || (x >= interval_low && x <= interval_high)) {
      x
    } else if (x < interval_low) {
      (log10(x / reducer) / reducer2 + interval_low)
    } else {
      log10((x - interval_high) / reducer + interval_high) / reducer2 + interval_high
    }
  })

  inv <- Vectorize(function(x) {
    if (is.na(x) || (-x >= interval_low && -x <= interval_high)) {
      -x
    } else if (-x < interval_low) {
      (10^(-(interval_low + x) * reducer2) * reducer)
    } else {
      interval_high + 10^(-reducer2 * (interval_high + x)) * reducer - interval_high * reducer
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
