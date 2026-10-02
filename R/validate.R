# Input validation for conf_dist().
#
# `validate_args()` only checks the arguments (it signals errors and warnings
# but never modifies anything). Normalising the arguments (e.g. dropping
# confidence levels that are too low, sorting `xlim`, ...) is done in
# conf_dist() after the validation.

valid_types <- c(
  "ttest",
  "linreg",
  "gammareg",
  "general_t",
  "logreg",
  "poisreg",
  "coxreg",
  "general_z",
  "pearson",
  "spearman",
  "kendall",
  "var",
  "prop",
  "propdiff"
)

# Scalar TRUE or FALSE
check_flag <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be either {.code TRUE} or {.code FALSE}.",
        x = "You supplied {.cls {class(x)}} of length {length(x)}."
      ),
      call = call
    )
  }

  invisible(x)
}

# Single finite number within [lower, upper] (or (lower, upper) if open)
check_number <- function(
  x,
  lower = -Inf,
  upper = Inf,
  lower_open = FALSE,
  upper_open = FALSE,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)

  if (ok) {
    ok <- (if (lower_open) x > lower else x >= lower) &&
      (if (upper_open) x < upper else x <= upper)
  }

  if (!ok) {
    range_text <- paste0(
      if (lower_open) "(" else "[",
      lower,
      ", ",
      upper,
      if (upper_open) ")" else "]"
    )

    cli::cli_abort(
      c(
        "{.arg {arg}} must be a single number in {range_text}.",
        x = "You supplied {.val {x}}."
      ),
      call = call
    )
  }

  invisible(x)
}

# Numeric vector without missing values (and without infinite values unless
# `finite = FALSE`). With `lower`, all values must be larger than that bound.
check_numeric <- function(
  x,
  finite = TRUE,
  lower = NULL,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!is.numeric(x) || length(x) == 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric vector.",
        x = "You supplied {.cls {class(x)}} of length {length(x)}."
      ),
      call = call
    )
  }

  if (anyNA(x) || (finite && !all(is.finite(x)))) {
    cli::cli_abort(
      "{.arg {arg}} must not contain missing{if (finite) ' or infinite'} values.",
      call = call
    )
  }

  if (!is.null(lower) && any(x <= lower)) {
    cli::cli_abort(
      c(
        "All values of {.arg {arg}} must be larger than {lower}.",
        x = "You supplied {.val {x}}."
      ),
      call = call
    )
  }

  invisible(x)
}

# Optional single positive whole number (e.g. `nrow` and `ncol`)
check_count <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (is.null(x)) {
    return(invisible(x))
  }

  if (
    !(is.numeric(x) &&
      length(x) == 1L &&
      is.finite(x) &&
      x >= 1 &&
      x == round(x))
  ) {
    cli::cli_abort(
      "{.arg {arg}} must be {.code NULL} or a single whole number >= 1.",
      call = call
    )
  }

  invisible(x)
}

validate_args <- function(
  estimate,
  n,
  df,
  stderr,
  tstat,
  type,
  plot_type,
  n_values,
  est_names,
  conf_level,
  null_values,
  alternative,
  log_yaxis,
  cut_logyaxis,
  xlab,
  xlim,
  together,
  plot_legend,
  same_color,
  nrow,
  ncol,
  plot_p_limit,
  plot_counternull,
  inverted,
  plot,
  call = rlang::caller_env()
) {
  #-----------------------------------------------------------------------------
  # Estimate and type
  #-----------------------------------------------------------------------------

  if (is.null(estimate)) {
    cli::cli_abort("Please provide an estimate.", call = call)
  }

  if (length(type) == 0L) {
    cli::cli_abort("Please provide the type of the estimate(s).", call = call)
  }

  if (
    !is.character(type) ||
      length(type) != 1L ||
      !tolower(type) %in% valid_types
  ) {
    cli::cli_abort(
      c(
        "{.arg type} must be one of {.val {valid_types}}.",
        x = "You supplied {.val {type}}."
      ),
      call = call
    )
  }

  type <- tolower(type)

  check_numeric(estimate, call = call)

  #-----------------------------------------------------------------------------
  # Options that are not specific to a type
  #-----------------------------------------------------------------------------

  check_flag(log_yaxis, call = call)
  check_flag(together, call = call)
  check_flag(plot_legend, call = call)
  check_flag(same_color, call = call)
  check_flag(plot_counternull, call = call)
  check_flag(inverted, call = call)
  check_flag(plot, call = call)

  check_number(n_values, lower = 2, call = call)
  check_number(
    plot_p_limit,
    lower = 0,
    upper = 1,
    upper_open = TRUE,
    call = call
  )
  check_number(
    cut_logyaxis,
    lower = 0,
    upper = 1,
    lower_open = TRUE,
    call = call
  )
  check_count(nrow, call = call)
  check_count(ncol, call = call)

  plot_p_limit <- round(plot_p_limit, 10)

  if (plot_p_limit == 0 && isTRUE(log_yaxis)) {
    cli::cli_abort("Cannot plot 0 on logarithmic axis.", call = call)
  }

  if (plot_p_limit >= 0.5 && alternative %in% "one_sided") {
    cli::cli_abort(
      "Plot limit must be below 0.5 for one-sided hypotheses.",
      call = call
    )
  }

  if (length(conf_level) > 0L) {
    check_numeric(conf_level, call = call)

    # Levels that are too low for the type of the alternative are dropped
    # in conf_dist(); levels of 1 or more can never be valid.
    if (any(conf_level >= 1)) {
      cli::cli_abort("All confidence levels must lie between 0 and 1.", call = call)
    }
  }

  if (!is.null(null_values)) {
    check_numeric(null_values, call = call)
  }

  if (!is.null(xlab) && (length(xlab) != 1L)) {
    cli::cli_abort("Length of x-axis label must be 1.", call = call)
  }

  if (!is.null(xlim)) {
    if (length(xlim) != 2L) {
      cli::cli_abort("Please provide two limits for the x-axis.", call = call)
    }

    if (!is.numeric(xlim) || anyNA(xlim) || any(!is.finite(xlim))) {
      cli::cli_abort(
        "Missing or infinite values are not allowed for x-axis limits (xlim). Please provide exactly two finite x-axis limits.",
        call = call
      )
    }
  }

  #-----------------------------------------------------------------------------
  # Required arguments by type
  #-----------------------------------------------------------------------------

  if (type %in% "ttest" && (is.null(tstat) || is.null(df))) {
    cli::cli_abort(
      "Please provide the t-statistic and the degrees of freedom of the t-test.",
      call = call
    )
  }

  if (
    type %in%
      c("linreg", "gammareg", "general_t") &&
      (is.null(df) || is.null(stderr))
  ) {
    cli::cli_abort(
      "Please provide the (residual) degrees of freedom and the standard error of the estimates.",
      call = call
    )
  }

  if (
    type %in% c("logreg", "poisreg", "coxreg", "general_z") && is.null(stderr)
  ) {
    cli::cli_abort(
      "Please provide the standard error of the estimates.",
      call = call
    )
  }

  if (type %in% c("pearson", "spearman", "kendall", "prop") && is.null(n)) {
    cli::cli_abort(
      "Please provide the sample size for correlations and proportions.",
      call = call
    )
  }

  if (type %in% "var" && is.null(n)) {
    cli::cli_abort(
      "Sample size must be given for variance estimates.",
      call = call
    )
  }

  #-----------------------------------------------------------------------------
  # Types of the arguments by type
  #-----------------------------------------------------------------------------

  if (type %in% c("ttest", "linreg", "gammareg", "general_t")) {
    check_numeric(df, finite = FALSE, lower = 0, call = call)
  }

  if (type %in% "ttest") {
    check_numeric(tstat, call = call)
  }

  if (
    type %in%
      c(
        "linreg",
        "gammareg",
        "general_t",
        "logreg",
        "poisreg",
        "coxreg",
        "general_z"
      )
  ) {
    check_numeric(stderr, lower = 0, call = call)
  }

  if (
    type %in% c("pearson", "spearman", "kendall", "var", "prop", "propdiff")
  ) {
    check_numeric(n, lower = 0, call = call)
  }

  #-----------------------------------------------------------------------------
  # Lengths
  #-----------------------------------------------------------------------------

  if (type %in% "ttest" && !same_length(df, estimate, tstat)) {
    cli::cli_abort(
      "Degrees of freedom (df) and t-statistics (tstat) must be the same length as estimates.",
      call = call
    )
  }

  if (
    type %in%
      c("linreg", "gammareg", "general_t") &&
      !same_length(stderr, df, estimate)
  ) {
    cli::cli_abort(
      "Standard errors (stderr) and degrees of freedom (df) must be the same length as estimates.",
      call = call
    )
  }

  if (
    type %in%
      c("coxreg", "logreg", "poisreg", "general_z") &&
      !same_length(stderr, estimate)
  ) {
    cli::cli_abort(
      "Standard errors (stderr) must be the same length as estimates.",
      call = call
    )
  }

  if (
    type %in%
      c("pearson", "spearman", "kendall", "var", "prop") &&
      !same_length(n, estimate)
  ) {
    cli::cli_abort(
      "Sample sizes (n) must be the same length as estimates.",
      call = call
    )
  }

  if (type %in% "propdiff" && ((length(estimate) != 2L) || (length(n) != 2L))) {
    cli::cli_abort(
      "Please provide exactly two estimates and two sample sizes (n) for a difference in proportions.",
      call = call
    )
  }

  if (
    !type %in% "propdiff" &&
      !is.null(est_names) &&
      (length(est_names) != length(estimate))
  ) {
    cli::cli_abort(
      "Length of estimates does not match length of estimate names.",
      call = call
    )
  }

  if (type %in% "propdiff" && !is.null(est_names) && (length(est_names) > 1L)) {
    cli::cli_abort(
      "Provide only one estimate name for a proportion difference.",
      call = call
    )
  }

  if (!is.null(est_names) && anyDuplicated(est_names) > 0L) {
    cli::cli_abort(
      c(
        "Estimate names must be unique.",
        x = "Duplicated name{?s}: {.val {unique(est_names[duplicated(est_names)])}}."
      ),
      call = call
    )
  }

  if (
    !type %in% "propdiff" &&
      isFALSE(together) &&
      (length(estimate) > 1L) &&
      !is.null(nrow) &&
      !is.null(ncol) &&
      (nrow * ncol < length(estimate))
  ) {
    cli::cli_abort(
      "nrow * ncol must be greater than or equal the number of estimates to be plotted if together = FALSE.",
      call = call
    )
  }

  #-----------------------------------------------------------------------------
  # Values by type
  #-----------------------------------------------------------------------------

  if (type %in% "ttest") {
    stderr_implied <- estimate / tstat

    if (!all(is.finite(stderr_implied)) || any(stderr_implied <= 0)) {
      cli::cli_abort(
        c(
          "The t-statistics must be non-zero and have the same sign as the estimates.",
          i = "The standard error is calculated as {.code estimate / tstat} and must be positive."
        ),
        call = call
      )
    }
  }

  if (type %in% "var" && any(estimate <= 0)) {
    cli::cli_abort("Variance estimates must be larger than 0.", call = call)
  }

  if (type %in% "var" && any(n <= 1)) {
    cli::cli_abort(
      "Sample sizes for variance estimates must be larger than 1.",
      call = call
    )
  }

  if (type %in% c("prop", "propdiff") && (any(estimate < 0) || any(estimate > 1))) {
    cli::cli_abort(
      "Please provide proportion estimates as decimals between 0 and 1.",
      call = call
    )
  }

  if (type %in% "propdiff" && !plot_type %in% c("p_val", "s_val")) {
    cli::cli_abort(
      "Currently, only P-value functions (p_val) and S-value functions (s_val) are allowed for difference in proportions.",
      call = call
    )
  }

  if (
    type %in%
      "propdiff" &&
      (((estimate[1] * n[1]) %% 1 >= 0.05) ||
        ((estimate[2] * n[2]) %% 1 >= 0.05))
  ) {
    cli::cli_warn(
      "Number of successes (i.e. estimate*n) of proportions not integer! The the number of successes was rounded (i.e. round(estimate*n))."
    )
  }

  if (type %in% c("pearson", "spearman", "kendall") && any(abs(estimate) >= 1)) {
    cli::cli_abort(
      "Correlation coefficients must lie strictly between -1 and 1.",
      call = call
    )
  }

  if (
    type %in%
      c("pearson", "spearman", "kendall") &&
      !is.null(null_values) &&
      any(abs(null_values) > 1)
  ) {
    cli::cli_abort(
      "Null values for correlations must lie between -1 and 1.",
      call = call
    )
  }

  if (
    type %in%
      "prop" &&
      !is.null(null_values) &&
      (any(null_values <= 0) || any(null_values >= 1))
  ) {
    cli::cli_abort(
      "Null values for proportions must lie between 0 and 1 (excluding).",
      call = call
    )
  }

  if (
    (type %in% "pearson" && any(n < 4)) ||
      (type %in% "spearman" && any(n < 4)) ||
      (type %in% "kendall" && any(n < 5))
  ) {
    cli::cli_abort(
      "Sample size must be at least 4 for Pearson and Spearman and at least 5 for Kendall's correlation.",
      call = call
    )
  }

  if (type %in% "spearman" && (any(estimate >= 0.9) || any(n < 10))) {
    cli::cli_warn(
      "Approximations for Spearman's correlation are only valid for r < 0.9 and n >= 10. Interpret with caution."
    )
  }

  if (type %in% "kendall" && any(estimate >= 0.8)) {
    cli::cli_warn(
      "Approximations for Kendall's correlation are only valid for r < 0.8. Interpret with caution."
    )
  }

  invisible(NULL)
}
