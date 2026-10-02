# Argument lists of the configurations stored in `regression.rds`.
#
# This file is sourced by `make-regression.R` (which creates the fixture) and
# by `helper-regression.R` (which re-runs the configurations in the tests). It
# must therefore not depend on testthat or on the package.

# Example of a user-supplied transformation (relative standard error)
rse_fun <- function(x) 100 * (1 - exp(x))
rse_fun_inv <- function(x) log(1 - (x / 100))

regression_configs <- function() {
  z <- list(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    n_values = 100L
  )

  list(
    z_basic = c(
      z,
      list(
        conf_level = c(0.95, 0.8),
        null_values = c(0, 0.3),
        plot_counternull = TRUE,
        xlim = c(-0.5, 1.5)
      )
    ),
    z_one_sided = c(
      z,
      list(
        alternative = "one_sided",
        conf_level = c(0.95, 0.8),
        null_values = 0,
        plot_counternull = TRUE
      )
    ),
    z_multi = list(
      estimate = c(0, 1, 2),
      stderr = c(1, 0.5, 2),
      type = "general_z",
      n_values = 100L,
      est_names = c("a", "b", "c"),
      conf_level = 0.95,
      null_values = 0,
      plot_counternull = TRUE
    ),
    linreg = list(
      estimate = -0.02143,
      df = 43,
      stderr = 0.02394,
      type = "linreg",
      n_values = 100L,
      conf_level = c(0.95, 0.9, 0.8),
      null_values = 0
    ),
    linreg_multi = list(
      estimate = c(1, 2),
      df = c(10, 20),
      stderr = c(0.5, 0.7),
      type = "linreg",
      n_values = 100L,
      conf_level = 0.95
    ),
    gammareg = list(
      estimate = 1,
      df = 15,
      stderr = 0.4,
      type = "gammareg",
      n_values = 100L,
      conf_level = 0.95
    ),
    ttest = list(
      estimate = 0.8,
      tstat = 2.1,
      df = 25,
      type = "ttest",
      n_values = 100L,
      conf_level = 0.95,
      null_values = 0,
      plot_counternull = TRUE
    ),
    logreg_exp = list(
      estimate = 0.804037549,
      stderr = 0.331819298,
      type = "logreg",
      n_values = 100L,
      conf_level = c(0.95, 0.9, 0.8),
      null_values = 0,
      trans = "exp",
      xlim = log(c(0.7, 5.2)),
      plot_counternull = TRUE
    ),
    coxreg_custom = list(
      estimate = log(0.72),
      stderr = 0.187618,
      type = "coxreg",
      n_values = 100L,
      est_names = "RSE",
      conf_level = c(0.95, 0.8, 0.5),
      null_values = rse_fun_inv(0),
      trans = rse_fun,
      xlim = rse_fun_inv(c(-30, 60))
    ),
    spearman = list(
      estimate = 0.5,
      n = 40,
      type = "spearman",
      n_values = 100L,
      conf_level = 0.95,
      null_values = 0
    ),
    kendall = list(
      estimate = 0.4,
      n = 40,
      type = "kendall",
      n_values = 100L,
      conf_level = 0.95
    ),
    var = list(
      estimate = 4,
      n = 20,
      type = "var",
      n_values = 100L,
      conf_level = c(0.95, 0.8),
      null_values = 3,
      plot_counternull = TRUE
    ),
    var_one_sided = list(
      estimate = c(4, 6),
      n = c(20, 30),
      type = "var",
      n_values = 100L,
      alternative = "one_sided",
      conf_level = 0.95
    ),
    prop = list(
      estimate = 0.3,
      n = 50,
      type = "prop",
      n_values = 100L,
      conf_level = 0.95,
      null_values = 0.5,
      plot_counternull = TRUE
    ),
    prop_multi = list(
      estimate = c(0.3, 0.6),
      n = c(50, 80),
      type = "prop",
      n_values = 100L,
      conf_level = 0.95
    ),
    propdiff = list(
      estimate = c(68 / 100, 98 / 150),
      n = c(100, 150),
      type = "propdiff",
      n_values = 200L,
      conf_level = c(0.95, 0.9, 0.8),
      null_values = 0
    ),
    propdiff_one_sided = list(
      estimate = c(0.4, 0.3),
      n = c(60, 60),
      type = "propdiff",
      n_values = 200L,
      alternative = "one_sided",
      conf_level = 0.95
    ),
    # Exact distribution for Pearson's correlation (slow)
    pearson = list(
      estimate = 0.3,
      n = 30,
      type = "pearson",
      n_values = 40L,
      conf_level = c(0.95, 0.8),
      null_values = 0,
      plot_counternull = TRUE
    ),
    pearson_one_sided = list(
      estimate = 0.3,
      n = 30,
      type = "pearson",
      n_values = 40L,
      alternative = "one_sided",
      conf_level = 0.95
    )
  )
}
