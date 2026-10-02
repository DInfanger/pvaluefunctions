#========================================================================
# Construct confidence distributions, densities and p-value functions
# Author: Denis Infanger
# Creation date (dd.mm.yyyy): 22.09.2018
#========================================================================

#' Create and Plot \emph{P}-Value Functions, S-Value Functions, Confidence Distributions and Confidence Densities
#'
#' The function \code{conf_dist} generates confidence distributions (cdf), confidence densities (pdf), Shannon suprisal (s-value) functions and \emph{p}-value functions for several commonly used estimates. In addition, counternulls (see Rosenthal et al. 1994), point estimates and the area under the confidence curve (AUCC) are calculated.
#'
#' \emph{P}-value functions and confidence intervals are calculated based on the \emph{t}-distribution for \emph{t}-tests, linear regression coefficients, and gamma regression models (GLM). The normal distribution is used for logistic regression, poisson regression and cox regression models. For correlation coefficients, Fisher's transform is used using the corresponding variances (see Bonett et al. 2000). \emph{P}-value functions and confidence intervals for variances are constructed using the Chi2 distribution. Finally, Wilson's score intervals are used for one proportion. For differences of proportions, the Wilson score interval with continuity correction is used (Newcombe 1998).
#'
#' @param estimate Numerical vector containing the estimate(s).
#' @param n Numerical vector containing the sample size(s). Required for correlations, variances, proportions and differences between proportions. Must be equal the number of estimates.
#' @param df Numerical vector containing the degrees of freedom. Required for statistics based on the \emph{t}-distribution (e.g. linear regression) and \emph{t}-tests. Must be equal the number of estimates.
#' @param stderr Numerical vector containing the standard error(s) of the estimate(s). Required for statistics based on the \emph{t}-distribution (e.g. linear regression) and the normal distribution (e.g. logistic regression). Must be equal the number of estimate(s).
#' @param tstat Numerical vector containing the \emph{t}-statistic(s). Required for \emph{t}-tests (means and mean differences). Must be equal the number of estimates.
#' @param type String indicating the type of the estimate. Must be one of the following: \code{ttest}, \code{linreg}, \code{gammareg}, \code{general_t}, \code{logreg}, \code{poisreg}, \code{coxreg}, \code{general_z}, \code{pearson}, \code{spearman}, \code{kendall}, \code{var}, \code{prop}, \code{propdiff}.
#' @param plot_type String indicating the type of plot. Must be one of the following: \code{cdf} (confidence distribution), \code{pdf} (confidence density), \code{p_val} (\emph{p}-value function, the default), \code{s_val} (Surprisal value functions). For differences between independent proportions, only \emph{p}-value functions and Surprisal values are available.
#' @param n_values (optional) Integer indicating the number of points that are used to generate the graphics. Must be at least 2. The higher this number, the higher the computation time and resolution.
#' @param est_names (optional) String vector indicating the names of the estimate(s). Must be equal the number of estimates.
#' @param conf_level (optional) Numerical vector indicating the confidence level(s). Must be between 0 and 1.
#' @param null_values (optional) Numerical vector indicating the null value(s) in the plot on the \emph{untransformed (original)} scale. For example: The null values for an odds ratio of 1 is 0 on the log-odds scale. If x limits are specified with \code{xlim}, all null values outside of the specified x limits are ignored for plotting and a message is printed.
#' @param trans (optional) String indicating the transformation function that will be applied to the estimates and confidence curves. For example: \code{"exp"} for an exponential transformation of the log-odds in logistic regression. Can be a custom function, given either as a function or as the name of a function (as a string). A custom function must be vectorized and receives the values to be transformed as its first argument.
#' @param alternative String indicating if the confidence level(s) are two-sided or one-sided. Must be one of the following: \code{two_sided}, \code{one_sided}.
#' @param log_yaxis Logical. Indicating if a portion of the y-axis should be displayed on the logarithmic scale.
#' @param cut_logyaxis Numerical value indicating the threshold below which the y-axis will be displayed logarithmically. Must lie between 0 and 1.
#' @param xlim (optional) Optional numerical vector of length 2 (x1, x2) indicating the limits of the x-axis on the \emph{untransformed} scale if \code{trans} is not \code{identity}. The scale of the x-axis set by \code{x_scale} does not affect the x limits. For example: If you want to plot \emph{p}-value functions for odds ratios from logistic regressions, the limits have to be given on the log-odds scale if \code{trans = "exp"}. Note that x1 > x2 is allowed but then x2 will be the left limit and x1 the right limit (i.e. the limits are sorted before plotting). Null values (specified in \code{null_values}) that are outside of the specified limits are ignored and a message is printed.
#' @param together Logical. Indicating if graphics for multiple estimates should be displayed together or on separate plots.
#' @param plot_legend Logical. Indicating if a legend should be plotted if multiple curves are plotted together with different colors (i.e. \code{together = TRUE} and \code{same_color = FALSE}).
#' @param same_color Logical. Indicating if curves should be distinguished using colors if they are plotted together (i.e. \code{together = TRUE}). Setting this to FALSE also disables the default behavior that the two halves of the curves are plotted in different colors for a one-sided alternative.
#' @param col String indicating the colour of the curves. Only relevant for single curves, multiple curves not plotted together (i.e. \code{together = FALSE}) and multiple curves plotted together but with the option \code{same_color} set to \code{TRUE}.
#' @param nrow (optional) Integer greater than 0 indicating the number of rows when \code{together = FALSE} is specified for multiple estimates. Used in \code{facet_wrap} in ggplot2.
#' @param ncol (optional) Integer greater than 0 indicating the number of columns when \code{together = FALSE} is specified for multiple estimates. Used in \code{facet_wrap} in ggplot2.
#' @param plot_p_limit Numerical value indicating the lower limit of the y-axis. Must be greater than 0 for a logarithmic scale (i.e. \code{log_yaxis = TRUE}). The default is to omit plotting \emph{p}-values smaller than 1 - 0.999 = 0.001.
#' @param plot_counternull Logical. Indicating if the counternull should be plotted as a point. Only available for \emph{p}-value functions and s-value functions. Counternull values that are outside of the plotted functions are not shown.
#' @param title (optional) String containing a title of the plot.
#' @param xlab (optional) String indicating the label of the x-axis.
#' @param ylab (optional) String indicating the title for the primary (left) y-axis.
#' @param ylab_sec (optional) String indicating the title for the secondary (right) y-axis.
#' @param inverted Logical. Indicating the orientation of the y-axis for the \emph{P}-value function (\code{p_val}), S-value function (\code{s_val}) and the confidence distribution (\code{cdf}). By default (i.e. \code{inverted = FALSE}) small \emph{P}-values are plotted at the bottom and large ones at the top so that the cusp of the \emph{P}-value function is at the top. By setting \code{inverted = TRUE}, the y-axis is inverted. Ignored for confidence densities.
#' @param x_scale String indicating the scaling of the x-axis. The default is to scale the x-axis logarithmically if the transformation specified in \code{trans} is "exp" (exponential) and linearly otherwise. The option \code{linear} (can be abbreviated) forces a linear scaling and the option \code{logarithm} (can be abbreviated) forces a logarithmic scaling, regardless what has been specified in \code{trans}.
#' @param plot Logical. Should a plot be created (\code{TRUE}, the default) or not (\code{FALSE}). \code{FALSE} can be useful if users want to create their own plots using the returned data from the function. If \code{FALSE}, no ggplot2 object is returned.

#' @return \code{conf_dist} returns four data frames and if \code{plot = TRUE} was specified, a ggplot2-plot object: \code{res_frame} (contains parameter values (e.g. mean differences, odds ratios etc.), \emph{p}-values (one- and two-sided), s-values, confidence distributions and densities, variable names and type of hypothesis), \code{conf_frame} (contains the specified confidence level(s) and the corresponding lower and upper limits as well as the corresponding variable name), \code{counternull_frame} (contains the counternull and the corresponding null values), \code{point_est} (contains the mean, median and mode point estimates) and if \code{plot = TRUE} was specified, \code{aucc_frame} contains the estimated AUCC (area under the confidence curve, see Berrar 2017) calculated by trapezoidal integration on the untransformed scale. Also provides the proportion of the aucc that lies above the null value(s) if they are provided. \code{plot} (a ggplot2 object).
#' @references Bender R, Berg G, Zeeb H. Tutorial: using confidence curves in medical research. \emph{Biom J.} 2005;47(2):237-247.
#'
#' Berrar D. Confidence curves: an alternative to null hypothesis significance testing for the comparison of classifiers. \emph{Mach Learn.} 2017;106:911-949.
#'
#' Bonett DG, Wright TA. Sample size requirements for estimating Pearson, Kendall and Spearman correlations. \emph{Psychometrika.} 2000;65(1):23-28.
#'
#' Cole SR, Edwards JK, Greenland S. Surprise! \emph{Am J Epidemiol.} 2021:190(2):191-193.
#'
#' Infanger D, Schmidt-Trucksäss A. \emph{P} value functions: An underused method to present research results and to promote quantitative reasoning. \emph{Stat Med.} 2019;38:4189-4197.
#'
#' Newcombe RG. Interval estimation for the difference between independent proportions: comparison of eleven methods. \emph{Stat Med.} 1998;17:873-890.
#'
#' Poole C. Confidence intervals exclude nothing. \emph{Am J Public Health.} 1987;77(4):492-493.
#'
#' Poole C. Beyond the confidence interval. \emph{Am J Public Health.} 1987;77(2):195-199.
#'
#' Rafi Z, Greenland S. Semantic and cognitive tools to aid statistical science: replace confidence and significance by compatibility and surprise. \emph{BMC Med Res Methodol} 2020;20:244.
#'
#' Rosenthal R, Rubin D. The counternull value of an effect size: a new statistic. \emph{Psychological Science.} 1994;5(6):329-334.
#'
#' Rothman KJ, Greenland S, Lash TL. Modern epidemiology. 3rd ed. Philadelphia, PA: Wolters Kluwer; 2008.
#'
#' Schweder T, Hjort NL. Confidence, likelihood, probability: statistical inference with confidence distributions. New York, NY: Cambridge University Press; 2016.
#'
#' Sullivan KM, Foster DA. Use of the confidence interval function. \emph{Epidemiology.} 1990;1(1):39-42.
#'
#' Xie Mg, Singh K. Confidence distribution, the frequentist distribution estimator of a parameter: A review. \emph{Internat Statist Rev.} 2013;81(1):3-39.
#'
#' @examples
#'
#' #======================================================================================
#' # Create a p-value function for an estimate using the normal distribution
#' #======================================================================================
#'
#' res <- conf_dist(
#'   estimate = c(-0.13)
#'   , stderr = c(0.224494)
#'   , type = "general_z"
#'   , plot_type = "p_val"
#'   , n_values = 1e4L
#'   , est_names = c("Parameter value")
#'   , log_yaxis = FALSE
#'   , cut_logyaxis = 0.05
#'   , conf_level = c(0.95)
#'   , null_values = c(0)
#'   , trans = "identity"
#'   , alternative = "two_sided"
#'   , xlab = "Var"
#'   , xlim = c(-1, 1)
#'   , together = TRUE
#'   , plot_p_limit = 1 - 0.9999
#'   , plot_counternull = TRUE
#'   , title = NULL
#'   , ylab = NULL
#'   , ylab_sec = NULL
#'   , inverted = FALSE
#'   , x_scale = "default"
#'   , plot = TRUE
#' )
#'
#' #======================================================================================
#' # P-value function for a single regression coefficient (Agriculture in the model below)
#' #======================================================================================
#'
#' mod <- lm(Infant.Mortality~Agriculture + Fertility + Examination, data = swiss)
#' summary(mod)
#'
#' res <- conf_dist(
#'   estimate = c(-0.02143)
#'   , df = c(43)
#'   , stderr = (0.02394)
#'   , type = "linreg"
#'   , plot_type = "p_val"
#'   , n_values = 1e4L
#'   , conf_level = c(0.95, 0.90, 0.80)
#'   , null_values = c(0)
#'   , trans = "identity"
#'   , alternative = "two_sided"
#'   , log_yaxis = TRUE
#'   , cut_logyaxis = 0.05
#'   , xlab = "Coefficient Agriculture"
#'   , together = FALSE
#'   , plot_p_limit = 1 - 0.999
#'   , plot_counternull = FALSE
#'   , title = NULL
#'   , ylab = NULL
#'   , ylab_sec = NULL
#'   , inverted = FALSE
#'   , x_scale = "default"
#'   , plot = TRUE
#' )
#'
#' #=======================================================================================
#' # P-value function for an odds ratio (logistic regression), plotted with inverted y-axis
#' #=======================================================================================
#'
#' res <- conf_dist(
#'   estimate = c(0.804037549)
#'   , stderr = c(0.331819298)
#'   , type = "logreg"
#'   , plot_type = "p_val"
#'   , n_values = 1e4L
#'   , est_names = c("GPA")
#'   , conf_level = c(0.95, 0.90, 0.80)
#'   , null_values = c(log(1)) # null value on the log-odds scale
#'   , trans = "exp"
#'   , alternative = "two_sided"
#'   , log_yaxis = FALSE
#'   , cut_logyaxis = 0.05
#'   , xlab = "Odds Ratio (GPA)"
#'   , xlim = log(c(0.7, 5.2)) # axis limits on the log-odds scale
#'   , together = FALSE
#'   , plot_p_limit = 1 - 0.999
#'   , plot_counternull = TRUE
#'   , title = NULL
#'   , ylab = NULL
#'   , ylab_sec = NULL
#'   , inverted = TRUE
#'   , x_scale = "default"
#'   , plot = TRUE
#' )
#'
#' #======================================================================================
#' # Difference between two independent proportions: Newcombe with continuity correction
#' #======================================================================================
#'
#' res <- conf_dist(
#'   estimate = c(68/100, 98/150)
#'   , n = c(100, 150)
#'   , type = "propdiff"
#'   , plot_type = "p_val"
#'   , n_values = 1e4L
#'   , conf_level = c(0.95, 0.90, 0.80)
#'   , null_values = c(0)
#'   , trans = "identity"
#'   , alternative = "two_sided"
#'   , log_yaxis = FALSE
#'   , cut_logyaxis = 0.05
#'   , xlab = "Difference between proportions"
#'   , together = FALSE
#'   , col = "#A52A2A" # Color curve in auburn
#'   , plot_p_limit = 1 - 0.9999
#'   , plot_counternull = FALSE
#'   , title = NULL
#'   , ylab = NULL
#'   , ylab_sec = NULL
#'   , inverted = FALSE
#'   , x_scale = "default"
#'   , plot = TRUE
#' )
#'
#' #======================================================================================
#' # Difference between two independent proportions: Agresti & Caffo
#' #======================================================================================
#'
#' # First proportion
#' x1 <- 8
#' n1 <- 40
#'
#' # Second proportion
#' x2 <- 11
#' n2 <- 30
#'
#' # Apply the correction
#' p1hat <- (x1 + 1)/(n1 + 2)
#' p2hat <- (x2 + 1)/(n2 + 2)
#'
#' # The original estimator
#' est0 <- (x1/n1) - (x2/n2)
#'
#' # The unmodified estimator and its standard error using the correction
#'
#' est <- p1hat - p2hat
#' se <- sqrt(((p1hat*(1 - p1hat))/(n1 + 2)) + ((p2hat*(1 - p2hat))/(n2 + 2)))
#'
#' res <- conf_dist(
#'   estimate = c(est)
#'   , stderr = c(se)
#'   , type = "general_z"
#'   , plot_type = "p_val"
#'   , n_values = 1e4L
#'   , log_yaxis = FALSE
#'   , cut_logyaxis = 0.05
#'   , conf_level = c(0.95, 0.99)
#'   , null_values = c(0, 0.3)
#'   , trans = "identity"
#'   , alternative = "two_sided"
#'   , xlab = "Difference of proportions"
#'   , together = FALSE
#'   , plot_p_limit = 1 - 0.9999
#'   , plot_counternull = FALSE
#'   , title = "P-value function for the difference of two independent proportions"
#'   , ylab = NULL
#'   , ylab_sec = NULL
#'   , inverted = FALSE
#'   , x_scale = "default"
#'   , plot = TRUE
#' )
#'
#' #========================================================================================
#' # P-value function and confidence distribution for the relative survival effect (1 - HR%)
#' # Replicating Figure 1 in Bender et al. (2005)
#' #========================================================================================
#'
#' # Define the transformation function and its inverse for the relative survival effect
#'
#' rse_fun <- function(x){ # x is the log-hazard ratio
#'   100 * (1 - exp(x))
#' }
#'
#' rse_fun_inv <- function(x){
#'   log(1 - (x / 100))
#' }
#'
#' res <- conf_dist(
#'   estimate = log(0.72)
#'   , stderr = 0.187618
#'   , type = "coxreg"
#'   , plot_type = "p_val"
#'   , n_values = 1e4L
#'   , est_names = c("RSE")
#'   , conf_level = c(0.95, 0.8, 0.5)
#'   , null_values = rse_fun_inv(0)
#'   , trans = "rse_fun"
#'   , alternative = "two_sided"
#'   , log_yaxis = FALSE
#'   , cut_logyaxis = 0.05
#'   , xlab = "Relative survival effect (1 - HR%)"
#'   , xlim = rse_fun_inv(c(-30, 60))
#'   , together = FALSE
#'   , plot_p_limit = 1 - 0.999
#'   , plot_counternull = TRUE
#'   , inverted = TRUE
#'   , title = "Figure 1 in Bender et al. (2005)"
#'   , x_scale = "default"
#'   , plot = TRUE
#' )
#'
#' @export

conf_dist <- function(
  estimate = NULL,
  n = NULL,
  df = NULL,
  stderr = NULL,
  tstat = NULL,
  type = NULL,
  plot_type = c("p_val", "s_val", "cdf", "pdf"),
  n_values = 1e4L,
  est_names = NULL,
  conf_level = NULL,
  null_values = NULL,
  trans = "identity",
  alternative = c("two_sided", "one_sided"),
  log_yaxis = FALSE,
  cut_logyaxis = 0.05,
  xlab = NULL,
  xlim = NULL,
  together = FALSE,
  plot_legend = TRUE,
  same_color = FALSE,
  col = "black",
  nrow = NULL,
  ncol = NULL,
  plot_p_limit = (1 - 0.999),
  plot_counternull = FALSE,
  title = NULL,
  ylab = NULL,
  ylab_sec = NULL,
  inverted = FALSE,
  x_scale = c("default", "linear", "logarithm"),
  plot = TRUE
) {
  #-----------------------------------------------------------------------------
  # Safety checks and clean ups
  #-----------------------------------------------------------------------------

  trans_caller_env <- parent.frame()
  alternative <- match.arg(alternative)
  plot_type <- match.arg(plot_type)
  x_scale <- match.arg(x_scale)

  validate_args(
    estimate = estimate,
    n = n,
    df = df,
    stderr = stderr,
    tstat = tstat,
    type = type,
    plot_type = plot_type,
    n_values = n_values,
    est_names = est_names,
    conf_level = conf_level,
    null_values = null_values,
    alternative = alternative,
    log_yaxis = log_yaxis,
    cut_logyaxis = cut_logyaxis,
    xlab = xlab,
    xlim = xlim,
    together = together,
    plot_legend = plot_legend,
    same_color = same_color,
    nrow = nrow,
    ncol = ncol,
    plot_p_limit = plot_p_limit,
    plot_counternull = plot_counternull,
    inverted = inverted,
    plot = plot
  )

  type <- tolower(type)
  plot_p_limit <- round(plot_p_limit, 10)
  cut_logyaxis <- round(cut_logyaxis, 10)

  # match.arg() returns the full option name but the code below uses "log"
  if (x_scale == "logarithm") {
    x_scale <- "log"
  }

  # Confidence levels that are too low for the alternative cannot be displayed
  if (!is.null(conf_level)) {
    conf_level_too_low <- switch(
      alternative,
      two_sided = (1 - conf_level) >= 1,
      one_sided = (2 - 2 * conf_level) >= 1
    )

    if (any(conf_level_too_low)) {
      cli::cli_inform(
        "Ignoring confidence levels that are too low for the {.val {alternative}} alternative: {.val {conf_level[conf_level_too_low]}}."
      )
      conf_level <- conf_level[!conf_level_too_low]
    }
  }

  if (
    type %in%
      c("pearson", "spearman", "kendall", "var", "prop", "propdiff") &&
      !is_identity_trans(trans)
  ) {
    trans <- "identity"
    cli::cli_inform("Transformation changed to identity.")
  }

  # `trans` is now the (lower-case) name, `trans_fun` the function to apply
  trans_spec <- resolve_trans(trans, env = trans_caller_env)
  trans <- trans_spec$name
  trans_fun <- trans_spec$fun

  if (!is.null(xlim)) {
    xlim <- sort(xlim, decreasing = FALSE)
  }

  if (is.null(est_names)) {
    est_names <- if (type %in% "propdiff") 1L else seq_along(estimate)
  }

  if ((trans %in% "exp") && (x_scale %in% "default")) {
    x_scale <- "log"
  }

  if (length(col) > 1) {
    cli::cli_warn("Only first color is used: {col[1]}")
    col <- col[1]
  }

  #-----------------------------------------------------------------------------
  # Calculate the confidence distributions/densities and p-value curves
  #-----------------------------------------------------------------------------

  if (type %in% "ttest") {
    stderr <- estimate / tstat

    res <- cdist_t(
      estimate = estimate,
      stderr = stderr,
      df = df,
      n_values = n_values,
      conf_level = conf_level,
      alternative = alternative,
      null_values = null_values
    )
  } else if (type %in% c("linreg", "gammareg", "general_t")) {
    res <- cdist_t(
      estimate = estimate,
      stderr = stderr,
      df = df,
      n_values = n_values,
      conf_level = conf_level,
      alternative = alternative,
      null_values = null_values
    )
  } else if (type %in% c("logreg", "poisreg", "coxreg", "general_z")) {
    res <- cdist_z(
      estimate = estimate,
      stderr = stderr,
      n_values = n_values,
      conf_level = conf_level,
      null_values = null_values,
      alternative = alternative
    )
  } else if (type %in% c("pearson", "spearman", "kendall")) {
    # Calculate approximate standard error for each type of correlation coefficient
    # Ref 1: Bonett & Wright (2000): Sample size requirements for estimating Pearson, Kendall and Spearman correlations
    # Ref 2: Fieller, Hartley, Pearson (1957): Tests for rank correlation coefficients I.

    if (type %in% c("spearman", "kendall")) {
      # Use approximations for Spearman and Kendall

      stderr <- switch(
        type,
        pearson = 1 / sqrt(n - 3),
        spearman = sqrt((1 + (estimate)^2 / 2) / (n - 3)),
        kendall = sqrt(0.437 / (n - 4))
      )

      res <- cdist_corr(
        estimate = estimate,
        stderr = stderr,
        n = n,
        n_values = n_values,
        conf_level = conf_level,
        null_values = null_values,
        alternative = alternative
      )
    } else if (type %in% "pearson") {
      # Use exact distribution for Pearson

      res <- cdist_corr_exact(
        estimate = estimate,
        n = n,
        n_values = n_values,
        conf_level = conf_level,
        null_values = null_values,
        alternative = alternative
      )
    }
  } else if (type %in% "var") {
    res <- cdist_var(
      estimate = estimate,
      n = n,
      n_values = n_values,
      conf_level = conf_level,
      null_values = null_values,
      alternative = alternative
    )
  } else if (type %in% "prop") {
    res <- cdist_prop1(
      estimate = estimate,
      n = n,
      n_values = n_values,
      conf_level = conf_level,
      null_values = null_values,
      alternative = alternative
    )
  } else if (type %in% "propdiff") {
    res <- cdist_propdiff(
      estimate = estimate,
      n = n,
      n_values = n_values,
      conf_level = conf_level,
      null_values = null_values,
      alternative = alternative
    )

    estimate <- estimate[1] - estimate[2]
  }

  #-----------------------------------------------------------------------------
  # Calculate Shannon-surprisal value (S-value)
  #-----------------------------------------------------------------------------

  res$res_frame$s_val <- -log2(res$res_frame$p_two)

  #-----------------------------------------------------------------------------
  # Calculate AUCC (area under the confidence curve), see Berrar (2017) Mach Learn 106:911-494
  # Also calculate the area above the null values
  #-----------------------------------------------------------------------------

  res$aucc_frame <- compute_aucc(res$res_frame, length(estimate), null_values)

  #-----------------------------------------------------------------------------
  # Add an indicator variable for the type of hypothesis for plotting
  #-----------------------------------------------------------------------------

  res$res_frame$hypothesis <- NA

  if (alternative %in% "one_sided") {
    if (type %in% "var") {
      for (i in seq_along(estimate)) {
        estimate[i] <- res$point_est$est_median[i]
      }
    }

    for (i in seq_along(estimate)) {
      res$res_frame$hypothesis[
        res$res_frame$variable == i & res$res_frame$values < estimate[i]
      ] <- 1 # greater
      res$res_frame$hypothesis[
        res$res_frame$variable == i & res$res_frame$values >= estimate[i]
      ] <- (-1) # less
    }
  }

  res$res_frame$hypothesis <- factor(
    res$res_frame$hypothesis,
    levels = c(-1, 1),
    labels = c("less", "greater")
  )

  #-----------------------------------------------------------------------------
  # Calculate the limits of the x-axis if not provided
  #-----------------------------------------------------------------------------

  # If there are limits given, take those in any case
  if (is.null(xlim)) {
    res_tmp <- res$res_frame
    res_tmp$values[res$res_frame$p_two < plot_p_limit] <- NA

    plot_range_tmp <- range(res_tmp$values, na.rm = TRUE)

    if (is.null(null_values)) {
      # If no limits and no null values given, take the plot_p_limit
      xlim <- plot_range_tmp
    } else {
      # If no limits but null values given, look if the plot_p_limits are outside the null_values
      xlim <- c(NA, NA)

      # If the smallest null value is outside of the plotting area, set the lower limit to that null value
      if (min(null_values, na.rm = TRUE) <= plot_range_tmp[1]) {
        xlim[1] <- min(null_values, na.rm = TRUE)
      } else {
        xlim[1] <- plot_range_tmp[1]
      }

      # If the largest null value is outside of the plotting area, set the upper limit to that null value
      if (max(null_values, na.rm = TRUE) >= plot_range_tmp[2]) {
        xlim[2] <- max(null_values, na.rm = TRUE)
      } else {
        xlim[2] <- plot_range_tmp[2]
      }
    }

    rm(res_tmp, plot_range_tmp)
  }

  #-----------------------------------------------------------------------------
  # Text frame coordinates and contents for plotting the confidence levels
  #-----------------------------------------------------------------------------

  if (!is.null(conf_level)) {
    text_frame <- make_text_frame(conf_level, alternative, trans)
  }

  #-----------------------------------------------------------------------------
  # Apply transformations if applicable
  #-----------------------------------------------------------------------------

  if (!trans %in% "identity") {
    res$res_frame$values <- trans_fun(res$res_frame$values)

    point_est_cols <- c("est_mean", "est_median", "est_mode")
    res$point_est[point_est_cols] <- lapply(
      res$point_est[point_est_cols],
      trans_fun
    )

    if (!is.null(conf_level)) {
      if (!trans %in% "exp") {
        text_frame$theor_values[is.finite(
          text_frame$theor_values
        )] <- trans_fun(
          text_frame$theor_values[is.finite(text_frame$theor_values)]
        )
      }
      res$conf_frame$lwr <- trans_fun(res$conf_frame$lwr)
      res$conf_frame$upr <- trans_fun(res$conf_frame$upr)
    }

    if (!is.null(null_values)) {
      res$counternull_frame$counternull <- trans_fun(
        res$counternull_frame$counternull
      )
      res$counternull_frame$null_value <- trans_fun(
        res$counternull_frame$null_value
      )
    }

    xlim <- sort(trans_fun(xlim), decreasing = FALSE)
  }

  #-----------------------------------------------------------------------------
  # Cutoff for nicer plotting
  #-----------------------------------------------------------------------------

  p_cutoff <- ifelse(
    alternative %in% c("two_sided"),
    plot_p_limit,
    plot_p_limit * 2
  )

  if (plot_type %in% c("p_val")) {
    res$res_frame$values[res$res_frame$p_two < p_cutoff] <- NA
    res$res_frame$p_two[res$res_frame$p_two < p_cutoff] <- NA
    res$res_frame$p_one[res$res_frame$p_one < p_cutoff] <- NA
  }

  if (plot_type %in% c("s_val")) {
    outside_ind <- which(res$res_frame$p_two < p_cutoff)

    res$res_frame$values[outside_ind] <- NA
    res$res_frame$p_two[outside_ind] <- NA
    res$res_frame$p_one[outside_ind] <- NA
    res$res_frame$s_val[outside_ind] <- NA
  }

  if (!is.null(conf_level) && any(text_frame$p_value < p_cutoff)) {
    text_frame <- text_frame[-which(text_frame$p_value < p_cutoff), ]
  }

  #-----------------------------------------------------------------------------
  # Add counternull values to the result frame for plotting
  #-----------------------------------------------------------------------------

  if (
    !is.null(null_values) &&
      isTRUE(plot_counternull) &&
      plot_type %in% c("p_val", "s_val")
  ) {
    res$res_frame$counternull <- NA

    null_values_trans <- trans_fun(null_values)

    for (i in seq_along(estimate)) {
      plot_range_tmp <- range(
        res$res_frame$values[res$res_frame$variable == i],
        na.rm = TRUE
      )

      for (j in seq_along(null_values)) {
        counternull_tmp <- res$counternull_frame$counternull[
          res$counternull_frame$variable == i &
            res$counternull_frame$null_value %in% null_values_trans[j]
        ]

        # Counternulls can be missing (e.g. null values in the "gap" of
        # differences of proportions), in which case there is nothing to plot
        counternull_tmp <- counternull_tmp[1]

        if (
          !is.na(counternull_tmp) &&
            (counternull_tmp <= plot_range_tmp[2]) &&
            (counternull_tmp >= plot_range_tmp[1])
        ) {
          counternull_index_tmp <- which.min(abs(
            res$res_frame$values[res$res_frame$variable == i] - counternull_tmp
          ))

          res$res_frame$counternull[res$res_frame$variable == i][
            counternull_index_tmp
          ] <- res$res_frame$p_two[
            res$res_frame$variable == i
          ][which.min(abs(
            res$res_frame$values[res$res_frame$variable == i] -
              null_values_trans[j]
          ))]
        }
      }
    }

    if (plot_type %in% "s_val") {
      res$res_frame$counternull <- -log2(res$res_frame$counternull)
    }
  }

  #-----------------------------------------------------------------------------
  # Assign estimate names
  #-----------------------------------------------------------------------------

  res$point_est$variable <- factor(res$point_est$variable, labels = est_names)
  res$res_frame$variable <- factor(res$res_frame$variable, labels = est_names)
  res$aucc_frame$variable <- est_names[res$aucc_frame$variable]

  if (!is.null(conf_level)) {
    res$conf_frame$variable <- factor(
      res$conf_frame$variable,
      labels = est_names
    )
  }
  if (!is.null(null_values)) {
    res$counternull_frame$variable <- factor(
      res$counternull_frame$variable,
      labels = est_names
    )
  }

  #-----------------------------------------------------------------------------
  # Plot using ggplot2
  #-----------------------------------------------------------------------------

  # Create custom y-axis scale (mixed linear and logarithmic)

  # Transform cutoff for log-y-axis if applicable

  if (alternative %in% "one_sided") {
    if (cut_logyaxis > 0.5) {
      cut_logyaxis <- 0.5
    }
    cut_logyaxis_one <- cut_logyaxis
    cut_logyaxis <- cut_logyaxis_one * 2
  } else {
    cut_logyaxis_one <- cut_logyaxis / 2
  }

  # Labeller functions for custom log-scale

  lab_onesided <- Vectorize(function(x) {
    if (!is.na(x) && (x < cut_logyaxis_one) && (round((x %% 1) * 10) == 0)) {
      sprintf("%.5g", x)
    } else {
      sprintf("%.2f", x)
    }
  })

  lab_twosided <- Vectorize(function(x) {
    if (!is.na(x) && (x <= cut_logyaxis) && (round((x %% 1) * 10) == 0)) {
      sprintf("%.5g", x)
    } else {
      sprintf("%.1f", x)
    }
  })

  # Start plotting

  # Set variable to be plotted on y-axis depending on the type of the plot

  y_var <- switch(
    plot_type,
    p_val = "p_two",
    cdf = "conf_dist",
    pdf = "conf_dens",
    s_val = "s_val"
  )

  # Set label of the x-axis if not provided

  if (is.null(xlab)) {
    xlab <- "Estimate"
  }

  # Set label of the primary and secondary y-axis depending on the type of the plot

  if (is.null(ylab)) {
    ylab <- switch(
      plot_type,
      p_val = expression(paste(
        italic("P"),
        "-value (two-sided) / Significance level" ~ alpha,
        sep = ""
      )),
      cdf = "Confidence distribution",
      pdf = "Confidence density",
      s_val = expression(paste(
        "Surprisal in bits (two-sided ",
        ~ italic("P"),
        "-value)",
        sep = ""
      ))
    )
  }

  if (is.null(ylab_sec)) {
    ylab_sec <- switch(
      plot_type,
      p_val = expression(paste(
        italic("P"),
        "-value (one-sided) / Significance level" ~ alpha,
        sep = ""
      )),
      s_val = expression(paste(
        "Surprisal in bits (one-sided ",
        ~ italic("P"),
        "-value)",
        sep = ""
      ))
    )
  }

  # Create a ggplot2-object with no geoms

  p <- ggplot(
    res$res_frame,
    aes(x = .data$values, y = .data[[y_var]], group = .data$variable)
  ) +
    theme_bw()

  # Multiple estimates that are plotted in the same graph
  multi_together <- isTRUE(together) && length(estimate) >= 2

  # If 2 or more estimates are plotted together, differentiate them by color (if user did not specify "same_color = TRUE")

  if (multi_together && !same_color) {
    p <- p + aes(colour = .data$variable)
  }

  if (alternative == "one_sided" && !multi_together) {
    # For only one one-sided p-value curve, set the colors to black and blue
    # (or to the specified color if "same_color = TRUE")
    curve_colours <- if (same_color) c(col, col) else c("black", "#08A9CF")

    p <- p +
      geom_line(aes(colour = .data$hypothesis), linewidth = 1.5) +
      scale_colour_manual(values = curve_colours) +
      theme(
        legend.position = "none"
      )
  } else if (
    alternative == "two_sided" && (!multi_together || isTRUE(same_color))
  ) {
    # For only one two-sided p-value curve, set the color to black
    p <- p + geom_line(linewidth = 1.5, colour = col)
  } else if (isTRUE(plot_legend)) {
    # For 2 or more estimates plotted together: set the colors according to "Set1" palette
    p <- p +
      geom_line(linewidth = 1.5) +
      scale_colour_brewer(palette = "Set1", name = "") +
      theme(
        legend.position = "top",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
      )
  } else {
    # Same, but without legend
    p <- p +
      geom_line(linewidth = 1.5) +
      scale_colour_brewer(palette = "Set1", name = "") +
      theme(
        legend.position = "none"
      )
  }

  # Add the labels for the axes

  p <- p + xlab(xlab) + ylab(ylab)

  #-----------------------------------------------------------------------------
  # y-axis
  #-----------------------------------------------------------------------------

  # For p-value curves: Set the left and right y-axes, possibly with a logarithmic part

  if (plot_type %in% "p_val") {
    if (isTRUE(log_yaxis) && (p_cutoff < cut_logyaxis)) {
      lower_ylim_two <- round(10^(ceiling(round(log10(p_cutoff), 5))), 10)
      lower_ylim_one <- round(10^(ceiling(round(log10(p_cutoff / 2), 5))), 10)

      if (
        (alternative == "two_sided" && (lower_ylim_two <= cut_logyaxis)) ||
          (alternative == "one_sided" && (lower_ylim_one <= cut_logyaxis * 2))
      ) {
        # Split the breaks into two parts: i) below the cutoff i.e. the logarithmic part and ii) the linear part above the cutoff

        breaks_two <- c(
          10^(seq(log10(lower_ylim_two), log10(cut_logyaxis), by = 1)),
          seq(ceiling(cut_logyaxis / 0.1) * 0.1, 1, by = 0.1)
        )
        breaks_one <- c(
          10^(seq(log10(lower_ylim_one), log10(cut_logyaxis_one), by = 1)),
          seq(ceiling(cut_logyaxis_one / 0.05) * 0.05, 0.5, 0.05)
        )
      } else {
        breaks_two <- c(
          10^(log10(seq(ceiling(cut_logyaxis / 0.1) * 0.1, 1, by = 0.1)))
        )
        breaks_one <- c(
          10^(log10(seq(ceiling(cut_logyaxis_one / 0.05) * 0.05, 0.5, 0.05)))
        )
      }

      # Remove possible duplicates

      breaks_two <- unique(breaks_two)
      breaks_one <- unique(breaks_one)

      # The inverted axis uses the reversed limits and the reversed transformation
      if (isTRUE(inverted)) {
        y_limits <- c(1, p_cutoff)
        y_transform <- magnify_trans_log_rev(
          interval_low = cut_logyaxis,
          interval_high = 1,
          reducer = cut_logyaxis,
          reducer2 = 8
        )
      } else {
        y_limits <- c(p_cutoff, 1)
        y_transform <- magnify_trans_log(
          interval_low = cut_logyaxis,
          interval_high = 1,
          reducer = cut_logyaxis,
          reducer2 = 8
        )
      }

      p <- p +
        scale_y_continuous(
          limits = y_limits,
          breaks = breaks_two,
          labels = lab_twosided,
          transform = y_transform,
          sec.axis = sec_axis(
            transform = ~ . * (1 / 2),
            name = ylab_sec,
            breaks = breaks_one,
            labels = lab_onesided
          )
        ) +
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = p_cutoff,
          ymax = cut_logyaxis,
          alpha = 0.1,
          colour = grey(0.9)
        )
    } else {
      y_scale <- if (isTRUE(inverted)) scale_y_reverse else scale_y_continuous

      p <- p +
        y_scale(
          breaks = seq(0, 1, 0.1),
          sec.axis = sec_axis(
            ~ . * (1 / 2),
            name = ylab_sec,
            breaks = seq(0, 1, 0.1) / 2
          )
        )
    }
  }

  # For s-value curves: inverted the y-axis and transform it according to log2

  if (plot_type %in% "s_val") {
    s_val_max <- max(res$res_frame$s_val, na.rm = TRUE)

    # The inverted axis is a regular axis with the limits in ascending order
    # while the default axis is a reversed axis
    if (isTRUE(inverted)) {
      y_scale <- scale_y_continuous
      y_limits <- c(0, s_val_max)
    } else {
      y_scale <- scale_y_reverse
      y_limits <- c(s_val_max, 0)
    }

    p <- p +
      y_scale(
        limits = y_limits,
        breaks = scales::pretty_breaks(n = 10)(c(0, s_val_max)),
        sec.axis = sec_axis(
          transform = ~ . + log2(2),
          name = ylab_sec,
          breaks = scales::pretty_breaks(n = 10)(c(
            1,
            -log2(min(res$res_frame$p_two, na.rm = TRUE) / 2)
          ))
        )
      )
  }

  # Y-axis formatting For confidence distributions and densities

  if (plot_type %in% c("cdf", "pdf")) {
    # For the confidence distributions (plot_type == "cdf"), inverted y-axis if specified

    if (plot_type %in% c("cdf") && isTRUE(inverted)) {
      p <- p +
        scale_y_continuous(
          transform = "reverse",
          breaks = scales::pretty_breaks(n = 10)
        )
    } else {
      p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    }
  }

  #-----------------------------------------------------------------------------
  # x-axis
  #-----------------------------------------------------------------------------

  if (x_scale %in% "log") {
    # Plot x-axis on a log-scale (can't set "xlim" here, because then the gray area cannot be added!)
    p <- p +
      scale_x_continuous(
        transform = "log",
        breaks = scales::pretty_breaks(n = 10)
      )

    # If y-axis is plotted on a log-scale, re-add the gray rectangle
    if (plot_type %in% c("p_val") && isTRUE(log_yaxis)) {
      p <- p +
        annotate(
          "rect",
          xmin = 0,
          xmax = 100,
          ymin = ifelse(
            alternative %in% "two_sided",
            plot_p_limit,
            plot_p_limit * 2
          ),
          ymax = cut_logyaxis,
          alpha = 0.1,
          colour = grey(0.9)
        )
    }
  } else {
    p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  }

  # Set x-limits now

  p <- p + coord_cartesian(xlim = xlim, expand = TRUE)

  #-----------------------------------------------------------------------------
  # Horizontal lines at the significance levels if specified
  #-----------------------------------------------------------------------------

  if (plot_type %in% c("p_val", "s_val", "cdf") && !is.null(conf_level)) {
    hlines_tmp <- text_frame$p_value

    if (plot_type %in% "s_val") {
      hlines_tmp <- -log2(hlines_tmp)
    }

    p <- p + geom_hline(yintercept = hlines_tmp, linetype = 2)
  }

  #-----------------------------------------------------------------------------
  # Vertical lines at the null values if specified
  #-----------------------------------------------------------------------------

  if (!is.null(null_values)) {
    x_range <- ggplot_build(p)$layout$panel_params[[1]]$x.range

    if (trans %in% "exp" && x_scale %in% "log") {
      # If the x-axis was log-transformed, we need to backtransforme the plotting limits because they are given on the log-scale
      plot_limits <- trans_fun(x_range)
    } else if (!trans %in% "exp" && x_scale %in% "log") {
      plot_limits <- exp(x_range)
    } else {
      plot_limits <- x_range
    }

    # Which null_values are outside of the plotting limits

    null_outside_plot <- which(
      (res$counternull_frame$null_value <= plot_limits[1]) |
        (res$counternull_frame$null_value >= plot_limits[2])
    )

    # Only add lines for those null values that are inside the plotting limits

    null_lines <- res$counternull_frame

    if (length(null_outside_plot) > 0) {
      # Print a message that shows which null values were outside of the x-axis limits.

      cli::cli_inform(paste0(
        "The following null values are outside of the specified x-axis range (xlim) and are not shown: ",
        paste(
          unique(res$counternull_frame$null_value[null_outside_plot]),
          collapse = ", "
        )
      ))

      null_lines <- null_lines[-null_outside_plot, ]
    }

    if (nrow(null_lines) > 0L) {
      p <- p +
        geom_vline(
          data = null_lines,
          aes(xintercept = .data$null_value),
          linetype = 1,
          linewidth = 0.5
        )
    }

    rm(x_range, plot_limits, null_outside_plot, null_lines)
  }

  #-----------------------------------------------------------------------------
  # Facets for multiple estimates not plotted together
  #-----------------------------------------------------------------------------

  if (length(estimate) >= 2 && isFALSE(together)) {
    if (plot_type %in% "pdf") {
      p <- p +
        facet_wrap(
          vars(.data$variable),
          nrow = nrow,
          ncol = ncol,
          scales = "free"
        )
    } else {
      p <- p +
        facet_wrap(
          vars(.data$variable),
          nrow = nrow,
          ncol = ncol,
          scales = "free_x"
        )
    }
  }

  #-----------------------------------------------------------------------------
  # Add text boxes at the specified significance levels, if any
  #-----------------------------------------------------------------------------

  if (!plot_type %in% "pdf" && !is.null(conf_level)) {
    if (plot_type %in% "s_val") {
      text_frame$p_value <- -log2(text_frame$p_value)
    }

    # Boxes without border; `label.size` was deprecated in ggplot2 4.0.0
    border_args <- if (utils::packageVersion("ggplot2") >= "4.0.0") {
      list(linewidth = 0)
    } else {
      list(label.size = NA)
    }

    p <- p +
      do.call(
        geom_label,
        c(
          list(
            data = text_frame,
            mapping = aes(
              x = .data$theor_values,
              y = .data$p_value,
              label = .data$label
            ),
            inherit.aes = FALSE,
            parse = TRUE,
            size = 5.5,
            hjust = "inward"
          ),
          border_args
        )
      )
  }

  #-----------------------------------------------------------------------------
  # Add points for the counternull if specified
  #-----------------------------------------------------------------------------

  if (
    !is.null(null_values) &&
      isTRUE(plot_counternull) &&
      (plot_type %in% c("p_val", "s_val")) &&
      !all(is.na(res$res_frame$counternull))
  ) {
    if (multi_together && !same_color) {
      p <- p +
        geom_point(
          aes(x = .data$values, y = .data$counternull, colour = .data$variable),
          size = 4,
          shape = 21,
          fill = "white",
          stroke = 1.7
        ) +
        guides(colour = guide_legend(override.aes = list(pch = NA)))
    } else {
      p <- p +
        geom_point(
          aes(x = .data$values, y = .data$counternull),
          colour = col,
          size = 4,
          shape = 21,
          fill = "white",
          stroke = 1.7
        )
    }
  }

  #-----------------------------------------------------------------------------
  # Make the plot prettier by increasing font size
  #-----------------------------------------------------------------------------

  p <- p +
    theme(
      axis.title.y.left = element_text(
        colour = "black",
        size = 17,
        hjust = 0.5,
        margin = margin(0, 10, 0, 0)
      ),
      axis.title.y.right = element_text(
        colour = "black",
        size = 17,
        hjust = 0.5,
        margin = margin(0, 0, 0, 10)
      ),
      axis.title.x = element_text(colour = "black", size = 17),
      axis.text.x = element_text(colour = "black", size = 15),
      axis.text.y = element_text(colour = "black", size = 15),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text.x = element_text(size = 15)
    )

  #-----------------------------------------------------------------------------
  # Add title and labels if specified
  #-----------------------------------------------------------------------------

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  #-----------------------------------------------------------------------------
  # Only print and return the ggplot2-object if requested
  #-----------------------------------------------------------------------------

  if (isTRUE(plot)) {
    res$plot <- p
    suppressWarnings(print(p))
  }

  #-----------------------------------------------------------------------------
  # Sort the data frame for convenience
  #-----------------------------------------------------------------------------

  res$res_frame <- res$res_frame[order(res$res_frame$values), ]

  res
}
