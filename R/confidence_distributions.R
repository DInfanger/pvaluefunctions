#========================================================================
# Construct confidence distributions, densities and p-value functions
# Author: Denis Infanger
# Creation date (dd.mm.yyyy): 22.09.2018
#========================================================================

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("values", "variable", "hypothesis", "null_value", "theor_values", "p_value", "label", "counternull"))
}

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
#' @param tstat Numerical vector contaiqning the \emph{t}-statistic(s). Required for \emph{t}-tests (means and mean differences). Must be equal the number of estimates.
#' @param type String indicating the type of the estimate. Must be one of the following: \code{ttest}, \code{linreg}, \code{gammareg}, \code{general_t}, \code{logreg}, \code{poisreg}, \code{coxreg}, \code{general_z}, \code{pearson}, \code{spearman}, \code{kendall}, \code{var}, \code{prop}, \code{propdiff}.
#' @param plot_type String indicating the type of plot. Must be one of the following: \code{cdf} (confidence distribution), \code{pdf} (confidence density), \code{p_val} (\emph{p}-value function, the default), \code{s_val} (Surprisal value functions). For differences between independent proportions, only \emph{p}-value functions and Surprisal values are available.
#' @param n_values (optional) Integer indicating the number of points that are used to generate the graphics. The higher this number, the higher the computation time and resolution.
#' @param est_names (optional) String vector indicating the names of the estimate(s). Must be equal the number of estimates.
#' @param conf_level (optional) Numerical vector indicating the confidence level(s). Bust be between 0 and 1.
#' @param null_values (optional) Numerical vector indicating the null value(s) in the plot on the \emph{untransformed (original)} scale. For example: The null values for an odds ratio of 1 is 0 on the log-odds scale. If x limits are specified with \code{xlim}, all null values outside of the specified x limits are ignored for plotting and a message is printed.
#' @param trans (optional) String indicating the transformation function that will be applied to the estimates and confidence curves. For example: \code{"exp"} for an exponential transformation of the log-odds in logistic regression. Can be a custom function.
#' @param alternative String indicating if the confidence level(s) are two-sided or one-sided. Must be one of the following: \code{two_sided}, \code{one_sided}.
#' @param log_yaxis Logical. Indicating if a portion of the y-axis should be displayed on the logarithmic scale.
#' @param cut_logyaxis Numerical value indicating the threshold below which the y-axis will be displayed logarithmically. Must lie between 0 and 1.
#' @param xlim (optional) Optional numerical vector of length 2 (x1, x2) indicating the limits of the x-axis on the \emph{untransformed} scale if \code{trans} is not \code{identity}. The scale of the x-axis set by \code{x_scale} does not affect the x limits. For example: If you want to plot \emph{p}-value functions for odds ratios from logistic regressions, the limits have to be given on the log-odds scale if \code{trans = "exp"}. Note that x1 > x2 is allowed but then x2 will be the left limit and x1 the right limit (i.e. the limits are sorted before plotting). Null values (specified in \code{null_values}) that are outside of the specified limits are ignored and a message is printed.
#' @param together Logical. Indicating if graphics for multiple estimates should be displayed together or on separate plots.
#' @param plot_legend Logical. Indicating if a legend should be plotted if multiple curves are plotted together with different colors (i.e. \code{together = TRUE)} and \code{same_color = FALSE}).
#' @param same_color Logical. Indicating if curves should be distinguished using colors if they are plotted together (i.e. \code{together = TRUE}).
#' @param col String indicating the colour of the curves. Only relevant for single curves, multiple curves not plotted together (i.e. \code{together = FALSE}) and multiple curves plotted together but with the option \code{same_color} set to \code{TRUE}.
#' @param nrow (optional) Integer greater than 0 indicating the number of rows when \code{together = FALSE} is specified for multiple estimates. Used in \code{facet_wrap} in ggplot2.
#' @param ncol (optional) Integer greater than 0 indicating the number of columns when \code{together = FALSE} is specified for multiple estimates. Used in \code{facet_wrap} in ggplot2.
#' @param plot_p_limit Numerical value indicating the lower limit of the y-axis. Must be greater than 0 for a logarithmic scale (i.e. \code{log_yaxis = TRUE}). The default is to omit plotting \emph{p}-values smaller than 1 - 0.999 = 0.001.
#' @param plot_counternull Logical. Indicating if the counternull should be plotted as a point. Only available for \emph{p}-value functions and s-value functions. Counternull values that are outside of the plotted functions are not shown.
#' @param title (optional) String containing a title of the plot.
#' @param xlab (optional) String indicating the label of the x-axis.
#' @param ylab (optional) String indicating the title for the primary (left) y-axis.
#' @param ylab_sec (optional) String indicating the title for the secondary (right) y-axis.
#' @param inverted Logical. Indicating the orientation of the y-axis for the \emph{P}-value function (\code{p_val}), S-value function (\code{s_val}) and the confidence distribution (\code{cdf}). By default (i.e. \code{inverted = FALSE}) small \emph{P}-values are plotted at the bottom and large ones at the top so that the cusp of the \emph{P}-value function is a the top. By setting \code{inverted = TRUE}, the y-axis is inverted. Ignored for confidence densities.
#' @param x_scale String indicating the scaling of the x-axis. The default is to scale the x-axis logarithmically if the transformation specified in \code{trans} is "exp" (exponential) and linearly otherwise. The option \code{linear} (can be abbreviated) forces a linear scaling and the option \code{logarithm} (can be abbreviated) forces a logarithmic scaling, regardless what has been specified in \code{trans}.
#' @param plot Logical. Should a plot be created (\code{TRUE}, the default) or not (\code{FALSE}). \code{FALSE} can be useful if users want to create their own plots using the returned data from the function. If \code{FALSE}, no ggplot2 object is returned.

#' @return \code{conf_dist} returns four data frames and if \code{plot = TRUE} was specified, a ggplot2-plot object: \code{res_frame} (contains parameter values (e.g. mean differences, odds ratios etc.), \emph{p}-values (one- and two-sided), s-values, confidence distributions and densities, variable names and type of hypothesis), \code{conf_frame} (contains the specified confidence level(s) and the corresponding lower and upper limits as well as the corresponding variable name), \code{counternull_frame} (contains the counternull and the corresponding null values), \code{point_est} (contains the mean, median and mode point estimates) and if \code{plot = TRUE} was specified, \code{aucc_frame} contains the estimated AUCC (area under the confidence curves) calculated by trapezoidal integration on the untransformed scale, \code{plot} (a ggplot2 object).
#' @references Bender R, Berg G, Zeeb H. Tutorial: using confidence curves in medical research. \emph{Biom J.} 2005;47(2):237-247.
#'
#' Berrar D. Confidence curves: an alternative to null hypothesis significance testing for the comparison of classifiers. \emph{Mach Learn.} 2017;106:911-949.
#'
#' Bonett DG, Wright TA. Sample size requirements for estimating Pearson, Kendall and Spearman correlations. \emph{Psychometrika.} 2000;65(1):23-28.
#'
#' Infanger D, Schmidt-Trucksäss A. \emph{P} value functions: An underused method to present research results and to promote quantitative reasoning. \emph{Stat Med.} 2019;38:4189-4197.
#'
#' Newcombe RG. Interval estimation for the difference between independent proportions: comparison of eleven methods. \emph{Stat Med.} 1998;17:873-890.
#'
#' Poole C. Confidence intervals exclude nothing. \emph{Am J Public Health.} 1987;77(4):492-493.
#'
#' Poole C. Beyond the confidence interval. \emph{Am J Public Health.} 1987;77(2):195-199.
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
#'   100*(1 - exp(x))
#' }
#'
#' rse_fun_inv <- function(x){
#'   log(1 - (x/100))
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
#' @import stats ggplot2 scales
#' @importFrom grDevices grey
#' @importFrom utils find
#' @export

conf_dist <- function(
  estimate = NULL
  , n = NULL
  , df = NULL
  , stderr= NULL
  , tstat = NULL
  , type = NULL
  , plot_type = c("p_val", "s_val", "cdf", "pdf")
  , n_values = 1e4L
  , est_names = NULL
  , conf_level = NULL
  , null_values = NULL
  , trans = "identity"
  , alternative = c("two_sided", "one_sided")
  , log_yaxis = FALSE
  , cut_logyaxis = 0.05
  , xlab = NULL
  , xlim = NULL
  , together = FALSE
  , plot_legend = TRUE
  , same_color = FALSE
  , col = "black"
  , nrow = NULL
  , ncol = NULL
  , plot_p_limit = (1 - 0.999)
  , plot_counternull = FALSE
  , title = NULL
  , ylab = NULL
  , ylab_sec = NULL
  , inverted = FALSE
  , x_scale = c("default", "linear", "logarithm")
  , plot = TRUE
) {

  #-----------------------------------------------------------------------------
  # Load required packages
  #-----------------------------------------------------------------------------

  # require("ggplot2")
  # require(reshape2)
  # require("scales")
  # require(RColorBrewer)
  # require(MASS)
  # require("zipfR")

  #-----------------------------------------------------------------------------
  # Safety checks and clean ups
  #-----------------------------------------------------------------------------

  # alternative <- gsub("[-|.]", "_", alternative)
  # alternative <- tolower(alternative)
  alternative <- match.arg(alternative)
  trans <- tolower(trans)
  type <- tolower(type)
  # plot_type <- gsub("[-|.]", "_", plot_type)
  # plot_type <- tolower(plot_type)
  plot_type <- match.arg(plot_type)
  plot_p_limit <- round(plot_p_limit, 10)
  cut_logyaxis <- round(cut_logyaxis, 10)
  # x_scale <- tolower(x_scale)
  x_scale <- match.arg(x_scale)

  # if (!x_scale %in% "default") {
  #   x_scale <- substr(x_scale, 1, 3) # abbreviate the string to the first 3 letters if not "default"
  # }

  if (is.null(estimate)) {stop("Please provide an estimate.")}

  if (length(type) == 0L) {stop("Please provide the type of the estimate(s).")}

  if (!alternative %in% c("one_sided", "two_sided")) {stop("Alternative must be either \"two_sided\" or \"one_sided\".")}

  if (plot_p_limit == 0 && isTRUE(log_yaxis)) {stop("Cannot plot 0 on logarithmic axis.")}

  if (plot_p_limit >= 0.5 && alternative %in% "one_sided") {stop("Plot limit must be below 0.5 for one-sided hypotheses.")}

  if (any((1 - conf_level) >= 1) && alternative %in% "two_sided") {
    conf_level <- conf_level[-which((1 - conf_level) >= 1)]
  }

  if (any((2 - 2*conf_level) >= 1) && alternative %in% "one_sided") {
    conf_level <- conf_level[-which((2 - 2*conf_level) >= 1)]
  }

  if (type %in% c("pearson", "spearman", "kendall", "var", "prop", "propdiff") && !trans %in% "identity") {
    trans <- "identity"
    cat("\nTransformation changed to identity.\n")
  }

  if (!trans %in% "identity" && length(find(trans)) == 0) {
    stop(paste0("Function ", trans, " was not found."))
  }

  # if (!plot_type %in% c("p_val", "cdf", "pdf", "s_val")) {
  #   stop("plot_type must be one of: p_val, cdf, pdf or s_val.")
  # }

  if (type %in% c("prop", "propdiff") && (any(estimate < 0) || any(estimate > 1))) {
    stop("Please provide proportion estimates as decimals between 0 and 1.")
  }

  if (type %in% "propdiff" && ((length(estimate) != 2L) || (length(n) != 2L))) {
    stop("Please provide exactly two estimates and two sample sizes (n) for a difference in proportions.")
  }

  if (type %in% "propdiff" && !plot_type %in% c("p_val", "s_val")) {
    stop("Currently, only P-value functions (p_val) and S-value functions (s_val) are allowed for difference in proportions.")
  }

  if (type %in% "propdiff" && (((estimate[1]*n[1])%%1 >= 0.05) || ((estimate[2]*n[2])%%1 >= 0.05))) {
    warning("Number of successes (i.e. estimate*n) of proportions not integer! The the number of successes was rounded (i.e. round(estimate*n)).")
  }

  if (!is.null(conf_level) && (any(conf_level <= 0) || any(conf_level >= 1))) {
    stop("All confidence levels must lie between 0 and 1.")
  }

  if (type %in% c("pearson", "spearman", "kendall") && !is.null(estimate) && any(abs(estimate) > 1)) {
    stop("Correlation coefficients must lie between -1 and 1.")
  }

  if (type %in% c("pearson", "spearman", "kendall") && !is.null(null_values) && any(abs(null_values) > 1)) {
    stop("Null values for correlations must lie between -1 and 1.")
  }

  if (type %in% c("prop") && !is.null(null_values) && (any(null_values <= 0) || any(null_values >= 1))) {
    stop("Null values for proportions must lie between 0 and 1 (excluding).")
  }

  if (!is.null(xlab) & (length(xlab) != 1L)) {
    stop("Length of x-axis label must be 1.")
  }

  if (!type %in% "propdiff" && !is.null(est_names) && (length(est_names) != length(estimate))) {
    stop("Length of estimates does not match length of estimate names.")
  }

  if (type %in% "propdiff" && !is.null(est_names) && (length(est_names) > 1L)) {
    stop("Provide only one estimate name for a proportion difference.")
  }

  if (!is.null(xlim)) {
    if ((length(xlim) != 2L)) {
      stop("Please provide two limits for the x-axis.")
    }

    if (any(is.na(xlim)) || any(!is.finite(xlim))) {
      stop("Missing or infinite values are not allowed for x-axis limits (xlim). Please provide exactly two finite x-axis limits.")
    }

    xlim <- sort(xlim, decreasing = FALSE)
  }

  # if ((length(type) == 0L) || (!type %in% c("ttest", "linreg", "gammareg", "general_t", "logreg", "poisreg", "coxreg", "general_z", "pearson", "spearman", "kendall", "var", "prop", "propdiff"))) {
  #   stop("\"type\" must be one of: ttest, linreg, gammareg, general_t, logreg, poisreg, coxreg, general_z, pearson, spearman, kendall, var, prop and propdiff.")
  # }

  if (is.null(est_names)) {
    if (type %in% "propdiff") {
      est_names <- (1L)
    } else {
      est_names <- (1L:length(estimate))
    }
  }

  if (type %in% "ttest" && (is.null(tstat) || is.null(df))){
    stop("Please provide the t-statistic and the degrees of freedom of the t-test.")
  }

  if (type %in% c("linreg", "gammareg", "general_t") && (is.null(df) || is.null(stderr))){
    stop("Please provide the (residual) degrees of freedom and the standard error of the estimates.")
  }

  if (type %in% c("logreg", "poisreg", "coxreg", "general_z") && is.null(stderr)){
    stop("Please provide the standard error of the estimates.")
  }

  if (type %in% c("pearson", "spearman", "kendall", "prop") && is.null(n)){
    stop("Please provide the sample size for correlations and proportions.")
  }

  if ((type %in% c("pearson", "spearman") && any(n < 4)) || (type %in% c("kendall") && any(n < 5))){
    stop("Sample size must be at least 4 for Pearson and Spearman and at least 5 for Kendall's correlation.")
  }

  if (type %in% c("var") && is.null(n)){
    stop("Sample size must be given for variance estimates.")
  }

  if (type %in% "ttest" && all(!is.null(df), !is.null(estimate), !is.null(tstat)) && !identical(length(df), length(estimate), length(tstat))) {
    stop("Degrees of freedom (df) and t-statistics (tstat) must be the same length as estimates.")
  }

  if (type %in% c("linreg", "gammareg", "general_t") && all(!is.null(stderr), !is.null(df), !is.null(estimate)) && !identical(length(stderr), length(df), length(estimate))) {
    stop("Standard errors (stderr) and degrees of freedom (df) must be the same length as estimates.")
  }

  if (type %in% c("coxreg", "logreg", "poisreg", "general_z") && all(!is.null(stderr), !is.null(estimate)) && !identical(length(stderr), length(estimate))) {
    stop("Standard errors (stderr) must be the same length as estimates.")
  }

  if (type %in% c("pearson", "spearman", "kendall", "var", "prop") && all(!is.null(n), !is.null(estimate)) && !identical(length(n), length(estimate))) {
    stop("Sample sizes (n) must be the same length as estimates.")
  }

  if (type %in% c("spearman") && (any(estimate >= 0.9) || any(n < 10))) {
    warning("Approximations for Spearman's correlation are only valid for r < 0.9 and n >= 10. Interpret with caution.")
  }

  if (type %in% c("kendall") && any(estimate >= 0.8)) {
    warning("Approximations for Kendall's correlation are only valid for r < 0.8. Interpret with caution.")
  }

  if (!type %in% "propdiff" && isFALSE(together) && (length(estimate) > 1L) && !is.null(nrow) && !is.null(ncol) && (nrow*ncol < length(estimate))) {
    stop("nrow * ncol must be greater than or equal the number of estimates to be plotted if together = FALSE.")
  }

  # if (!x_scale %in% c("default", "lin", "log")) {
  #   stop("x_scale must be: default, linear (lin), logarithm (log).")
  # }

  if ((trans %in% "exp") && (x_scale %in% "default")) {
    x_scale <- "log"
  }

  if (length(col) > 1) {
    warning(paste0("Only first color is used: ", col[1]))
    col <- col[1]
  }

  #-----------------------------------------------------------------------------
  # Calculate the confidence distributions/densities and p-value curves
  #-----------------------------------------------------------------------------

  if (type %in% "ttest") {

    stderr <- estimate/tstat

    res <- cdist_t(
      estimate = estimate
      , stderr = stderr
      , df = df
      , n_values = n_values
      , conf_level = conf_level
      , alternative = alternative
      , null_values = null_values
    )

  } else if (type %in% c("linreg", "gammareg", "general_t")) {

    res <- cdist_t(
      estimate = estimate
      , stderr = stderr
      , df = df
      , n_values = n_values
      , conf_level = conf_level
      , alternative = alternative
      , null_values = null_values
    )

  } else if (type %in% c("logreg", "poisreg", "coxreg", "general_z")) {

    res <- cdist_z(
      estimate = estimate
      , stderr = stderr
      , n_values = n_values
      , conf_level = conf_level
      , null_values = null_values
      , alternative = alternative
    )

  } else if (type %in% c("pearson", "spearman", "kendall")) {

    # Calculate approximate standard error for each type of correlation coefficient
    # Ref 1: Bonett & Wright (2000): Sample size requirements for estimating Pearson, Kendall and Spearman correlations
    # Ref 2: Fieller, Hartley, Pearson (1957): Tests for rank correlation coefficients I.

    stderr <- switch(
      type
      , pearson = 1/sqrt(n - 3)
      , spearman = sqrt((1 + (estimate)^2/2)/(n - 3))
      , kendall = sqrt(0.437/(n - 4))
    )

    res <- cdist_corr(
      estimate = estimate
      , stderr = stderr
      , n = n
      , n_values = n_values
      , conf_level = conf_level
      , null_values = null_values
      , alternative = alternative
    )

  } else if (type %in% "var") {

    res <-  cdist_var(
      estimate = estimate
      , n = n
      , n_values = n_values
      , conf_level = conf_level
      , null_values = null_values
      , alternative = alternative
    )

  } else if (type %in% "prop") {

    res <-  cdist_prop1(
      estimate = estimate
      , n = n
      , n_values = n_values
      , conf_level = conf_level
      , null_values = null_values
      , alternative = alternative
    )

  } else if (type %in% "propdiff") {

    res <-  cdist_propdiff(
      estimate = estimate
      , n = n
      , n_values = n_values
      , conf_level = conf_level
      , null_values = null_values
      , alternative = alternative
    )

    estimate <- estimate[1] - estimate[2]

  }

  #-----------------------------------------------------------------------------
  # Calculate Shannon-surprisal value (S-value)
  #-----------------------------------------------------------------------------

  res$res_frame$s_val <- -log2(res$res_frame$p_two)

  #-----------------------------------------------------------------------------
  # Assign estimate names
  #-----------------------------------------------------------------------------

  if (!is.null(est_names)) {

    res$point_est$variable <- factor(res$point_est$variable, labels = est_names)
    res$res_frame$variable <- factor(res$res_frame$variable, labels = est_names)

    if (!is.null(conf_level)) {
      res$conf_frame$variable <- factor(res$conf_frame$variable, labels = est_names)
    }
    if (!is.null(null_values)) {
      res$counternull_frame$variable <- factor(res$counternull_frame$variable, labels = est_names)
    }
  }

  #-----------------------------------------------------------------------------
  # Calculate AUCC (area under the confidence curve), see Berrar (2017) Mach Learn 106:911-494
  #-----------------------------------------------------------------------------

  res$aucc_frame <- data.frame(
    variable = est_names
    , aucc = NA
  )

  for(i in seq_along(estimate)) {

    x_tmp <- res$res_frame$values[res$res_frame$variable %in% est_names[i]]
    y_tmp <- res$res_frame$p_two[res$res_frame$variable %in% est_names[i]]

    nona_ind <- which(!is.na(y_tmp) & !is.na(x_tmp))

    order_tmp <- order(res$res_frame$values[res$res_frame$variable %in% est_names[i]][nona_ind], decreasing = FALSE)

    res$aucc_frame$aucc[res$aucc_frame$variable %in% est_names[i]] <- pracma::trapz(
      x = x_tmp[nona_ind][order_tmp]
      , y = y_tmp[nona_ind][order_tmp]
    )

    rm(x_tmp, y_tmp, nona_ind, order_tmp)

  }

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
      res$res_frame$hypothesis[res$res_frame$variable %in% est_names[i] & res$res_frame$values < estimate[i]] <- 1     # greater
      res$res_frame$hypothesis[res$res_frame$variable %in% est_names[i] & res$res_frame$values >= estimate[i]] <- (-1) # less
    }
  }

  res$res_frame$hypothesis <- factor(res$res_frame$hypothesis, levels = c(-1, 1), labels = c("less", "greater"))

  #-----------------------------------------------------------------------------
  # Calculate the limits of the x-axis if not provided
  #-----------------------------------------------------------------------------

  # If there are limits given, take those in any case
  if(is.null(xlim)) {

    xlim <- c(NA, NA)

    if (is.null(null_values)) {

      # If no limits and no null values given, take the plot_p_limit
      res_tmp <- res$res_frame
      res_tmp$values[res$res_frame$p_two < plot_p_limit] <- NA

      xlim <- range(res_tmp$values, na.rm = TRUE)

      rm(res_tmp)

    } else if (!is.null(null_values)) {

      # If no limits but null values given, look if the plot_p_limits are outside the null_values
      res_tmp <- res$res_frame
      res_tmp$values[res$res_frame$p_two < plot_p_limit] <- NA

      plot_range_tmp <- range(res_tmp$values, na.rm = TRUE)

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

      rm(res_tmp, plot_range_tmp)

    }
  }

  # if (is.null(xlim) & plot_type %in% c("p_val", "s_val")) {
  #   res_tmp <- res$res_frame
  #   res_tmp$values[res$res_frame$p_two < plot_p_limit] <- NA # Set all values were the p-value is below the specified value to missing
  #   # res_tmp$p_two[res$res_frame$p_two < plot_p_limit] <- NA
  #
  #   xlim <- range(res_tmp$values, na.rm = TRUE)
  #
  #   rm(res_tmp)
  # }

  #-----------------------------------------------------------------------------
  # Text frame coordinates and contents for plotting the confidence levels
  #-----------------------------------------------------------------------------

  if (!is.null(conf_level)) {
    if ((isTRUE(together) & (length(estimate) >= 2)) | plot_type %in% "cdf") {

      min_theor_values <- min(tapply(res$res_frame$values, res$res_frame$variable, min, na.rm = TRUE), na.rm = TRUE)
      max_theor_values <- max(tapply(res$res_frame$values, res$res_frame$variable, max, na.rm = TRUE), na.rm = TRUE)

      theor_val_tmp <- switch(
        alternative
        , two_sided = rep(ifelse(trans %in% "exp", 0, -Inf), each = length(conf_level))
        , one_sided = rep(Inf, each = length(conf_level))
      )

      p_val_tmp <- switch(
        alternative
        , two_sided = round((1 - conf_level), 10)
        , one_sided = round((2 - 2*conf_level), 10)
      )

      text_frame <- data.frame(
        label = (1 - conf_level)
        , theor_values = theor_val_tmp
        , p_value = p_val_tmp
        # , variable = factor(rep(est_names, length(conf_level)))
        # , variable = rep(1, length(conf_level))
      )

      # if (!is.null(xlim)) {
      #   text_frame$theor_values <- ifelse(alternative %in% "two_sided", min(xlim), max(xlim))
      # }

    } else if (isFALSE(together) | (isTRUE(together) & (length(estimate) < 2))) {

      min_theor_values <- tapply(res$res_frame$values, res$res_frame$variable, min, na.rm = TRUE)
      max_theor_values <- tapply(res$res_frame$values, res$res_frame$variable, max, na.rm = TRUE)

      theor_val_tmp <- switch(
        alternative
        , two_sided = rep(ifelse(trans %in% "exp", 0, -Inf), each = length(conf_level))
        , one_sided = rep(Inf, each = length(conf_level))
      )

      p_val_tmp <- switch(
        alternative
        , two_sided = round((1 - conf_level), 10)
        , one_sided = round((2 - 2*conf_level), 10)
      )

      text_frame <- data.frame(
        label = (1 - conf_level)
        , theor_values = theor_val_tmp
        , p_value = p_val_tmp
        # , variable = factor(rep(est_names, each = length(conf_level)))
      )

    }
  }

  #-----------------------------------------------------------------------------
  # Apply transformations if applicable
  #-----------------------------------------------------------------------------

  if (!trans %in% "identity") {

    res$res_frame$values <- do.call(trans, list(x = res$res_frame$values))

    res$point_est[, c("est_mean", "est_median", "est_mode")] <- do.call(trans, list(x = res$point_est[, c("est_mean", "est_median", "est_mode")]))

    if (!is.null(conf_level)) {
      if (!trans %in% "exp") {
        text_frame$theor_values[is.finite(text_frame$theor_values)] <- do.call(trans, list(x = text_frame$theor_values[is.finite(text_frame$theor_values)]))
      }
      res$conf_frame$lwr <- do.call(trans, list(x = res$conf_frame$lwr))
      res$conf_frame$upr <- do.call(trans, list(x = res$conf_frame$upr))
    }

    if (!is.null(null_values)) {
      res$counternull_frame$counternull <- do.call(trans, list(x = res$counternull_frame$counternull))
      res$counternull_frame$null_value <- do.call(trans, list(x = res$counternull_frame$null_value))
    }

    xlim <- sort(do.call(trans, list(x = xlim)), decreasing = FALSE)

  }

  #-----------------------------------------------------------------------------
  # Cutoff for nicer plotting
  #-----------------------------------------------------------------------------

  p_cutoff <- ifelse(alternative %in% c("two_sided"), plot_p_limit, plot_p_limit*2)

  if (plot_type %in% c("p_val")) {
    res$res_frame$values[res$res_frame$p_two < p_cutoff] <- NA
    res$res_frame$p_two[res$res_frame$p_two < p_cutoff] <- NA
    res$res_frame$p_one[res$res_frame$p_one < p_cutoff] <- NA
  }

  if (plot_type %in% c("s_val")) {

    # outside_ind <- which(res$res_frame$values < min(xlim) | res$res_frame$values > max(xlim))
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

  if (!is.null(null_values) && isTRUE(plot_counternull) && plot_type %in% c("p_val", "s_val")) {

    res$res_frame$counternull <- NA

    null_values_trans <- do.call(trans, list(x = null_values))

    for (i in seq_along(estimate)) {

      plot_range_tmp <- range(res$res_frame$values[res$res_frame$variable %in% est_names[i]], na.rm = TRUE)

      for (j in seq_along(null_values)) {

        counternull_tmp <- res$counternull_frame$counternull[res$counternull_frame$variable %in% est_names[i] & res$counternull_frame$null_value %in% null_values_trans[j]]

        if ((counternull_tmp <= plot_range_tmp[2]) & (counternull_tmp >= plot_range_tmp[1])) {

          counternull_index_tmp <- which.min(abs(res$res_frame$values[res$res_frame$variable %in% est_names[i]] - counternull_tmp))

          res$res_frame$counternull[res$res_frame$variable %in% est_names[i]][counternull_index_tmp] <- res$res_frame$p_two[res$res_frame$variable %in% est_names[i]][which.min(abs(res$res_frame$values[res$res_frame$variable %in% est_names[i]] - null_values_trans[j]))]

        }
      }
    }

    if (plot_type %in% "s_val") {
      res$res_frame$counternull <- -log2(res$res_frame$counternull)
    }
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
    cut_logyaxis <- cut_logyaxis_one*2
  } else {
    cut_logyaxis_one <- cut_logyaxis/2
  }

  # Labeller functions for custom log-scale

  lab_onesided <- Vectorize(function(x){
    if(!is.na(x) && (x < cut_logyaxis_one) & (round((x %% 1)*10) == 0)) {sprintf("%.5g", x)}
    else {sprintf("%.2f", x)}
  })

  lab_twosided <- Vectorize(function(x){
    if(!is.na(x) && (x <= cut_logyaxis) & (round((x %% 1)*10) == 0)) {sprintf("%.5g", x)}
    else {sprintf("%.1f", x)}
  })

  # Start plotting

  theme_set(theme_bw())

  # Set variable to be plotted on y-axis depending on the type of the plot

  y_var <- switch(
    plot_type
    , p_val = "p_two"
    , cdf = "conf_dist"
    , pdf = "conf_dens"
    , s_val = "s_val"
  )

  # Set label of the x-axis if not provided

  if (is.null(xlab)) {
    xlab <- "Estimate"
  }

  # Set label of the primary and secondary y-axis depending on the type of the plot

  if (is.null(ylab)) {
    ylab <- switch(
      plot_type
      , p_val = expression(paste(italic("P"), "-value (two-sided) / Significance level"~alpha, sep = ""))
      , cdf = "Confidence distribution"
      , pdf = "Confidence density"
      , s_val = expression(paste("Surprisal in bits (two-sided ",~italic("P"), "-value)", sep = ""))
    )
  }

  if (is.null(ylab_sec)) {

    ylab_sec <- switch(
      plot_type
      , p_val = expression(paste(italic("P"), "-value (one-sided) / Significance level"~alpha, sep = ""))
      , s_val = expression(paste("Surprisal in bits (one-sided ",~italic("P"), "-value)", sep = ""))
    )
  }

  # Create a ggplot2-object with no geoms

  p <- ggplot(res$res_frame, aes(x = values, y = eval(parse(text = y_var)), group = variable))

  # If 2 or more estimates are plotted together, differentiate them by color (if user did not specify "same_color = TRUE")

  if ((length(estimate) >= 2) & isTRUE(together) & isFALSE(same_color)) {
    p <- p + aes(colour = variable)
  }

  # For only one one-sided p-value curve, set the colors to black and blue

  if (alternative %in% "one_sided" & (isFALSE(together) | (isTRUE(together) & length(estimate) < 2)) & isFALSE(same_color)) {
    p <- p + geom_line(aes(colour = hypothesis), size = 1.5) +
      scale_colour_manual(values = c("black", "#08A9CF"))  +
      theme(
        legend.position="none"
      )

  } else if (alternative %in% "one_sided" & (isFALSE(together) | (isTRUE(together) & length(estimate) < 2)) & isTRUE(same_color)) {

    p <- p + geom_line(aes(colour = hypothesis), size = 1.5) +
      scale_colour_manual(values = c(col, col))  +
      theme(
        legend.position="none"
      )

    # For only one two-sided p-value curve, set the color to black

  } else if (alternative %in% "two_sided" & (isFALSE(together) | (isTRUE(together) & length(estimate) < 2) | (isTRUE(together) & isTRUE(same_color)))) {
    p <- p + geom_line(size = 1.5, colour = col)

    # For 2 or more estimates plotted together: set the colors according to "Set1" palette

  } else if ((alternative %in% "two_sided" & (length(estimate) >= 2) & isTRUE(together) & isFALSE(same_color) & isTRUE(plot_legend)) |
             (alternative %in% "one_sided" & (length(estimate) >= 2) & isTRUE(together) & isTRUE(plot_legend))) { # Plot legend
    p <- p + geom_line(size = 1.5) +
      scale_colour_brewer(palette = "Set1", name = "") +
      theme(
        legend.position="top"
        , legend.text=element_text(size=15)
        , legend.title=element_text(size=15)
      )

  } else if((alternative %in% "two_sided" & (length(estimate) >= 2) & isTRUE(together) & isFALSE(same_color) & isFALSE(plot_legend)) |
            (alternative %in% "one_sided" & (length(estimate) >= 2) & isTRUE(together) & isFALSE(plot_legend))) { # Plot no legend
    p <- p + geom_line(size = 1.5) +
      scale_colour_brewer(palette = "Set1", name = "") +
      theme(
        legend.position="none"
      )

  }

  # Add the labels for the axes

  p <- p + xlab(xlab) +
    ylab(ylab)

  #-----------------------------------------------------------------------------
  # y-axis
  #-----------------------------------------------------------------------------

  # For p-value curves: Set the left and right y-axes, possibly with a logarithmic part

  if (plot_type %in% "p_val") {
    if (isTRUE(log_yaxis) & (p_cutoff < cut_logyaxis)) {

      lower_ylim_two <- round(10^(ceiling(round(log10(ifelse(alternative %in% "two_sided", plot_p_limit, plot_p_limit*2)), 5))), 10)
      lower_ylim_one <- round(10^(ceiling(round(log10(ifelse(alternative %in% "two_sided", plot_p_limit/2, plot_p_limit)), 5))), 10)

      if ((alternative %in% "two_sided" & (lower_ylim_two <= cut_logyaxis)) | (alternative %in% "one_sided" & (lower_ylim_one <= cut_logyaxis*2))) {

        # Split the breaks into two parts: i) below the cutoff i.e. the logarithmic part and ii) the linear part above the cutoff

        breaks_two <- c(10^(seq(log10(lower_ylim_two), log10(cut_logyaxis), by = 1)), seq(ceiling(cut_logyaxis/0.1)*0.1, 1, by = 0.1))
        breaks_one <- c(10^(seq(log10(lower_ylim_one), log10(cut_logyaxis_one), by = 1)), seq(ceiling(cut_logyaxis_one/0.05)*0.05, 0.5, 0.05))
      } else {
        breaks_two <- c(10^(log10(seq(ceiling(cut_logyaxis/0.1)*0.1, 1, by = 0.1))))
        breaks_one <- c(10^(log10(seq(ceiling(cut_logyaxis_one/0.05)*0.05, 0.5, 0.05))))
      }

      # Remove possible duplicates

      breaks_two <- unique(breaks_two)
      breaks_one <- unique(breaks_one)

      if (isTRUE(inverted)) {

        p <- p +
          scale_y_continuous(
            limits = rev(c(ifelse(alternative %in% "two_sided", plot_p_limit, plot_p_limit*2), 1))
            , breaks = breaks_two
            , labels = lab_twosided
            , trans = magnify_trans_log_rev(interval_low = cut_logyaxis, interval_high = 1, reducer = cut_logyaxis, reducer2 = 8)
            , sec.axis = sec_axis(
              trans = ~.*(1/2)
              , name = ylab_sec
              , breaks = breaks_one
              , labels = lab_onesided
            )
          ) +
          annotate("rect", xmin=-Inf, xmax=Inf, ymin=ifelse(alternative %in% "two_sided", plot_p_limit, plot_p_limit*2), ymax=cut_logyaxis, alpha=0.1, colour = grey(0.9))

      } else {
        p <- p +
          scale_y_continuous(
            limits = c(ifelse(alternative %in% "two_sided", plot_p_limit, plot_p_limit*2), 1)
            , breaks = breaks_two
            , labels = lab_twosided
            , trans = magnify_trans_log(interval_low = cut_logyaxis, interval_high = 1, reducer = cut_logyaxis, reducer2 = 8)
            , sec.axis = sec_axis(
              trans = ~.*(1/2)
              , name = ylab_sec
              , breaks = breaks_one
              , labels = lab_onesided
            )
          ) +
          annotate("rect", xmin=-Inf, xmax=Inf, ymin=ifelse(alternative %in% "two_sided", plot_p_limit, plot_p_limit*2), ymax=cut_logyaxis, alpha=0.1, colour = grey(0.9))
      }

    } else if (isFALSE(log_yaxis) | (p_cutoff >= cut_logyaxis)) {

      if (isTRUE(inverted)) {
        p <- p +
          scale_y_reverse(
            # limits = c(plot_p_limit, 1)
            breaks = seq(0, 1, 0.1)
            , sec.axis = sec_axis(
              ~.*(1/2)
              , name = ylab_sec
              , breaks = seq(0, 1, 0.1)/2
            )
          )
      } else {
        p <- p +
          scale_y_continuous(
            # limits = c(plot_p_limit, 1)
            breaks = seq(0, 1, 0.1)
            , sec.axis = sec_axis(
              ~.*(1/2)
              , name = ylab_sec
              , breaks = seq(0, 1, 0.1)/2
            )
          )
      }
    }
  }

  # For s-value curves: inverted the y-axis and transform it according to log2

  if (plot_type %in% "s_val") {

    if (isTRUE(inverted)) {
      p <- p + scale_y_continuous(
        limits = rev(c(max(res$res_frame$s_val, na.rm = TRUE), 0))
        , breaks = scales::pretty_breaks(n = 10)(c(0, max(res$res_frame$s_val, na.rm = TRUE)))
        , sec.axis = sec_axis(
          trans = ~. + log2(2)
          , name = ylab_sec
          , breaks = scales::pretty_breaks(n = 10)(c(1, -log2(min(res$res_frame$p_two, na.rm = TRUE)/2)))
        )
      )
    } else {
      p <- p + scale_y_reverse(
        limits = c(max(res$res_frame$s_val, na.rm = TRUE), 0)
        , breaks = scales::pretty_breaks(n = 10)(c(0, max(res$res_frame$s_val, na.rm = TRUE)))
        , sec.axis = sec_axis(
          trans = ~. + log2(2)
          , name = ylab_sec
          , breaks = scales::pretty_breaks(n = 10)(c(1, -log2(min(res$res_frame$p_two, na.rm = TRUE)/2)))
        )
      )
    }
  }

  # Y-axis formatting For confidence distributions and densities

  if (plot_type %in% c("cdf", "pdf")) {

    # For the confidence distributions (plot_type == "cdf"), inverted y-axis if specified

    if (plot_type %in% c("cdf") && isTRUE(inverted)) {
      p <- p + scale_y_continuous(trans = "reverse", breaks = scales::pretty_breaks(n = 10))
    } else {
      p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    }
  }

  #-----------------------------------------------------------------------------
  # x-axis
  #-----------------------------------------------------------------------------

  # if (trans %in% "exp") {
  if (x_scale %in% "log") {

    # Plot x-axis on a log-scale (can't set "xlim" here, because then the gray area cannot be added!)
    p <- p + scale_x_continuous(trans = "log", breaks = scales::pretty_breaks(n = 10))

    # If y-axis is plotted on a log-scale, re-add the gray rectangle
    if (plot_type %in% c("p_val") && isTRUE(log_yaxis)) {
      p <- p + annotate(
        "rect"
        , xmin = 0
        , xmax = 100
        , ymin = ifelse(alternative %in% "two_sided", plot_p_limit, plot_p_limit*2)
        , ymax = cut_logyaxis
        , alpha = 0.1
        , colour = grey(0.9)
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

    if (trans %in% "exp" && x_scale %in% "log") { # If the x-axis was log-transformed, we need to backtransforme the plotting limits because they are given on the log-scale
      plot_limits <- do.call(trans, list(x = ggplot_build(p)$layout$panel_params[[1]]$x.range))
    } else if (!trans %in% "exp" && x_scale %in% "log") {
      plot_limits <- exp(ggplot_build(p)$layout$panel_params[[1]]$x.range)
    } else {
      plot_limits <- ggplot_build(p)$layout$panel_params[[1]]$x.range
    }

    # Which null_values are outside of the plotting limits

    null_outside_plot <- which((res$counternull_frame$null_value <= plot_limits[1]) | (res$counternull_frame$null_value >= plot_limits[2]))

    # Only add lines for those null values that are inside the plotting limits

    if (length(null_outside_plot) > 0) {

      # Print a message that shows which null values were outside of the x-axis limits.

      message(paste0("The following null values are outside of the specified x-axis range (xlim) and are not shown: ", paste(unique(res$counternull_frame$null_value[null_outside_plot]), collapse = ", ")))

      if ((length(null_outside_plot) < length(null_values))) { # There are some null-values to plot

        p <- p + geom_vline(data = res$counternull_frame, aes(xintercept = null_value), linetype = 1, size = 0.5)

      } else if ((length(null_outside_plot) == length(null_values))) { # All null-values outside of plotting area

        p <- p + geom_vline(data = res$counternull_frame[-null_outside_plot, ], aes(xintercept = null_value), linetype = 1, size = 0.5)

      }

    } else if (length(null_outside_plot) == 0) { # No null-values outside of plot

      p <- p + geom_vline(data = res$counternull_frame, aes(xintercept = null_value), linetype = 1, size = 0.5)

    }

    rm(plot_limits, null_outside_plot)
  }

  # if (trans %in% "exp" && plot_type %in% "p_val") {
  #
  #   xlim_new <- xlim
  #
  #   # curr_x_limits <- ggplot_build(p)$layout$panel_params[[1]]$x.range
  #
  #   p <- p + scale_x_continuous(trans = "log", breaks = scales::pretty_breaks(n = 10))
  #
  #   if (log_yaxis == TRUE) {
  #     p <- p + annotate("rect", xmin=0, xmax=100, ymin=ifelse(alternative %in% "two_sided", plot_p_limit, plot_p_limit*2), ymax=cut_logyaxis, alpha=0.1, colour = grey(0.9))
  #   }
  #
  #   if (exists("null_values_trans") && xlim[1] <= min(null_values_trans)) {
  #     xlim_new[1] <- xlim[1]
  #   } else if (exists("null_values_trans") && xlim[1] > min(null_values_trans)) {
  #     xlim_new[1] <- min(res$res_frame$values, na.rm = TRUE)
  #   }
  #
  #   if (exists("null_values_trans") && xlim[2] >= max(null_values_trans)) {
  #     xlim_new[2] <- xlim[2]
  #   } else if (exists("null_values_trans") && xlim[2] < max(null_values_trans)) {
  #     xlim_new[2] <- max(res$res_frame$values, na.rm = TRUE)
  #   }
  #
  #   p <- p + coord_cartesian(xlim = xlim_new, expand = TRUE)
  #
  #   if (exists("hlines_tmp")) {
  #     p <- p + geom_hline(yintercept = hlines_tmp, linetype = 2)
  #   }
  #
  # } else if (trans %in% "exp" && plot_type %in% c("cdf", "pdf", "s_val")){
  #
  #   p <- p +
  #     scale_x_continuous(trans = "log", breaks = scales::pretty_breaks(n = 10), limits = xlim)
  #
  # } else {
  #
  #   xlim_new <- xlim
  #
  #   if (!is.null(null_values) && xlim[1] <= min(null_values)) {
  #     xlim_new[1] <- xlim[1]
  #   } else if (!is.null(null_values) && xlim[1] > min(null_values)) {
  #     xlim_new[1] <- min(res$res_frame$values, na.rm = TRUE)
  #   }
  #
  #   if (!is.null(null_values) && xlim[2] >= max(null_values)) {
  #     xlim_new[2] <- xlim[2]
  #   } else if (!is.null(null_values) && xlim[2] < max(null_values)) {
  #     xlim_new[2] <- max(res$res_frame$values, na.rm = TRUE)
  #   }
  #
  #   p <- p +
  #     scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = xlim_new)
  #
  # }
  # }

  #-----------------------------------------------------------------------------
  # Facets for multiple estimates not plotted together
  #-----------------------------------------------------------------------------

  if (length(estimate) >= 2 & isFALSE(together)) {
    if (plot_type %in% "pdf"){
      p <- p + facet_wrap(vars(variable), nrow = nrow, ncol = ncol, scales = "free")
    } else {
      p <- p + facet_wrap(vars(variable), nrow = nrow, ncol = ncol, scales = "free_x")
    }
  }

  #-----------------------------------------------------------------------------
  # Add text boxes at the specified significance levels, if any
  #-----------------------------------------------------------------------------

  if (!plot_type %in% "pdf" && !is.null(conf_level)) {
    if (plot_type %in% "s_val") {
      text_frame$p_value <- -log2(text_frame$p_value)
    }

    p <- p + geom_label(
      data = text_frame
      , mapping = aes(x = theor_values, y = p_value, label = label)
      , inherit.aes = FALSE
      , label.size = NA
      , parse = TRUE
      , size = 5.5
      , hjust = "inward"
    )
  }

  #-----------------------------------------------------------------------------
  # Add points for the counternull if specified
  #-----------------------------------------------------------------------------

  if (!is.null(null_values) && isTRUE(plot_counternull) && (plot_type %in% c("p_val", "s_val")) && !all(is.na(res$res_frame$counternull))) {

    if (isTRUE(together) & (length(estimate) >= 2) & isFALSE(same_color)) {

      p <- p + geom_point(aes(x = values, y = counternull, colour = variable), size = 4, pch = 21, fill = "white", stroke = 1.7) +
        guides(colour = guide_legend(override.aes = list(pch = NA)))

    } else if (isTRUE(together) & (length(estimate) >= 2) & isTRUE(same_color)) {

      p <- p + geom_point(aes(x = values, y = counternull), colour = col, size = 4, pch = 21, fill = "white", stroke = 1.7)

    } else if (isFALSE(together) | (isTRUE(together) & (length(estimate) < 2))) {

      p <- p + geom_point(aes(x = values, y = counternull), colour = col, size = 4, pch = 21, fill = "white", stroke = 1.7)

    }
  }

  #-----------------------------------------------------------------------------
  # Make the plot prettier by increasing font size
  #-----------------------------------------------------------------------------

  p <- p + theme(
    axis.title.y.left=element_text(colour = "black", size = 17, hjust = 0.5, margin = margin(0, 10, 0, 0))
    , axis.title.y.right=element_text(colour = "black", size = 17, hjust = 0.5, margin = margin(0, 0, 0, 10))
    , axis.title.x=element_text(colour = "black", size = 17)
    # , axis.title.y=element_text(size=15,hjust=0.5, vjust=1)
    , axis.text.x=element_text(colour = "black", size=15)
    , axis.text.y=element_text(colour = "black", size=15)
    # , plot.margin=unit(c(2,2,2,2,2),"line")
    , panel.grid.minor.y = element_blank()
    # , panel.grid.major = element_line(colour=grey(0.8), size=0.5)
    , plot.title = element_text(face = "bold")
    # , strip.background=element_rect(fill="white")
    , strip.text.x=element_text(size=15)
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

  return(res)

}


cdist_t <- function(
  estimate = NULL
  , stderr = NULL
  , df = NULL
  , n_values = NULL
  , conf_level = NULL
  , null_values = NULL
  , alternative = NULL
){

  eps <- 1e-10

  res_mat <- matrix(NA, nrow = 0, ncol = 6)

  conf_mat <- matrix(NA, nrow = 0, ncol = 4)

  counternull_mat <- matrix(NA, nrow = 0, ncol = 3)

  for (i in seq_along(estimate)) {

    limits <- c(qt(eps, df = df[i])*(stderr[i]) + estimate[i], qt(1 - eps, df = df[i])*(stderr[i]) + estimate[i])

    x_calc <- c(estimate[i], null_values, seq(limits[1], limits[2], length.out = n_values))

    res_mat_tmp <- matrix(NA, nrow = length(x_calc), ncol = 6)

    res_mat_tmp[, 1] <- x_calc
    res_mat_tmp[, 2] <- pt((x_calc - estimate[i])/stderr[i], df = df[i])
    # res_mat_tmp[, 3] <- t_dens(x = x_calc, estimate = estimate[i], df = df[i], stderr = stderr[i])
    res_mat_tmp[, 3] <- dt((x_calc - estimate[i])/stderr[i], df = df[i])*(1/stderr[i])
    res_mat_tmp[, 4] <- 1 - 2*abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 5] <- (1/2) - abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 6] <- rep(i, length(x_calc))

    res_mat <- rbind(res_mat, res_mat_tmp)

    # Confidence intervals

    if (!is.null(conf_level)) {

      quants_tmp <- switch(
        alternative
        , two_sided = c(1 - (conf_level + 1)/2, (conf_level + 1)/2)
        , one_sided = c((1 - conf_level), conf_level)
      )

      # quants_tmp <- c(1 - (conf_level + 1)/2, (conf_level + 1)/2)

      conf_mat_tmp <- matrix(NA, ncol = 4, nrow = length(conf_level))

      limits_tmp <- matrix(qt(quants_tmp, df = df[i])*stderr[i] + estimate[i], ncol = 2)

      conf_mat_tmp[, 1] <- conf_level
      conf_mat_tmp[, 2] <- limits_tmp[, 1]
      conf_mat_tmp[, 3] <- limits_tmp[, 2]
      conf_mat_tmp[, 4] <- rep(i, length(conf_level))

      conf_mat <- rbind(conf_mat, conf_mat_tmp)

    }

    # Counternulls

    if (!is.null(null_values)) {

      counternull_mat_tmp <- matrix(NA, ncol = 3, nrow = length(null_values))

      counternull_mat_tmp[, 1] <- null_values
      counternull_mat_tmp[, 2] <- 2*estimate[i] - null_values
      counternull_mat_tmp[, 3] <- rep(i, length(null_values))

      counternull_mat <- rbind(counternull_mat, counternull_mat_tmp)

    }

  }

  res_frame <- as.data.frame(res_mat)

  names(res_frame) <- c("values", "conf_dist", "conf_dens", "p_two", "p_one", "variable")

  if (!is.null(conf_level)) {

    conf_frame <- as.data.frame(conf_mat)

    names(conf_frame) <- c("conf_level", "lwr", "upr", "variable")

  } else {
    conf_frame <- NULL
  }

  if (!is.null(null_values)) {

    counternull_frame <- as.data.frame(counternull_mat)

    names(counternull_frame) <- c("null_value", "counternull", "variable")

  } else {
    counternull_frame <- NULL
  }

  # Point estimators

  point_est_frame <- data.frame(
    est_mean = rep(NA, length(estimate))
    , est_median = rep(NA, length(estimate))
    , est_mode = rep(NA, length(estimate))
    , variable = seq(1, i, 1)
  )

  mean_fun <- function(x, df, estimate, stderr){
    x*dt(((x - estimate)/stderr), df = df)*(1/stderr)
  }

  for (i in seq_along(estimate)) {

    point_est_frame$est_mean[i] <- integrate(mean_fun, lower = -Inf, upper = Inf, df = df[i], estimate = estimate[i], stderr = stderr[i], rel.tol = 1e-10)$value # Mean
    point_est_frame$est_median[i] <- res_frame$values[res_frame$variable %in% i][which.min(abs(res_frame$conf_dist[res_frame$variable %in% i][-1] - 0.5)) + 1] # Median
    point_est_frame$est_mode[i] <- res_frame$values[res_frame$variable %in% i][which.max(res_frame$conf_dens[res_frame$variable %in% i])] # Mode

  }

  return(list(
    res_frame = res_frame
    , conf_frame = conf_frame
    , counternull_frame = counternull_frame
    , point_est = point_est_frame
  ))

}

cdist_z <- function(
  estimate = NULL
  , stderr = NULL
  , n_values = NULL
  , conf_level = NULL
  , null_values = NULL
  , alternative = NULL
){

  eps <- 1e-10

  res_mat <- matrix(NA, nrow = 0, ncol = 6)

  conf_mat <- matrix(NA, nrow = 0, ncol = 4)

  counternull_mat <- matrix(NA, nrow = 0, ncol = 3)

  for (i in seq_along(estimate)) {

    limits <- c(qnorm(eps)*stderr[i] + estimate[i], qnorm(1 - eps)*stderr[i] + estimate[i])

    x_calc <- c(estimate[i], null_values, seq(limits[1], limits[2], length.out = n_values))

    res_mat_tmp <- matrix(NA, nrow = length(x_calc), ncol = 6)

    res_mat_tmp[, 1] <- x_calc
    res_mat_tmp[, 2] <- pnorm((x_calc - estimate[i])/stderr[i])
    # res_mat_tmp[, 3] <- z_dens(x = x_calc, estimate = estimate[i], stderr = stderr[i])
    res_mat_tmp[, 3] <- dnorm((x_calc - estimate[i])/stderr[i])*(1/stderr[i])
    res_mat_tmp[, 4] <- 1 - 2*abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 5] <- (1/2) - abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 6] <- rep(i, length(x_calc))

    res_mat <- rbind(res_mat, res_mat_tmp)

    # Confidence intervals

    if (!is.null(conf_level)) {

      quants_tmp <- switch(
        alternative
        , two_sided = c(1 - (conf_level + 1)/2, (conf_level + 1)/2)
        , one_sided = c((1 - conf_level), conf_level)
      )

      # quants_tmp <- c(1 - (conf_level + 1)/2, (conf_level + 1)/2)

      conf_mat_tmp <- matrix(NA, ncol = 4, nrow = length(conf_level))

      limits_tmp <- matrix(qnorm(quants_tmp)*stderr[i] + estimate[i], ncol = 2)

      conf_mat_tmp[, 1] <- conf_level
      conf_mat_tmp[, 2] <- limits_tmp[, 1]
      conf_mat_tmp[, 3] <- limits_tmp[, 2]
      conf_mat_tmp[, 4] <- rep(i, length(conf_level))

      conf_mat <- rbind(conf_mat, conf_mat_tmp)

    }

    # Counternulls

    if (!is.null(null_values)) {

      counternull_mat_tmp <- matrix(NA, ncol = 3, nrow = length(null_values))

      counternull_mat_tmp[, 1] <- null_values
      counternull_mat_tmp[, 2] <- 2*estimate[i] - null_values
      counternull_mat_tmp[, 3] <- rep(i, length(null_values))

      counternull_mat <- rbind(counternull_mat, counternull_mat_tmp)

    }

  }

  res_frame <- as.data.frame(res_mat)

  names(res_frame) <- c("values", "conf_dist", "conf_dens", "p_two", "p_one", "variable")

  if (!is.null(conf_level)) {

    conf_frame <- as.data.frame(conf_mat)

    names(conf_frame) <- c("conf_level", "lwr", "upr", "variable")

  } else {
    conf_frame <- NULL
  }

  if (!is.null(null_values)) {

    counternull_frame <- as.data.frame(counternull_mat)

    names(counternull_frame) <- c("null_value", "counternull", "variable")

  } else {
    counternull_frame <- NULL
  }

  # Point estimators

  point_est_frame <- data.frame(
    est_mean = rep(NA, length(estimate))
    , est_median = rep(NA, length(estimate))
    , est_mode = rep(NA, length(estimate))
    , variable = seq(1, i, 1)
  )

  mean_fun <- function(x, estimate, stderr){
    x*dnorm((x - estimate)/stderr)*(1/stderr)
  }

  for (i in seq_along(estimate)) {

    point_est_frame$est_mean[i] <- integrate(mean_fun, lower = -Inf, upper = Inf, stderr = stderr[i], estimate = estimate[i], rel.tol = 1e-10)$value # Mean
    point_est_frame$est_median[i] <- res_frame$values[res_frame$variable %in% i][which.min(abs(res_frame$conf_dist[res_frame$variable %in% i][-1] - 0.5)) + 1] # Median
    point_est_frame$est_mode[i] <- res_frame$values[res_frame$variable %in% i][which.max(res_frame$conf_dens[res_frame$variable %in% i])] # Mode

  }

  return(list(
    res_frame = res_frame
    , conf_frame = conf_frame
    , counternull_frame = counternull_frame
    , point_est = point_est_frame
  ))

}


cdist_corr <- function(
  estimate = NULL
  , stderr = NULL
  , n = NULL
  , n_values = NULL
  , conf_level = NULL
  , null_values = NULL
  , alternative = NULL
){

  res_mat <- matrix(NA, nrow = 0, ncol = 6)

  conf_mat <- matrix(NA, nrow = 0, ncol = 4)

  counternull_mat <- matrix(NA, nrow = 0, ncol = 3)

  for (i in seq_along(estimate)) {

    limits <- c(-1, 1)

    x_calc <- c(estimate[i], null_values, seq(limits[1], limits[2], length.out = (n_values + 2)))
    x_calc <- x_calc[!(x_calc %in% c(-1, 1))]

    res_mat_tmp <- matrix(NA, nrow = length(x_calc), ncol = 6)

    res_mat_tmp[, 1] <- x_calc
    res_mat_tmp[, 2] <- 1 - pnorm((1/stderr[i])*(atanh(estimate[i]) - atanh(x_calc)))
    # res_mat_tmp[, 3] <- corr_dens(x = x_calc, estimate = estimate[i], stderr = stderr[i])
    res_mat_tmp[, 3] <- (-dnorm((1/stderr[i])*(atanh(estimate[i]) - atanh(x_calc)))*(-1/(stderr[i]*(1 - x_calc^2))))
    res_mat_tmp[, 4] <- 1 - 2*abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 5] <- (1/2) - abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 6] <- rep(i, length(x_calc))

    res_mat <- rbind(res_mat, res_mat_tmp)

    # Confidence intervals

    if (!is.null(conf_level)) {

      quants_tmp <- switch(
        alternative
        , two_sided = c(1 - (conf_level + 1)/2, (conf_level + 1)/2)
        , one_sided = c((1 - conf_level), conf_level)
      )

      # quants_tmp <- c(1 - (conf_level + 1)/2, (conf_level + 1)/2)

      conf_mat_tmp <- matrix(NA, ncol = 4, nrow = length(conf_level))

      limits_tmp <- matrix(tanh(qnorm(quants_tmp)*stderr[i] + atanh(estimate[i])), ncol = 2)

      conf_mat_tmp[, 1] <- conf_level
      conf_mat_tmp[, 2] <- limits_tmp[, 1]
      conf_mat_tmp[, 3] <- limits_tmp[, 2]
      conf_mat_tmp[, 4] <- rep(i, length(conf_level))

      conf_mat <- rbind(conf_mat, conf_mat_tmp)

    }

    # Counternulls

    if (!is.null(null_values)) {

      counternull_mat_tmp <- matrix(NA, ncol = 3, nrow = length(null_values))

      counternull_mat_tmp[, 1] <- null_values
      counternull_mat_tmp[, 2] <- tanh(2*atanh(estimate[i]) - atanh(null_values))
      counternull_mat_tmp[, 3] <- rep(i, length(null_values))

      counternull_mat <- rbind(counternull_mat, counternull_mat_tmp)

    }

  }

  res_frame <- as.data.frame(res_mat)

  names(res_frame) <- c("values", "conf_dist", "conf_dens", "p_two", "p_one", "variable")

  if (!is.null(conf_level)) {

    conf_frame <- as.data.frame(conf_mat)

    names(conf_frame) <- c("conf_level", "lwr", "upr", "variable")

  } else {
    conf_frame <- NULL
  }

  if (!is.null(null_values)) {

    counternull_frame <- as.data.frame(counternull_mat)

    names(counternull_frame) <- c("null_value", "counternull", "variable")

  } else {
    counternull_frame <- NULL
  }

  # Point estimators

  point_est_frame <- data.frame(
    est_mean = rep(NA, length(estimate))
    , est_median = rep(NA, length(estimate))
    , est_mode = rep(NA, length(estimate))
    , variable = seq(1, i, 1)
  )

  mean_fun <- function(x, stderr, estimate){
    x*((-dnorm((1/stderr)*(atanh(estimate) - atanh(x)))*(-1/(stderr*(1 - x^2)))))
  }

  for (i in seq_along(estimate)) {

    point_est_frame$est_mean[i] <- integrate(mean_fun, lower = -1, upper = 1, stderr = stderr[i], estimate = estimate[i], rel.tol = 1e-10)$value # Mean
    point_est_frame$est_median[i] <- res_frame$values[res_frame$variable %in% i][which.min(abs(res_frame$conf_dist[res_frame$variable %in% i][-1] - 0.5)) + 1] # Median
    point_est_frame$est_mode[i] <- res_frame$values[res_frame$variable %in% i][which.max(res_frame$conf_dens[res_frame$variable %in% i])] # Mode

  }

  return(list(
    res_frame = res_frame
    , conf_frame = conf_frame
    , counternull_frame = counternull_frame
    , point_est = point_est_frame
  ))

}

cdist_var <- function(
  estimate = NULL
  , n = NULL
  , n_values = NULL
  , conf_level = NULL
  , null_values = NULL
  , alternative = NULL
){

  eps <- 1e-10

  df <- (n - 1)

  res_mat <- matrix(NA, nrow = 0, ncol = 6)

  conf_mat <- matrix(NA, nrow = 0, ncol = 4)

  counternull_mat <- matrix(NA, nrow = 0, ncol = 3)

  for (i in seq_along(estimate)) {

    limits <- c(estimate[i]*df[i]/qchisq(eps, df = df[i], lower.tail = FALSE), estimate[i]*df[i]/qchisq(1 - eps, df = df[i], lower.tail = FALSE))

    x_calc <- c(estimate[i], null_values, seq(limits[1], limits[2], length.out = n_values))

    res_mat_tmp <- matrix(NA, nrow = length(x_calc), ncol = 6)

    res_mat_tmp[, 1] <- x_calc
    res_mat_tmp[, 2] <- (1 - pchisq(df[i]*estimate[i]/x_calc, df = df[i]))
    # res_mat_tmp[, 3] <- var_dens(x = x_calc, estimate = estimate[i], df = df[i])
    res_mat_tmp[, 3] <- (-dchisq(df[i]*estimate[i]/x_calc, df = df[i])*(-((df[i]*estimate[i])/x_calc^2)))
    res_mat_tmp[, 4] <- 1 - 2*abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 5] <- (1/2) - abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 6] <- rep(i, length(x_calc))

    res_mat <- rbind(res_mat, res_mat_tmp)

    # Confidence intervals

    if (!is.null(conf_level)) {

      quants_tmp <- switch(
        alternative
        , two_sided = c(1 - (conf_level + 1)/2, (conf_level + 1)/2)
        , one_sided = c((1 - conf_level), conf_level)
      )

      # quants_tmp <- c(1 - (conf_level + 1)/2, (conf_level + 1)/2)

      conf_mat_tmp <- matrix(NA, ncol = 4, nrow = length(conf_level))

      limits_tmp <- matrix((estimate[i]*df[i])/qchisq(quants_tmp, df = df[i], lower.tail = FALSE), ncol = 2)

      conf_mat_tmp[, 1] <- conf_level
      conf_mat_tmp[, 2] <- limits_tmp[, 1]
      conf_mat_tmp[, 3] <- limits_tmp[, 2]
      conf_mat_tmp[, 4] <- rep(i, length(conf_level))

      conf_mat <- rbind(conf_mat, conf_mat_tmp)

    }

    # Counternulls

    if (!is.null(null_values)) {

      counter_tmp <- 1 - (1 - pchisq(df[i]*estimate[i]/null_values, df = df[i]))

      counternull_mat_tmp <- matrix(NA, ncol = 3, nrow = length(null_values))

      counternull_mat_tmp[, 1] <- null_values
      counternull_mat_tmp[, 2] <- estimate[i]*df[i]/qchisq(counter_tmp, df = df[i], lower.tail = FALSE)
      counternull_mat_tmp[, 3] <- rep(i, length(null_values))

      counternull_mat <- rbind(counternull_mat, counternull_mat_tmp)

    }

  }

  res_frame <- as.data.frame(res_mat)

  names(res_frame) <- c("values", "conf_dist", "conf_dens", "p_two", "p_one", "variable")

  if (!is.null(conf_level)) {

    conf_frame <- as.data.frame(conf_mat)

    names(conf_frame) <- c("conf_level", "lwr", "upr", "variable")

  } else {
    conf_frame <- NULL
  }

  if (!is.null(null_values)) {

    counternull_frame <- as.data.frame(counternull_mat)

    names(counternull_frame) <- c("null_value", "counternull", "variable")

  } else {
    counternull_frame <- NULL
  }

  # Point estimators

  point_est_frame <- data.frame(
    est_mean = rep(NA, length(estimate))
    , est_median = rep(NA, length(estimate))
    , est_mode = rep(NA, length(estimate))
    , variable = seq(1, i, 1)
  )

  point_est_frame$est_mean <- exp((log(estimate) + log(df) + lgamma(-1 + (df/2))) - (log(2) + lgamma(df/2)))
  point_est_frame$est_median <- estimate*df/(2*zipfR::Rgamma.inv(df/2, 1/2))
  point_est_frame$est_mode <- exp((log(estimate) + log(df)) - (log(2 + df)))


  return(list(
    res_frame = res_frame
    , conf_frame = conf_frame
    , counternull_frame = counternull_frame
    , point_est = point_est_frame
  ))

}


wilson_ci <- function(
  estimate
  , n
  , conf_level
  , alternative
) {

  # if (alternative %in% "two_sided") {
  #   z <- qnorm((conf_level + 1)/2)
  # } else if (alternative %in% "one_sided") {
  #   z <- qnorm(conf_level)
  # }

  z <- qnorm((conf_level + 1)/2)

  p1 <- estimate + (1/2)*z^2/n
  p2 <- z*sqrt((estimate*(1 - estimate) + (1/4)*z^2/n)/n)
  p3 <- 1 + z^2/n

  c((p1 - p2)/p3, (p1 + p2)/p3)

}

wilson_cicc <- function(
  estimate
  , n
  , conf_level
  , alternative
) {

  z <- qnorm((conf_level + 1)/2)
  x <- round(estimate*n) # To get number of successes/failures
  estimate_compl <- (1 - estimate) # Complement of estimate

  lower <- max(0, (2*x + z^2 - 1 - z*sqrt(z^2 - 2 - 1/n + 4*estimate*(n*estimate_compl + 1)))/(2*(n + z^2)), na.rm = TRUE)
  upper <- min(1, (2*x + z^2 + 1 + z*sqrt(z^2 + 2 - 1/n + 4*estimate*(n*estimate_compl - 1)))/(2*(n + z^2)), na.rm = TRUE)

  c(lower, upper)

}

wilson_cicc_diff <- function(
  estimate
  , n
  , conf_level
  , alternative
) {

  est_diff <- (estimate[1] - estimate[2])

  res1 <- wilson_cicc(
    estimate = estimate[1]
    , n = n[1]
    , conf_level = conf_level
    , alternative = alternative
  )

  res2 <- wilson_cicc(
    estimate = estimate[2]
    , n = n[2]
    , conf_level = conf_level
    , alternative = alternative
  )

  l1 <- res1[1]
  u1 <- res1[2]
  l2 <- res2[1]
  u2 <- res2[2]

  lim1 <- max(-1, est_diff - sqrt((estimate[1] - l1)^2 + (u2 - estimate[2])^2))
  lim2 <- min(1, est_diff + sqrt((u1 - estimate[1])^2 + (estimate[2] - l2)^2))

  sort(c(lim1, lim2), decreasing = FALSE)

}

cdist_prop1 <- function(
  estimate = NULL
  , n = NULL
  , n_values = NULL
  , conf_level = NULL
  , null_values = NULL
  , alternative = NULL
){

  # Auxilliary functions

  cdf_fun <- function(x, n, p) {
    x <- as.complex(x)
    n <- as.complex(n)
    p <- as.complex(p)

    -Re((1i*sqrt(n)*(p - x))/(sqrt(x - 1)*sqrt(x)))
  }

  deriv_fun <- function(x, n, p) {
    x <- as.complex(x)
    n <- as.complex(n)
    p <- as.complex(p)

    Re((1i*sqrt(n)*(-x + p*(-1 + 2*x)))/(2*(-1 + x)^(3/2)*x^(3/2)))
  }

  eps <- 1e-10

  res_mat <- matrix(NA, nrow = 0, ncol = 6)

  conf_mat <- matrix(NA, nrow = 0, ncol = 4)

  counternull_mat <- matrix(NA, nrow = 0, ncol = 3)

  for (i in seq_along(estimate)) {

    limits <- wilson_ci(estimate = estimate[i], n = n[i], conf_level = (1 - eps), alternative = alternative)

    x_calc <- c(estimate[i], null_values, seq(limits[1], limits[2], length.out = n_values))

    res_mat_tmp <- matrix(NA, nrow = length(x_calc), ncol = 6)

    res_mat_tmp[, 1] <- x_calc
    res_mat_tmp[, 2] <- pnorm(cdf_fun(x = x_calc, n = n[i], p = estimate[i]))
    res_mat_tmp[, 3] <- dnorm(cdf_fun(x = x_calc, n = n[i], p = estimate[i]))*deriv_fun(x = x_calc, n = n[i], p = estimate[i])
    res_mat_tmp[, 4] <- 1 - 2*abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 5] <- (1/2) - abs(res_mat_tmp[, 2] - (1/2))
    res_mat_tmp[, 6] <- rep(i, length(x_calc))

    res_mat <- rbind(res_mat, res_mat_tmp)

    # Confidence intervals

    if (!is.null(conf_level)) {

      quants_tmp <- switch(
        alternative
        , two_sided = c(1 - (conf_level + 1)/2, (conf_level + 1)/2)
        , one_sided = c((1 - conf_level), conf_level)
      )

      # quants_tmp <- c(1 - (conf_level + 1)/2, (conf_level + 1)/2)

      conf_mat_tmp <- matrix(NA, ncol = 4, nrow = length(conf_level))

      limits_tmp <- matrix(wilson_ci(estimate = estimate[i], n = n[i], conf_level = conf_level, alternative = alternative), ncol = 2)

      conf_mat_tmp[, 1] <- conf_level
      conf_mat_tmp[, 2] <- limits_tmp[, 1]
      conf_mat_tmp[, 3] <- limits_tmp[, 2]
      conf_mat_tmp[, 4] <- rep(i, length(conf_level))

      conf_mat <- rbind(conf_mat, conf_mat_tmp)

    }

    # Counternulls

    if (!is.null(null_values)) {

      counter_tmp <- 1 - (pnorm(cdf_fun(x = null_values, n = n[i], p = estimate[i])))

      counternull_mat_tmp <- matrix(NA, ncol = 3, nrow = length(null_values))

      counternull_mat_tmp[, 1] <- null_values
      counternull_mat_tmp[, 2] <- sapply(counter_tmp, function(x, cdf, q){x[which.min(abs(cdf - q))]}, x = x_calc, cdf = res_mat_tmp[, 2])
      counternull_mat_tmp[, 3] <- rep(i, length(null_values))

      counternull_mat <- rbind(counternull_mat, counternull_mat_tmp)

    }

  }

  res_frame <- as.data.frame(res_mat)

  names(res_frame) <- c("values", "conf_dist", "conf_dens", "p_two", "p_one", "variable")

  if (!is.null(conf_level)) {

    conf_frame <- as.data.frame(conf_mat)

    names(conf_frame) <- c("conf_level", "lwr", "upr", "variable")

  } else {
    conf_frame <- NULL
  }

  if (!is.null(null_values)) {

    counternull_frame <- as.data.frame(counternull_mat)

    names(counternull_frame) <- c("null_value", "counternull", "variable")

  } else {
    counternull_frame <- NULL
  }

  # Point estimators

  point_est_frame <- data.frame(
    est_mean = rep(NA, length(estimate))
    , est_median = rep(NA, length(estimate))
    , est_mode = rep(NA, length(estimate))
    , variable = seq(1, i, 1)
  )

  mean_fun <- function(x, estimate, n){
    x*dnorm(cdf_fun(x = x, n = n, p = estimate))*deriv_fun(x = x, n = n, p = estimate)
  }

  for (i in seq_along(estimate)) {

    point_est_frame$est_mean[i] <- integrate(mean_fun, lower = 0, upper = 1, n = n[i], estimate = estimate[i], rel.tol = 1e-10)$value # Mean
    point_est_frame$est_median[i] <- res_frame$values[res_frame$variable %in% i][which.min(abs(res_frame$conf_dist[res_frame$variable %in% i][-1] - 0.5)) + 1] # Median
    point_est_frame$est_mode[i] <- res_frame$values[res_frame$variable %in% i][which.max(res_frame$conf_dens[res_frame$variable %in% i])] # Mode

  }

  return(list(
    res_frame = res_frame
    , conf_frame = conf_frame
    , counternull_frame = counternull_frame
    , point_est = point_est_frame
  ))

}

cdist_propdiff <- function(
  estimate = NULL
  , n = NULL
  , n_values = NULL
  , conf_level = NULL
  , null_values = NULL
  , alternative = NULL
){

  eps <- 1e-15

  res_mat <- matrix(NA, nrow = 0, ncol = 6)

  conf_mat <- matrix(NA, nrow = 0, ncol = 4)

  counternull_mat <- matrix(NA, nrow = 0, ncol = 3)

  # for (i in seq_along(estimate)) {

  # limits <- wilson_cicc_diff(estimate = estimate, n = n, conf_level = (1 - eps), alternative = alternative)

  conf_levels <- seq(eps, 1 - eps, length.out = ceiling(n_values/2))
  x_calc <- sapply(conf_levels, wilson_cicc_diff, estimate = estimate, n = n, alternative = alternative)

  val_min <- wilson_cicc_diff(estimate, n, conf_level = eps)
  val_between <- seq(min(val_min), max(val_min), length.out = 100)

  res_mat_tmp <- matrix(NA, nrow = length(x_calc) + 100, ncol = 6)

  res_mat_tmp[, 1] <- c(x_calc[1, ], x_calc[2, ], val_between)
  is.na(res_mat_tmp[, 2]) <- TRUE # No confidence distribution for this one
  is.na(res_mat_tmp[, 3]) <- TRUE # No confidence density for this one
  res_mat_tmp[, 4] <- c(1 - conf_levels, 1 - conf_levels, rep(NA, 100))
  res_mat_tmp[, 5] <- c((1 - conf_levels)/2, (1 - conf_levels)/2, rep(NA, 100))
  res_mat_tmp[, 6] <- rep(1, times = (length(x_calc) + 100))

  res_mat <- rbind(res_mat, res_mat_tmp)

  # Confidence intervals

  if (!is.null(conf_level)) {

    conf_tmp <- switch(
      alternative
      , two_sided = conf_level
      , one_sided = 2*conf_level - 1
    )

    conf_mat_tmp <- matrix(NA, ncol = 4, nrow = length(conf_level))

    limits_tmp <- matrix(NA, ncol = 2, nrow = length(conf_tmp))

    for (j in seq_along(conf_tmp)) {
      limits_tmp[j, ] <- wilson_cicc_diff(estimate = estimate, n = n, conf_level = conf_tmp[j], alternative = alternative)
    }

    conf_mat_tmp[, 1] <- conf_level
    conf_mat_tmp[, 2] <- limits_tmp[, 1]
    conf_mat_tmp[, 3] <- limits_tmp[, 2]
    conf_mat_tmp[, 4] <- rep(1, length(conf_level))

    conf_mat <- rbind(conf_mat, conf_mat_tmp)

  }

  # Counternulls

  if (!is.null(null_values)) {

    cnull_tmp <- rep(NA, length(null_values))

    tmp_fun_up <- function(conf_level, estimate, n, null_values) {
      wilson_cicc_diff(estimate = estimate, n = n, conf_level = conf_level)[2] - null_values
    }

    tmp_fun_low <- function(conf_level, estimate, n, null_values) {
      wilson_cicc_diff(estimate = estimate, n = n, conf_level = conf_level)[1] - null_values
    }

    for (j in seq_along(null_values)) {

      if ((null_values[j] > max(res_mat_tmp[, 1])) || (null_values[j] < min(res_mat_tmp[, 1]))) {
        is.na(cnull_tmp[j]) <- TRUE # Set values in the "gap" to missing for plotting
      } else {

        if (null_values[j] > -diff(estimate)) {

          null_conf <- uniroot(
            tmp_fun_up
            , lower = 1e-15
            , upper = 1 - 1e-15
            , null_values = null_values[j]
            , estimate = estimate
            , n = n
          )$root

          cnull_tmp[j] <- wilson_cicc_diff(estimate = estimate, n = n, conf_level = null_conf)[1]

        } else if (null_values[j] < -diff(estimate)) {

          null_conf <- uniroot(
            tmp_fun_low
            , lower = 1e-15
            , upper = 1 - 1e-15
            , null_values = null_values[j]
            , estimate = estimate
            , n = n
          )$root

          cnull_tmp[j] <- wilson_cicc_diff(estimate = estimate, n = n, conf_level = null_conf)[2]
        }
      }
    }


    counternull_mat_tmp <- matrix(NA, ncol = 3, nrow = length(null_values))
    counternull_mat_tmp[, 1] <- null_values
    counternull_mat_tmp[, 2] <- cnull_tmp
    counternull_mat_tmp[, 3] <- rep(1, length(null_values))
    counternull_mat <- rbind(counternull_mat, counternull_mat_tmp)

  }

  # }

  res_frame <- as.data.frame(res_mat)

  names(res_frame) <- c("values", "conf_dist", "conf_dens", "p_two", "p_one", "variable")

  if (!is.null(conf_level)) {

    conf_frame <- as.data.frame(conf_mat)

    names(conf_frame) <- c("conf_level", "lwr", "upr", "variable")

  } else {
    conf_frame <- NULL
  }

  if (!is.null(null_values)) {

    counternull_frame <- as.data.frame(counternull_mat)

    names(counternull_frame) <- c("null_value", "counternull", "variable")

  } else {
    counternull_frame <- NULL
  }

  # Point estimators

  point_est_frame <- data.frame(
    est_mean = rep(NA, 1)
    , est_median = rep(NA, 1)
    , est_mode = rep(NA, 1)
    , variable = seq(1, 1, 1)
  )

  return(list(
    res_frame = res_frame
    , conf_frame = conf_frame
    , counternull_frame = counternull_frame
    , point_est = point_est_frame
  ))

}

# Define new mixed scale: log for p <= 0.05, else linear

magnify_trans_log <- function(interval_low = 0.05, interval_high = 1,  reducer = 0.05, reducer2 = 8) {

  trans <- Vectorize(function(x, i_low = interval_low, i_high = interval_high, r = reducer, r2 = reducer2) {
    if(is.na(x) || (x >= i_low & x <= i_high)) {
      x
    } else if(x < i_low & !is.na(x)) {
      (log10(x / r)/r2 + i_low)
    } else {
      log10((x - i_high) / r + i_high)/r2
    }
  })

  inv <- Vectorize(function(x, i_low = interval_low, i_high = interval_high, r = reducer, r2 = reducer2) {
    if(is.na(x) || (x >= i_low & x <= i_high)) {
      x
    } else if(x < i_low & !is.na(x)) {
      10^(-(i_low - x)*r2)*r
    } else {
      i_high + 10^(x*r2)*r - i_high*r
    }
  })

  trans_new(name = 'customlog', transform = trans, inverse = inv, domain = c(1e-16, Inf))
}

magnify_trans_log_rev <- function(interval_low = 0.05, interval_high = 1,  reducer = 0.05, reducer2 = 8) {

  trans <- Vectorize(function(x, i_low = interval_low, i_high = interval_high, r = reducer, r2 = reducer2) {
    -if(is.na(x) || (x >= i_low & x <= i_high)) {
      x
    } else if(x < i_low & !is.na(x)) {
      (log10(x / r)/r2 + i_low)
    } else {
      log10((x - i_high)/r + i_high)/r2 + i_high
    }
  })

  inv <- Vectorize(function(x, i_low = interval_low, i_high = interval_high, r = reducer, r2 = reducer2) {
    if(is.na(x) || (-x >= i_low & -x <= i_high)) {
      -x
    } else if(-x < i_low & !is.na(x)) {
      (10^(-(i_low + x)*r2)*r)
    } else {
      i_high + 10^(-r2*(i_high + x))*r - i_high*r
    }
  })

  trans_new(name = 'customlog_rev', transform = trans, inverse = inv, domain = c(1e-16, Inf))
}

# # Surprisal transformation
#
# trans_surprisal <- function() {
#
#   trans <- Vectorize(function(x) {
#     if(!is.na(x)) {
#       -log2(x)
#     }
#   })
#
#   inv <- Vectorize(function(x) {
#     if(!is.na(x)) {
#       2^(-x)
#     }
#   })
#
#   trans_new(name = 'surprisal', transform = trans, inverse = inv, domain = c(1e-16, Inf))
# }
