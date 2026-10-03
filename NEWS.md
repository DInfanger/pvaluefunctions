## NEWS and changes for the pvaluefunctions package

Development version
-------------

Bug fixes

  * Point estimates for types based on the t- and normal distribution: the median is now exactly the estimate. Previously it was read off a coarse grid and could be far off, especially for small degrees of freedom. The mean is `NA` for `df <= 1` because it does not exist.
  * The shaded logarithmic part of the y-axis (`log_yaxis = TRUE`) now spans the whole logarithmic x-axis. Previously it stopped at x = 100.
  * The warnings for the approximations of Spearman's and Kendall's correlation now also apply to negative estimates.
  * Variances: the mean of the confidence distribution is now `NA` for `n <= 3` because it does not exist. Previously it was `Inf` for `n = 3` and a meaningless finite value for `n = 2`.
  * Differences of proportions: null values close to the estimated difference (inside the gap left by the continuity correction) no longer fail with "f() values at end points not of opposite sign". Their counternull is `NA`.
  * A logarithmic y-axis (`log_yaxis = TRUE`) no longer fails with "wrong sign in 'by' argument" if `plot_p_limit` is close to `cut_logyaxis` (e.g. one-sided with `plot_p_limit = 0.02`).
  * Pearson's correlation coefficient: the exact confidence distribution no longer exceeds 1 because of integration error close to r = 1. This gave negative *p*-values, `NaN` s-values and the warning "NaNs produced".

Changes

  * `conf_dist()` now returns its results invisibly. The plot is still printed if `plot = TRUE`, but calling `conf_dist()` without assigning the result no longer prints all data frames to the console.
  * With `plot = FALSE`, the plot is no longer built in the background. As a consequence, the message about null values outside of the x-axis range only appears if a plot is created.

Documentation

  * The help page of `conf_dist()` now lists all returned elements. `aucc_frame` is always returned, not only if `plot = TRUE`.
  * The dependencies listed in the README and the vignette are up to date.

Internal changes

  * The plotting code of `conf_dist()` is now in a separate internal function.
  * `devtools` was removed from the suggested packages.
  * Generated vignette files are no longer tracked in git and are excluded from the package build.

1.7.0
-------------

Bug fixes

  * `x_scale = "logarithm"` is now honored. Previously only the automatic logarithmic scale for `trans = "exp"` worked.
  * Pearson's correlation coefficient: the exact distribution no longer fails with `n` larger than about 172 (overflow of the gamma function), and the counternull is now computed accurately for extreme tail probabilities. Several Pearson estimates in one call (`n` and `estimate` of length > 1) no longer fail.
  * The lengths of `stderr`, `df`, `tstat` and `estimate` are now checked correctly. Previously a mismatch of the third argument went unnoticed.
  * `trans` can now be a function as documented, not only the name of a function. Names are no longer forced to lower case when looking up the function (a lower case name is still tried as a fallback), functions defined inside the calling function are found, and the transformation function does not need an argument named `x`.
  * An invalid `type` now gives an informative error message.
  * One-sided confidence limits for proportions (`type = "prop"`) now use the matching two-sided level (`2 * conf_level - 1`) like all other types. Previously, the two-sided interval of level `conf_level` was returned.
  * The labels for the confidence levels no longer trigger a deprecation warning with ggplot2 4.0.0 or later (`label.size`).
  * Null values that are outside of the plotting area or whose counternull is missing no longer cause errors when plotting.

Input validation

  * Inputs are now validated up front and give informative errors, e.g. for missing values, non-positive standard errors, degrees of freedom or sample sizes, non-logical flags, an invalid `n_values`, `plot_p_limit`, `cut_logyaxis`, `nrow` or `ncol`, and duplicated `est_names`.
  * Correlation coefficients of exactly -1 or 1 are rejected. The minimum sample size is 4 for Pearson's correlation (n = 3 never worked).
  * t-statistics must be non-zero and have the same sign as the estimates.
  * A message is printed if confidence levels are ignored because they are too low for the alternative.

Internal changes

  * Code is split into several files, duplicated code was removed and imports are now explicit. `ggplot2 (>= 3.5.0)` and `scales (>= 1.3.0)` are now required (`transform` instead of `trans`, `new_transform()`).
  * `README.Rmd` and tooling files are now correctly excluded from the package build.
  * Added a testthat test suite (closed-form checks against `stats`, input validation, plot properties and regression tests against the results of version 1.6.3).

1.6.3
-------------

  * Deprecated arguments in ggplot2 were changed, such as "size" in "geom_line". No more messages should appear now.
  * Small changes in the references in the vignette.

1.6.2
-------------

  * The calculations for Pearson's correlation coefficient now use the exact distribution described in a preprint by Gunnar Taraldsen (2020). This makes the use of the gsl package necessary because the calculations involve the Gaussian hypergeometric function (2F1). Warning: This can make the calculations drastically longer if a high value of `n_values` is used.
  * The function `conf_dist` now also calculates the proportion of the AUCC that lies above any specified null values.
  * Small changes in the help page for `conf_dist`.

1.6.1
-------------

  * Fixed a bug concerning the calculation of Newcombe's Wilson score interval with continuity correction for the difference of two proportions.
  * Removed the "cairo" device from pngs in the vignette.

1.6.0
-------------

  * The returned data frame is now sorted for convenience; this is purely cosmetic.
  * Dependence on R increased to R version 3.5.0
  * Added an option `same_color` to specify whether curves should be distinguished by colors or not if they are plotted together in the same graph. Can be useful if there are many curves plotted together.
  * Added an option `plot_legend` to specify whether a legend should be drawn if multiple curves are plotted together and distinguished by color (i.e. `same_color = FALSE` and `together = TRUE`).
  * Added an option `col` to specify the color of the curves if they are not to be distinguished by color (i.e. `same_color = FALSE`).
  * Areas under the confidence curves (AUCC) according to Berrar (2017) are now calculated and returned. They offer a way to compare multiple estimates with respect to their precision. The AUCC is calculated on the untransformed scale using numerical integration (trapezoidal integration) implemented in the `pracma` package which is now imported.
  * Removed data link to external data (UCLA) in an example (odds ratio). Certificates for this site apparently expired.


1.5.0
-------------

  * Vignette, README and DESCRIPTION updated to reflect that the newest version of ggplot2 (3.2.1) fixes the former bug with `sec_axis`.
  * Added an option `plot` to `conf_dist` that controls whether a plot is created or not. If users want to create their own plots, they can set this option to `FALSE` and use the returned data (`res_frame`) which is the basis for the plots to create them.
  

1.4.0
-------------
  
  * Added option `inverted` which allows users to plot p-value functions, s-value functions and confidence distributions with the y-axis inverted.
  * Added a new example showing a *p*-value function for an odds ratio with an inverted y-axis (cf. Bender et al. 2005).
  * The option `xlim` is now strictly enforced: Any null values that are outside of the specified x-axis-limits are not plotted and a corresponding message is printed out as information.
  * Added new option `x_scale` to manually force the scaling of the x-axis.
  * Added two more examples in the vignette replicating Figure 1 and Figure 2 from Bender et al. (2005).
  * Various smaller bug fixes and improvements: Fixed some checks, fixed plotting vertical lines for null values.

1.3.0
-------------

  * NEWS file was converted to an .md file.
  * Users can now provide a title for the plot (option `title`).
  * Users can now provide titles for the primary and secondary y-axis (options `ylab` and `ylab_sec`).
  * Users can now specify the number of rows `nrow` and columns `ncol` to be used in `facet_wrap` (ggplot2) when multiple estimates are plotted separately (option `together = FALSE`).
  * Changed some of the examples: Changed option `log_yaxis = TRUE` to `log_yaxis = FALSE`.
  * Code: Changed logical checks `x == TRUE` and `x == FALSE` to `isTRUE(x)` and to `isFALSE(x)`.
  * Code: Checked and improved/fixed some of the initial consistency/input checks that are performed at the beginning of the function.

1.2.0
-------------

  * Added option `plot_counternull` to plot the counternull value(s) on the graphics if applicable and possible.
  * Effect sizes on the log-scale (e.g. Odds ratio, Hazard ratio, Incidence rate ratio) are now plotted on a logarithmic x-axis so that the p-value function is symmetric around the point estimate.
  * The default is now to omit a logarithmic part of the y-axis.
  * Warnings from ggplot2 are now suppressed during the function call.
  * Not interesting for practitioners: Source code for creating the plots was cleaned up and more comments were added.

1.1.0
-------------

  * Improved vignette with reduced figure size and several orthographic errors fixed. Added a new example (difference between proportions)
  * New estimate type implemented: Difference between two independent proportions (`type = "propdiff"`) based on Wilson's interval (see Newcombe (1998) for details)

1.0.0
-------------

  * Initial release
