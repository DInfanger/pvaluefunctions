## NEWS and changes for the pvaluefunctions package

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
