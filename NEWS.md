## NEWS and changes for the pvaluefunctions package

1.4.0
-------------
  
  * Added option `inverted` which allows users to plot p-value functions, s-value functions and confidence distributions with the y-axis inverted.
  * Added a new example.

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
