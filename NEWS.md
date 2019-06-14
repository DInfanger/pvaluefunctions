## NEWS for the pvaluefunctions package

1.3.0
-------------

  * NEWS file was converted to .md file.
  * 

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
