pvaluefunctions
================

  - [*P*-value functions
    <img src="man/figures/logo3.svg" align="right" width="120" />](#p-value-functions)
      - [Accompanying paper](#accompanying-paper)
          - [Recreation of the figures in the
            paper](#recreation-of-the-figures-in-the-paper)
      - [Overview](#overview)
      - [Installation](#installation)
      - [Dependencies](#dependencies)
      - [Examples](#examples)
      - [References](#references)
      - [Contact](#contact)
      - [Session info](#session-info)
      - [License](#license)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# *P*-value functions <img src="man/figures/logo3.svg" align="right" width="120" />

[![CRAN\_Status\_Badge](https://cranlogs.r-pkg.org:443/badges/grand-total/pvaluefunctions)](https://cran.r-project.org/package=pvaluefunctions)
[![downloads](https://cranlogs.r-pkg.org/badges/pvaluefunctions)](https://cran.r-project.org/package=pvaluefunctions)
[![total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/pvaluefunctions)](http://cranlogs.r-pkg.org/badges/grand-total/pvaluefunctions)
[![Rdoc](https://www.rdocumentation.org/packages/pvaluefunctions)](http://www.rdocumentation.org/packages/pvaluefunctions)

## Accompanying paper

We published an [accompanying paper](https://doi.org/10.1002/sim.8293)
to illustrate the use of *p*-value functions:

Infanger D, Schmidt-Trucksäss A. (2019): *P* value functions: An
underused method to present research results and to promote quantitative
reasoning. *Statistics in Medicine.* **38**: 4189-4197. doi:
10.1002/sim.8293.

### Recreation of the figures in the paper

The code and instructions to reproduce all graphics in our paper can be
found in the following GitHub repository:
<https://github.com/DInfanger/pvalue_functions>

## Overview

This is the repository for the R-package
[`pvaluefunctions`](https://cran.r-project.org/package=pvaluefunctions).
The package contains R functions to create graphics of *p*-value
functions, confidence distributions, confidence densities, or the
[Surprisal value
(S-value)](http://www.umsl.edu/~fraundorfp/egsurpri.html) (Greenland
2019).

## Installation

You can install the package directly from CRAN by typing
`install.packages("pvaluefunctions")`. After installation, load it in R
using `library(pvaluefunctions)`.

## Dependencies

The function depends on the following R packages, which need to be
installed beforehand:

  - [ggplot2](https://cran.r-project.org/package=ggplot2)
  - [scales](https://cran.r-project.org/package=scales)
  - [zipfR](https://cran.r-project.org/package=zipfR)
  - [pracma](https://cran.r-project.org/package=pracma)

Use the command `install.packages(c("ggplot2", "scales", "zipfR",
"pracma"))` in R to install those packages.

## Examples

For more examples and code, see the
[vignette](https://CRAN.R-project.org/package=pvaluefunctions/vignettes/pvaluefun.html).

<img src="man/figures/README-ttest_pval-1.png" width="70%" style="display: block; margin: auto auto auto 0;" />

<img src="man/figures/README-ttest_sval-1.png" width="70%" style="display: block; margin: auto auto auto 0;" />

<img src="man/figures/README-benderfig1-1.png" width="70%" style="display: block; margin: auto auto auto 0;" />

## References

Bender R, Berg G, Zeeb H. (2005): Tutorial: using confidence curves in
medical research. *Biom J.* 47(2): 237-47.

Berrar D (2017): Confidence Curves: an alternative to null hypothesis
significance testing for the comparison of classifiers. *Mach Learn.*
106:911-949.

Fraser D. A. S. (2019): The *p*-value function and statistical
inference. *Am Stat.* 73:sup1, 135-147.

Greenland S (2019): Valid *P*-Values Behave Exactly as They Should: Some
Misleading Criticisms of *P*-Values and Their Resolution with
*S*-Values. *Am Stat.* 73sup1, 106-114.

Infanger D, Schmidt-Trucksäss A. (2019): *P* value functions: An
underused method to present research results and to promote quantitative
reasoning. *Stat Med.* 38, 4189-4197. doi: 10.1002/sim.8293.

Poole C. (1987a): Beyond the confidence interval. *Am J Public Health.*
77(2): 195-9.

Poole C. (1987b): Confidence intervals exclude nothing. *Am J Public
Health.* 77(4): 492-3.

Rafi Z, Greenland S. (2020): Semantic and cognitive tools to aid
statistical science: replace confidence and significance by
compatibility and surprise. *BMC Med Res Methodol.* 20, 244. doi:
10.1186/s12874-020-01105-9.

Rosenthal R, Rubin DB. (1994): The counternull value of an effect size:
A new statistic. *Psychol Sci.* 5(6): 329-34.

Schweder T, Hjort NL. (2016): Confidence, likelihood, probability:
statistical inference with confidence distributions. New York, NY:
Cambridge University Press.

Xie M, Singh K, Strawderman WE. (2011): Confidence Distributions and a
Unifying Framework for Meta-Analysis. *J Am Stat Assoc.* 106(493):
320-33. doi: 10.1198/jasa.2011.tm09803.

Xie Mg, Singh K. (2013): Confidence distribution, the frequentist
distribution estimator of a parameter: A review. *Internat Statist Rev.*
81(1): 3-39.

## Contact

[Denis Infanger](https://dsbg.unibas.ch/de/personen/denis-infanger/)

## Session info

    #> R Under development (unstable) (2020-12-07 r79587)
    #> Platform: x86_64-w64-mingw32/x64 (64-bit)
    #> Running under: Windows 10 x64 (build 19042)
    #> 
    #> Matrix products: default
    #> 
    #> locale:
    #> [1] LC_COLLATE=German_Switzerland.1252  LC_CTYPE=German_Switzerland.1252   
    #> [3] LC_MONETARY=German_Switzerland.1252 LC_NUMERIC=C                       
    #> [5] LC_TIME=German_Switzerland.1252    
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] pvaluefunctions_1.6.1
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] knitr_1.30       magrittr_2.0.1   munsell_0.5.0    colorspace_2.0-0
    #>  [5] R6_2.5.0         rlang_0.4.9      stringr_1.4.0    tools_4.1.0     
    #>  [9] grid_4.1.0       gtable_0.3.0     xfun_0.19        htmltools_0.5.0 
    #> [13] ellipsis_0.3.1   yaml_2.2.1       digest_0.6.27    tibble_3.0.4    
    #> [17] lifecycle_0.2.0  crayon_1.3.4     farver_2.0.3     ggplot2_3.3.2   
    #> [21] vctrs_0.3.5      glue_1.4.2       evaluate_0.14    rmarkdown_2.5   
    #> [25] pracma_2.2.9     stringi_1.5.3    compiler_4.1.0   pillar_1.4.7    
    #> [29] scales_1.1.1     pkgconfig_2.0.3

## License

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
