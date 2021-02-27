# Rphylopars
Rphylopars is an R package for conducting multivariate phylogenetic comparative analyses on datasets with missing observations and missing data. Rphylopars uses a fast linear-time algorithm and incorporate a variety of evolutionary models, including estimation of tree transformation parameters (Early-Burst, Ornstein-Uhlenbeck, lambda, kappa, delta) as well as the multivariate Ornstein-Uhlenbeck model.

For download information and tutorials, visit the [Rphylopars wiki](https://github.com/ericgoolsby/Rphylopars/wiki).

## Current versions available
### CRAN: 0.3.2
### GitHub master branch - 0.3.2
### GitHub devel branch - 0.3.6

## Version notes
* Vesrion 0.3.6 - New experimental internal function `Rphylopars:::get_cov_CIs` to estimate 95% CIs for trait covariance parameters from a fitted `phylopars` object.
* Vesrion 0.3.5 - Checks that all species names in `trait_data$species` have an exact match in `tree$tip.label`
* Version 0.3.4 - Fixes printing error when only one variable is included in `phylopars()`
* Vesion 0.3.3 - Prevents R from crashing when is.na(species) (e.g. trailing blank rows at the end of a spreadsheet).
* Version 0.3.1 and 0.3.2 - Removed broken links from fast.SSC.Rd and Rphylopars-package.Rd help files (resolves CRAN notes)
* Version 0.3.0 - R 4.0 compatibility
* Version 0.2.11 - Fixed confidence intervals (previously too large) and installation issues
