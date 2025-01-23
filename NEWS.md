Version 2.0.1
=============

* Added VA implementation of Tweedie
* Added "band matrix" as a sparsity pattern for phylogenetic models, as alternative to the NNGP
* Added 'tick.length' argument for phyplot.gllvm to control the tick length in the species-specific random effects plot

## Bug Fixes
* Fixed a bug in the calculation of standard errors for models involving traits and a ZIP/ZINB response distribution. See #204.
* Fixed a bug that prevented from incorporating correlated random effects via "formula"
* Fixed a bug that prevented successfully incorporating phylogenetic random effects with traits
* Fixed a bug in the confidence level of phyloplot.gllvm

Version 2.0
=============

* For CRAN release 2.0 see updates for versions 1.4.4 - 1.4.9

Version 1.4.9
=============

* Row.eff can now be used for community-level (species-common) effect
* Both fixed and random at the same time (i.e., a mixed effects formula)
  * Does not allow for a single random intercept
  * Does not yet allow for between random effect correlation
* New formula interface for phylogenetic model adapted to trait model too
* New phyplot.gllvm function for plotting the phylogenetic random effects
* Minor adjustment in the behavior of 'caption' in plot.gllvm

Version 1.4.8
=============

* Added functionality for correlated random canonical coefficients
* Changed "site.index" argument in getResidualCov.gllvm to "x", in line with getEnvironCov.gllvm
* New vignette for the correlation structures of random effects and latent variables.

## Bug Fixes
* Bug fixed for calculating residual covariances of quadratic concurrent ordination

Version 1.4.7
=============

* Added "fungi" dataset by Abrego et al. 2022
* Added "kelpforest" dataset by Reed and Miller 2023
* New vignette for phylogenetic random effects
* New vignette for percent cover data analysis
* Function for calculating and plotting variance partitioning (varPartitioning.gllvm and plotVP)

Version 1.4.6
=============

* Added a 'getLoadings' function for retrieving species' loadings
* Added 'fac.center' argument in ordiplot to plot canonical coefficients of binary variables as points
* Added a simple plotting function for the gllvm summary
* Improved scaling for ordiplot with quadratic model and with biplot = FALSE
* optima.gllvm and tolerances.gllvm for num.lv now correctly provide tolerances w.r.t. the scaled LV
* Improved starting values for models with 'randomB'
* 'which.Xcoef' in coefplot.gllvm now also works for fourth-corner models
* Added intercept if beta0com=TRUE to coefplot.gllvm for fourth-corner models

## Bug Fixes
* Bug fixed that prevented increasing he point size of sites in ordiplot with symbols = TRUE
* Bug fixed in optima.gllvm for models with a single LV

Version 1.4.5
=============
* Separated "n.init" functionality into gllvm.iter.R
  * Prep for parallelisation
  * Enabled parallelisation (see TMB::openmp)
* Largely vectorized "residuals.gllvm", and residuals in "gllvm.aux"
* Added covariance of random effects to summary
* In preparation of emmeans support: moved the design matrix in "lv.X" to "lv.X.design". "lv.X" now stores the original supplied data.frame

## Bug Fixes
* Bug in ZINB fixed

Version 1.4.4
=============
* Removed "dependent.row" feature
* Added possibility for multiple random row intercepts
* Added possibility for (correlated) random species random effects
  * Can be plotted with "randomCoefPlot"
* Added possibility to Phylogenetically structure the random species effects
  * Phylogenetic signal parameter is included as object$params$rho.sp
  * Can be covariate specific
* num.RR and num.lv.c can now be larger than the number of predictors if randomB!=FALSE
* Added "iid" option for "randomB"
* Added a "getEnvironCov" function to extract species associations due to random covariate effects

Version 1.4.3
==============
* For CRAN release 1.4.3 see updates for versions 1.4.2 and 1.4.3

## Bug Fixes
* Bug in correlated row effects fixed
* Bug in getPredictErr for models fitted with LA fixed, and it returns now prediction errors for random slopes of X covariates as well
* Bug in randomCoefplot fixed

Version 1.4.2
==============

### New Features
* Added a correction factor to the second partial derivatives of the canonical coefficients for concurrent and constrained ordination
* Added `randomCoefPlot` functionality of constrained and concurrent ordination models with random slopes. Currently not supported for models with quadratic responses
* Summary now provides the possibility to calculate wald statistics across LVs or predictors for concurrent and constrained ordination
* `coef` now renames parameter estimates with more intuitive names and allows to subset the parameter list with names
* Tweedie power parameter is estimated now if set to NULL in `gllvm. 
* VA support for Zero-inflated poisson distribution
* Zero-inflated negative-binomial distribution added
* Binomial (Ntrials>1) support added (previously only Bernoulli)
* Now allowed to have (some) NAs in the response data

## Bug Fixes
* Fixed an issue with structured row-effects in concurrent and constrained ordination
* Fixed a bug that prevented plotting prediction regions for constrained ordination with structured row-effects
* No standard errors should be returned by optima.gllvm and tolerances.gllvm with randomB != FALSE
* Species names were in the original order with order = TRUE in RandomCoefPlot
* Fixed an issue that arose when {0,1} bounded parameters reached the bounds
* Various bug fixes for constrained/concurrent ordination with random intercepts and random slopes
* Bug in predictions with structured row intercepts was fixed, see issue #86

Version 1.4.1
==============

### New Features
* Computational stability of random slopes for constr. and concr. ordination significantly improved
* Computational stability of quadratic model significantly improved
* Unstructured VA covariance matrix for quadratic models with random intercepts
* Added example for se.gllvm

## Bug Fixes
* Bugfix in random slopes for concr. ordination with LV-specific variances and random row intercepts
* Bugfix for quadratic model with Poisson, NB, gamma, or exponential responses
* Bugfix in starting values for constrained and concurrent quadratic model

### Bug Fixes
* Valgrind error fixed

Version 1.4.0
==============

### New Features
* For CRAN release 1.4.0's new features see features described for versions 1.3.2-1.3.3

### Bug Fixes
* For bug fixes to CRAN release 1.4.0 see versions 1.3.2-1.3.3

Version 1.3.3
==============

### New Features
* The n.init option has been improved, so that it stops if no improved fit has been found after n.init.max (defaults to 10) iterations.
* Row names from the data now carry over to the site scores, so that they can be displayed in ordiplot

### Bug Fixes
* Memory allocation problem in development version fixed

* Diagonal elements of loading matrix 'theta' fixed for fourth corner model

* Bug in 'predict' for random slopes fixed, occurred when new x-covariate values were given

Version 1.3.2
==============

### New Features
* Ordination with predictors (num.RR,num.lv.c) is now implemented with constrained optimization routines (alabama,nloptr) as long as the canonical coefficients are treated as fixed-effects. This follows from the necessary identifiability constraints. 

* The reduced-rank approximated predictor slopes of a multivariate regression can now be plotted (with confidence intervals) using coefplot. Not available yet for quadratic effects.

* Separate checks are put in place to warn users if the constraints on the canonical coefficients (orthogonality of the columns) have not converged.

* Separate checks are put in place to warn users if the coefficients of a quadratic model have not converged

* Canonical coefficients in ordination with predictors (num.RR,num.lv.c) can now be treated as random-effects using the 'randomB' argument. For the moment, all need to be either random or fixed, no mixing. Prediction intervals can be retrieved with the getPredictErr function.

* An extended version of the spider dataset has been made available

* Added an option to magnify the x-axis labels in coefplot

* Site names present as row labels in the response data are now shown in the ordination plot


### Bug Fixes
* The order of the quadratic coefficients was  wrong when num.RR, num.lv, and num.lv.c were all used in the same model.

*  Fixed a bug in the calculation of starting values for constrained ordination (num.RR) where the residuals were not re-calculated if num.lv.c>0

*  Fixed a bug in coefplot for when only one predictor was included in the model

*  Fixed a bug that would prevent using a gllvm with quadratic response model as starting values for another model

* Changed import/export of various functions as requested in github issue #65

* Various minor tweaks to the summary function

Version 1.3.1
==============

### New Features

* Structured row parameters are implemented, including a possibility for between or within group correlations for random row effects.

* Constrained ordination model is implemented.

* NB and binomial (with probit and logit) response model implemented using extended variational approximation method.

### Bug Fixes

* Vignettes are removed from the CRAN version of the package, can be seen at the package's website only.


Version 1.3.0
==============

### New Features

* Quadratic latent variables allowed, that is term - u_i'D_j u_i can be included in the model using 'quadratic = TRUE'. 
  In addition, functions 'optima()', 'tolerances()' and 'gradient.length()' included.

* Beta response distribution implemented using Laplace approximation and extended variational approximation method.

* Tweedie response model implemented using extended variational approximation method.

* Ordinal model works now for 'num.lv=0'.

* Residual covariance adjustment added for gaussian family.

### Bug Fixes

* Estimation of the variances of random slopes of the X covariates didn't work properly when 'row.eff = FALSE' or 'row.eff = "fixed"'.

* Problems occurred in calculation of the starting values for ordinal model.

* Problems occurred in predict() and residuals(), when random slopes for X covariates were included.

* Problems occurred in predict() when new X covariates were given.

* Problems occurred in predictLVs() for fourth corner models.

Version 1.3.1
==============

### New Features

* Structured row parameters are implemented, including a possibility for between or within group correlations for random row effects.

* Constrained ordination model is implemented.

* NB and binomial (with probit and logit) response model implemented using extended variational approximation method.

### Bug Fixes

* Vignettes are removed from the CRAN version of the package, can be seen at the package's website only.


Version 1.3.0
==============

### New Features

* Quadratic latent variables allowed, that is term - u_i'D_j u_i can be included in the model using 'quadratic = TRUE'. 
  In addition, functions 'optima()', 'tolerances()' and 'gradient.length()' included.

* Beta response distribution implemented using Laplace approximation and extended variational approximation method.

* Tweedie response model implemented using extended variational approximation method.

* Ordinal model works now for 'num.lv=0'.

* Residual covariance adjustment added for gaussian family.

### Bug Fixes

* Estimation of the variances of random slopes of the X covariates didn't work properly when 'row.eff = FALSE' or 'row.eff = "fixed"'.

* Problems occurred in calculation of the starting values for ordinal model.

* Problems occurred in predict() and residuals(), when random slopes for X covariates were included.

* Problems occurred in predict() when new X covariates were given.

* Problems occurred in predictLVs() for fourth corner models.

