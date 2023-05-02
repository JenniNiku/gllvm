Version 1.4.2
* Added a correction factor to the second partial derivatives of the canonical coefficients for concurrent and constrained ordination
* Added `randomCoefPlot` functionality of constrained and concurrent ordination models with random slopes. Currently not supported for models with quadratic responses
* `coef` now renames parameter estimates with more intuitive names and allows to subset the parameter list with names
* Tweedie power parameter is estimated now if set to NULL in `gllvm. 

## Bug Fixes
* Fixed an issue with structured row-effects in concurrent and constrained ordination
* Fixed a bug that prevented plotting prediction regions for constrained ordination with structured row-effects

Version 1.4.1
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

