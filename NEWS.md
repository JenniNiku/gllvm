Version 1.3.dev
==============

### New Features

* Structured row parameters are implemented, including a possibility for between or within group correlations for random row effects.

* Constrained ordination model is implemented.

* NB and binomial (with probit and logit) response model implemented using extended variational approximation method.

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

