Version 1.2.4
==============

### New Features

* Quadratic latent variables allowed, that is term - u_i'D_j u_i can be included in the model. 
  In addition, functions 'optima()', 'tolerances()' and 'gradient.length()' included.

* Beta distribution implemented for the responses.

* Ordinal model works now for 'num.lv=0'.

* Residual covariance adjustment added for gaussian family.

### Bug Fixes

* Estimation of the variances of random slopes of the X covariates didn't work properly when 'row.eff = FALSE' or 'row.eff = "fixed"'.

* Problems occurred in calculation of the starting values for ordinal model.

