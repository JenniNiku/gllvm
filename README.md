# gllvm

`gllvm` is an R package for analysing multivariate ecological data with Generalized Linear Latent Variable Models (GLLVM).
Estimation is performed using maximum likelihood estimation, together with either variational approximation (VA) or Laplace approximation (LA) method to approximate the marginal likelihood.

# Installation

From CRAN you can install the package using:
```
install.packages("gllvm")
```
Or the development version of `gllvm` from github with the help of `devtools` package using:
```
devtools::install_github("JenniNiku/gllvm")
```

# Getting started

For getting started with `gllvm` we recommend to read vignette [Analysing multivariate abundance data using gllvm](https://jenniniku.github.io/gllvm/articles/vignette1.html)
or introductions for using `gllvm` for [ordination](https://jenniniku.github.io/gllvm/articles/vignette3.html) and for [analysing species correlations](https://jenniniku.github.io/gllvm/articles/vignette4.html).


<!-- badges: start -->
[![R-CMD-check](https://github.com/BertvanderVeen/gllvm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BertvanderVeen/gllvm/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/gllvm)](https://CRAN.R-project.org/package=gllvm)
<!-- badges: end -->
