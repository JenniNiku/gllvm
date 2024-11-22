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

Other available vignettes are:  [Analysing microbial community data](https://CRAN.R-project.org/package=gllvm/vignettes/vignette2.html),
[How to use the quadratic response model](https://CRAN.R-project.org/package=gllvm/vignettes/vignette5.html),
[Ordination with predictors](https://CRAN.R-project.org/package=gllvm/vignettes/vignette6.html), [Analysing percent cover data](https://jenniniku.github.io/gllvm/articles/vignette8.html) and 
[Structured and correlated random effects and latent variables](https://jenniniku.github.io/gllvm/articles/vignette9.html).

# Citation
The `citation` function in R provides information on how to cite the methods in this package. Please remember to cite the software (version) separately from any relevent research articles to provide the appropriate credit to all associated contributors. The reference for the software package is: Niku, J., Brooks, W., Herliansyah, R., Hui, F. K. C., Korhonen, P., Taskinen, S., van der Veen, B., and Warton, D. I.
  (YYYY). gllvm: Generalized Linear Latent Variable Models.R package version XXX, where YYYY represents the publication date of the used version of the package represented by XXX.

## References

[Hui, F.K.C., Warton, D., Ormerod, J., Haapaniemi, V., & Taskinen, S. (2017). Variational approximations for generalized linear latent variable models. Journal of Computational and Graphical Statistics, 26(1), 35 - 43.](https://www.tandfonline.com/doi/abs/10.1080/10618600.2016.1164708)

[Niku, J., Warton, D., Hui, F.K.C., & Taskinen, S. (2017). Generalized linear latent variable models for multivariate count and biomass data in ecology. Journal of Agricultural, Biological and Environmental Statistics, 22(4), 498 - 522.](https://link.springer.com/article/10.1007/s13253-017-0304-7)

[Niku, J., Hui, F.K.C., Taskinen, S., & Warton, D. (2019). gllvm: Fast analysis of multivariate abundance data with generalized linear latent variable models in r. Methods in Ecology and Evolution, 10(12), 2173 - 2182.](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13303)

[Niku, J., Brooks, W., Herliansyah, R., Hui, F.K.C., Taskinen, S., & Warton, D. (2019). Efficient estimation of generalized linear latent variable models. PloS one, 14(5), e0216129.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0216129)

[Niku, J., Hui, F. K. C., Taskinen, S., and Warton, D. I. (2021). Analyzing environmental-trait interactions in ecological communities with fourth-corner latent variable models.       Environmetrics, 32(6), 1-17.](https://doi.org/10.1002/env.2683)

[van der Veen, B., Hui, F.K.C., Hovstad, K.A., Solbu, E.B., & O'Hara, R.B. (2021). Model-based ordination for species with unequal niche widths. Methods in Ecology and Evolution, 12(7), 1288 - 1300.](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13595)

[van der Veen, B., Hui, F. K. C., Hovstad, K.A., and O'Hara, R.B. (2023). Concurrent ordination: simultaneous unconstrained and constrained latent variable modelling. Methods in Ecology and Evolution, 14(2), 683-695.](https://doi.org/10.1111/2041-210X.14035)

[van der Veen, B. and O'Hara, R.B. (2024). Fast fitting of phylogenetic mixed effects models. arxiv.](https://www.arxiv.org/abs/2408.05333)

[Korhonen, P., Hui, F. K. C., Niku, J., and Taskinen, S. (2023). Fast and universal estimation of latent variable models using extended variational approximations. Statistics and Computing, 33(1), 1-16.](https://doi.org/10.1007/s11222-022-10189-w)
