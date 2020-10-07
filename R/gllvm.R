#' @title Generalized Linear Latent Variable Models
#' @description Fits generalized linear latent variable model for multivariate data. The model can be fitted using Laplace approximation method or variational
#' approximation method.
#'
#' @param y (n x m) matrix of responses.
#' @param X matrix or data.frame of environmental covariates.
#' @param TR matrix or data.frame of trait covariates.
#' @param data data in long format, that is, matrix of responses, environmental and trait covariates and row index named as ’id’. When used, model needs to be defined using formula. This is alternative data input for y, X and TR.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param num.lv  number of latent variables, d, in gllvm model. Non-negative integer, less than number of response variables (m). Defaults to 2.
#' @param family  distribution function for responses. Options are \code{poisson(link = "log")}, \code{"negative.binomial"} (with log link), \code{binomial(link = "probit")} (and also \code{binomial(link = "logit")} when \code{method = "LA"}), zero inflated poisson (\code{"ZIP"}), \code{gaussian(link = "identity")}, \code{"gamma"} (with log link), \code{"exponential"} (with log link), Tweedie (\code{"tweedie"}) (with log link, only with \code{"LA"}-method) and \code{"ordinal"} (only with \code{"VA"}-method).
#' @param method  model can be fitted using Laplace approximation method (\code{method = "LA"}) or variational approximation method (\code{method = "VA"}). Defaults to \code{"VA"}.
#' @param row.eff  \code{FALSE}, \code{fixed} or \code{"random"}, Indicating whether row effects are included in the model as a fixed or as a random effects. Defaults to \code{FALSE} when row effects are not included.
#' @param sd.errors  logical. If \code{TRUE} (default) standard errors for parameter estimates are calculated.
#' @param offset vector or matrix of offset terms.
#' @param la.link.bin link function for binomial family if \code{method = "LA"}. Options are "logit" and "probit.
#' @param Power fixed power parameter in Tweedie model. Scalar from interval (1,2). Defaults to 1.1.
#' @param seed a single seed value, defaults to \code{NULL}.
#' @param plot  logical, if \code{TRUE} ordination plots will be printed in each iteration step when \code{TMB = FALSE}. Defaults to \code{FALSE}.
#' @param zeta.struc Structure for cut-offs in the ordinal model. Either "common", for the same cut-offs for all species, or "species" for species-specific cut-offs. For the latter, classes are arbitrary per species, each category per species needs to have at least one observations. Defaults to "species".
#' @param randomX  formula for species specific random effects of environmental variables in fourth corner model. Defaults to \code{NULL}, when random slopes are not included.
#' @param dependent.row logical, whether or not random row effects are correlated (dependent) with the latent variables. Defaults to \code{FALSE} when correlation terms are not included.
#' @param beta0com logical, if \code{FALSE} column-specific intercepts are assumed. If \code{TRUE}, a common intercept is used which is allowed only for fourth corner models.
#' @param scale.X if \code{TRUE}, covariates are scaled when fourth corner model is fitted.
#' @param control A list with the following arguments controlling the optimization:
#' \itemize{
#'  \item{\emph{reltol}: }{ convergence criteria for log-likelihood, defaults to 1e-8.}
#'  \item{\emph{TMB}: }{ logical, if \code{TRUE} model will be fitted using Template Model Builder (TMB). TMB is always used if \code{method = "LA"}.  Defaults to \code{TRUE}.}
#'  \item{\emph{optimizer}: }{ if \code{TMB=TRUE}, log-likelihood can be optimized using \code{"\link{optim}"} (default) or \code{"\link{nlminb}"}.}
#'  \item{\emph{max.iter}: }{ maximum number of iterations when \code{TMB = FALSE}, defaults to 200.}
#'  \item{\emph{maxit}: }{ maximum number of iterations within \code{optim} function, defaults to 1000.}
#'  \item{\emph{trace}: }{ logical, if \code{TRUE} in each iteration step information on current step will be printed. Defaults to \code{FALSE}. Only with \code{TMB = FALSE}.}
#' }
#' @param control.va A list with the following arguments controlling the variational approximation method:
#' \itemize{
#'  \item{\emph{Lambda.struc}: }{ covariance structure of VA distributions for latent variables when \code{method = "VA"}, "unstructured" or "diagonal".}
#'  \item{\emph{Ab.struct}: }{ covariance structure of VA distributions for random slopes when \code{method = "VA"}, "unstructured" or "diagonal".}
#'  \item{\emph{diag.iter}: }{ non-negative integer which is used to speed up the updating of variational (covariance) parameters in VA method. Defaults to 5 if \code{TMB = FALSE}. If \code{TMB = TRUE} either 0 or 1.}
#'  \item{\emph{Ab.diag.iter}: }{ As above, but for variational covariance of random slopes.}
#'  \item{\emph{Lambda.start}: }{ starting values for variances in VA distributions for latent variables, random row effects and random slopes in variational approximation method. Defaults to 0.1.}
#' }
#' @param control.start A list with the following arguments controlling the starting values:
#' \itemize{
#'   \item{\emph{starting.val}: }{ starting values can be generated by fitting model without latent variables, and applying factorial analysis to residuals to get starting values for latent variables and their coefficients (\code{starting.val = "res"}). Another options are to use zeros as a starting values (\code{starting.val = "zero"}) or initialize starting values for latent variables with (n x num.lv) matrix. Defaults to \code{"res"}, which is recommended.}
#'   \item{\emph{n.init}: }{ number of initial runs. Uses multiple runs and picks up the one giving highest log-likelihood value. Defaults to 1.}
#'   \item{\emph{start.fit}: }{ object of class 'gllvm' which can be given as starting parameters for count data (poisson, NB, or ZIP).}
#'   \item{\emph{start.lvs}: }{ initialize starting values for latent variables with (n x num.lv) matrix. Defaults to \code{NULL}.}
#'   \item{\emph{jitter.var}: }{ jitter variance for starting values of latent variables. Defaults to 0, meaning no jittering.}
#'   \item{\emph{randomX.start}: }{ Starting value method for the random slopes. Options are \code{"zero"} and \code{"res"}. Defaults to \code{"res"}.}
#' }
#'
#' @details
#' Fits generalized linear latent variable models as in Hui et al. (2015 and 2017) and Niku et al. (2017).
#' Method can be used with two types of latent variable models depending on covariates. If only
#' site related environmental covariates are used, the expectation of response \eqn{Y_{ij}} is determined by
#'
#' \deqn{g(\mu_{ij}) = \eta_{ij} = \alpha_i + \beta_{0j} + x_i'\beta_j + u_i'\theta_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{u_i} are \eqn{d}-variate latent variables (\eqn{d}<<\eqn{m}), \eqn{\alpha_i} is an optional row effect
#' at site \eqn{i}, and it can be fixed or random effect, \eqn{\beta_{0j}} is an intercept term for species \eqn{j}, \eqn{\beta_j} and \eqn{\theta_j} are column
#' specific coefficients related to covariates and the latent variables, respectively.
#'
#' An alternative model is the fourth corner model (Brown et al., 2014, Warton et al., 2015) which will be fitted if also trait covariates
#' are included. The expectation of response \eqn{Y_{ij}} is
#'
#' \deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + x_i'(\beta_x + b_j) + TR_j'\beta_t + vec(B)*kronecker(TR_j,X_i) + u_i'\theta_j}
#'
#' where g(.), \eqn{u_i}, \eqn{\beta_{0j}} and \eqn{\theta_j} are defined as above. Vectors \eqn{\beta_x} and \eqn{\beta_t} are the main effects
#' or coefficients related to environmental and trait covariates, respectively, matrix \eqn{B} includes interaction terms. Vectors \eqn{b_j} are 
#' optional species-specific random slopes for environmental covariates.
#' The interaction/fourth corner terms are optional as well as are the main effects of trait covariates.
#'
#'
#' The method is sensitive for the choices of initial values of the latent variables. Therefore it is
#' recommendable to use multiple runs and pick up the one giving the highest log-likelihood value.
#' However, sometimes this is computationally too demanding, and default option
#' \code{starting.val = "res"} is recommended. For more details on different starting value methods, see Niku et al., (2018).
#'
#' Models are implemented using TMB (Kristensen et al., 2015) applied to variational approximation (Hui et al., 2017) and Laplace approximation (Niku et al., 2017).
#'
#' With ordinal family response classes must start from 0 or 1.
#'
#' \subsection{Distributions}{
#'
#'   Mean and variance for distributions are defined as follows.
#'\itemize{
#'   \item{For count data \code{family = poisson()}:} {Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}}, or}
#'   \item{ \code{family = "negative.binomial"}:}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}+\mu_{ij}^2\phi_j}, or}
#'   \item{ \code{family = "ZIP"}:}{ Expectation \eqn{E[Y_{ij}] = (1-p)\mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}(1-p)(1+\mu_{ij}p)}.}
#'
#'   \item{For binary data \code{family = binomial()}:}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}(1-\mu_{ij})}.}
#'
#'   \item{For positive continuous data \code{family = "gamma"}:}{Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}^2/\phi_j}, where \eqn{\phi_j} is species specific shape parameter.}
#'   
#'   \item{For non-negative  continuous data \code{family = "exponential"}:}{Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}^2}.}
#'   
#'   \item{For non-negative continuous or biomass data\code{family = "tweedie"}}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \phi_j*\mu_{ij}^\nu}, where \eqn{\nu} is a power parameter of Tweedie distribution. See details Dunn and Smyth (2005).}
#'
#'   \item{For ordinal data \code{family = "ordinal"}:}{ Cumulative probit model, see Hui et.al. (2016).}
#'   
#'   \item{For normal distributed data \code{family = gaussian()}:}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(y_{ij}) = \phi_j^2.}}
#' }
#' }
#'
#'@note If function gives warning: 'In f(x, order = 0) : value out of range in 'lgamma'', optimizer have visited an area where gradients become too big. It is automatically fixed by trying another step in the optimization process, and can be ignored if errors do not occur.
#'
#' @return An object of class "gllvm" includes the following components:
#'
#'
#'  \item{call }{function call}
#'  \item{logL }{log likelihood}
#'  \item{lvs }{latent variables}
#'  \item{params}{list of parameters
#'  \itemize{
#'    \item{theta }{ coefficients related to latent variables}
#'    \item{beta0 }{ column specific intercepts}
#'    \item{Xcoef }{ coefficients related to environmental covariates X}
#'    \item{B }{ coefficients in fourth corner model}
#'    \item{row.params }{ row-specific intercepts}
#'    \item{phi }{ dispersion parameters \eqn{\phi} for negative binomial or Tweedie family, probability of zero inflation for ZIP family, standard deviation for gaussian family or shape parameter for gamma family}
#'    \item{inv.phi }{ dispersion parameters \eqn{1/\phi} for negative binomial}
#'    }}
#'  \item{Power }{ power parameter \eqn{\nu} for Tweedie family}
#'  \item{sd }{ list of standard errors of parameters}
#'  \item{prediction.errors }{ list of prediction covariances for latent variables and variances for random row effects when method \code{"LA"} is used}
#'  \item{A, Ar }{ covariance matrices for variational densities of latent variables and variances for random row effects}
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Wesley Brooks, Riki Herliansyah, Francis K.C. Hui, Sara Taskinen, David I. Warton, Bert van der Veen
#' @references
#' Brown, A. M., Warton, D. I., Andrew, N. R., Binns, M., Cassis, G., and Gibb, H. (2014). The fourth-corner solution - using predictive models to understand how species traits interact with the environment. Methods in Ecology and Evolution, 5:344-352.
#'
#' Dunn, P. K. and Smyth, G. K. (2005).  Series evaluation of tweedie exponential dispersion model densities. Statistics and Computing, 15:267-280.
#'
#' Hui, F. K. C., Taskinen, S., Pledger, S., Foster, S. D., and Warton, D. I. (2015).  Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6:399-411.
#'
#' Hui, F. K. C., Warton, D., Ormerod, J., Haapaniemi, V., and Taskinen, S. (2017).  Variational approximations for generalized linear latent variable models. Journal of Computational and Graphical Statistics. Journal of Computational and Graphical Statistics, 26:35-43.
#'
#' Kasper Kristensen, Anders Nielsen, Casper W. Berg, Hans Skaug, Bradley M. Bell (2016). TMB: Automatic Differentiation and Laplace Approximation. Journal of Statistical Software, 70(5), 1-21.
#'
#' Niku, J., Warton,  D. I., Hui, F. K. C., and Taskinen, S. (2017). Generalized linear latent variable models for multivariate count and biomass data in ecology. Journal of Agricultural, Biological, and Environmental Statistics, 22:498-522.
#'
#' Niku, J., Brooks, W., Herliansyah, R., Hui, F. K. C., Taskinen, S., and Warton,  D. I. (2018). Efficient estimation of generalized linear latent variable models. PLoS One, 14(5):1-20.
#'
#' Warton, D. I., Guillaume Blanchet, F., O'Hara, R. B., Ovaskainen, O., Taskinen, S., Walker, S. C. and Hui, F. K. C. (2015). So many variables: Joint modeling in community ecology. Trends in Ecology & Evolution, 30:766-779.
#'
#'@seealso  \code{\link{coefplot.gllvm}}, \code{\link{confint.gllvm}}, \code{\link{ordiplot.gllvm}}, \code{\link{plot.gllvm}}, \code{\link{residuals.gllvm}}, \code{\link{summary.gllvm}}.
#' @examples
#'# Extract subset of the microbial data to be used as an example
#'data(microbialdata)
#'X <- microbialdata$Xenv
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, X, formula = ~ pH + Phosp, family = poisson())
#'ordiplot(fit)
#'coefplot(fit)
#'
#' \donttest{
#'## Load a dataset from the mvabund package
#'library(mvabund)
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- as.matrix(antTraits$env)
#'TR <- antTraits$traits
#'# Fit model with environmental covariates Bare.ground and Shrub.cover
#'fit <- gllvm(y, X, formula = ~ Bare.ground + Shrub.cover,
#'             family = poisson())
#'ordiplot(fit)
#'coefplot(fit)
#'
#'## Example 1: Fit model with two latent variables
#'# Using variational approximation:
#'fitv0 <- gllvm(y, family = "negative.binomial", method = "VA")
#'ordiplot(fitv0)
#'plot(fitv0, mfrow = c(2,2))
#'summary(fitv0)
#'confint(fitv0)
#'# Using Laplace approximation: (this line may take about 30 sec to run)
#'fitl0 <- gllvm(y, family = "negative.binomial", method = "LA")
#'ordiplot(fitl0)
#'
#'# Poisson family:
#'fit.p <- gllvm(y, family = poisson(), method = "LA")
#'ordiplot(fit.p)
#'# Use poisson model as a starting parameters for ZIP-model, this line may take few minutes to run
#'fit.z <- gllvm(y, family = "ZIP", method = "LA", control.start =list(start.fit = fit.p))
#'ordiplot(fit.z)
#'
#'
#'## Example 2: gllvm with environmental variables
#'# Fit model with two latent variables and all environmental covariates,
#'fitvX <- gllvm(formula = y ~ X, family = "negative.binomial")
#'ordiplot(fitvX, biplot = TRUE)
#'coefplot(fitvX)
#'# Fit model with environmental covariates Bare.ground and Shrub.cover
#'fitvX2 <- gllvm(y, X, formula = ~ Bare.ground + Shrub.cover,
#'  family = "negative.binomial")
#'ordiplot(fitvX2)
#'coefplot(fitvX2)
#'# Use 5 initial runs and pick the best one
#'fitvX_5 <- gllvm(y, X, formula = ~ Bare.ground + Shrub.cover,
#'  family = "negative.binomial", control.start=list(n.init = 5, jitter.var = 0.1))
#'ordiplot(fitvX_5)
#'coefplot(fitvX_5)
#'
#'## Example 3: Data in long format
#'# Reshape data to long format:
#'datalong <- reshape(data.frame(cbind(y,X)), direction = "long",
#'                    varying = colnames(y), v.names = "y")
#'head(datalong)
#'fitvLong <- gllvm(data = datalong, formula = y ~ Bare.ground + Shrub.cover,
#'                family = "negative.binomial")
#'
#'## Example 4: Fourth corner model
#'# Fit fourth corner model with two latent variables
#'fitF1 <- gllvm(y = y, X = X, TR = TR, family = "negative.binomial")
#'coefplot(fitF1)
#'# Fourth corner can be plotted also with next lines
#'#fourth = fitF1$fourth.corner
#'#library(lattice)
#'#a = max( abs(fourth) )
#'#colort = colorRampPalette(c("blue","white","red"))
#'#plot.4th = levelplot(t(as.matrix(fourth)), xlab = "Environmental Variables",
#'#              ylab = "Species traits", col.regions = colort(100),
#'#              at = seq( -a, a, length = 100), scales = list( x = list(rot = 45)))
#'#print(plot.4th)
#'
#'# Specify model using formula
#'fitF2 <- gllvm(y = y, X = X, TR = TR,
#'  formula = ~ Bare.ground + Canopy.cover * (Pilosity + Webers.length),
#'  family = "negative.binomial")
#'ordiplot(fitF2)
#'coefplot(fitF2)
#'
#'## Include species specific random slopes to the fourth corner model
#'fitF3 <- gllvm(y = y, X = X, TR = TR,
#'  formula = ~ Bare.ground + Canopy.cover * (Pilosity + Webers.length),
#'  family = "negative.binomial", randomX = ~ Bare.ground + Canopy.cover, 
#'  control.start = list(n.init = 3))
#'ordiplot(fitF3)
#'coefplot(fitF3)
#'
#'
#'## Example 5: Fit Tweedie model
#'# Load coral data
#'data(tikus)
#'ycoral <- tikus$abund
#'# Let's consider only years 1981 and 1983
#'ycoral <- ycoral[((tikus$x$time == 81) + (tikus$x$time == 83)) > 0, ]
#'# Exclude species which have observed at less than 4 sites
#'ycoral <- ycoral[-17, (colSums(ycoral > 0) > 4)]
#'# Fit Tweedie model for coral data (this line may take few minutes to run)
#'fit.twe <- gllvm(y = ycoral, family = "tweedie", method = "LA")
#'ordiplot(fit.twe)
#'
#'## Example 6: Random row effects
#'fitRand <- gllvm(y, family = "negative.binomial", row.eff = "random")
#'ordiplot(fitRand, biplot = TRUE)
#'}
#' @export
#'
#'@useDynLib gllvm, .registration = TRUE
#'@importFrom TMB MakeADFun
#'@importFrom graphics abline axis par plot segments text points boxplot panel.smooth lines polygon
#'@importFrom grDevices rainbow
#'@importFrom stats dnorm pnorm qnorm rnorm dbinom pbinom rbinom pnbinom rnbinom pexp rexp pgamma rgamma ppois rpois runif pchisq qchisq qqnorm lm AIC binomial constrOptim factanal glm model.extract model.frame model.matrix model.response nlminb optim optimHess reshape residuals terms BIC qqline sd formula ppoints quantile gaussian cov
#'@importFrom Matrix bdiag chol2inv diag
#'@importFrom MASS ginv polr
#'@importFrom mgcv gam predict.gam
#'@importFrom mvtnorm rmvnorm

gllvm <- function(y = NULL, X = NULL, TR = NULL, data = NULL, formula = NULL,
                  num.lv = 2, family, row.eff = FALSE,
                  offset = NULL, sd.errors = TRUE, method = "VA",
                  randomX = NULL, dependent.row = FALSE, beta0com = FALSE, zeta.struc="species",
                  plot = FALSE, la.link.bin = "probit",
                  Power = 1.1, seed = NULL, scale.X = TRUE, 
                  control = list(reltol = 1e-8, TMB = TRUE, optimizer = "optim", max.iter = 200, maxit = 1000, trace = FALSE), 
                  control.va = list(Lambda.struc = "unstructured", Ab.struct = "unstructured", diag.iter = 5, Ab.diag.iter=0, Lambda.start = c(0.1, 0.1, 0.1)),
                  control.start = list(starting.val = "res", n.init = 1, jitter.var = 0, start.fit = NULL, start.lvs = NULL, randomX.start = "res")
                  ) {
    constrOpt <- FALSE
    restrict <- 30
    term <- NULL
    datayx <- NULL
    
    fill_control = function(x){
      if (!("reltol" %in% names(x))) 
        x$reltol = 1e-8
      if (!("TMB" %in% names(x))) 
        x$TMB = TRUE
      if (!("optimizer" %in% names(x))) 
        x$optimizer = "optim"
      if (!("max.iter" %in% names(x))) 
        x$max.iter = 200
      if (!("maxit" %in% names(x))) 
        x$maxit = 1000
      if (!("trace" %in% names(x))) 
        x$trace = FALSE
      x
    }
    fill_control.va = function(x){
      if (!("Lambda.struc" %in% names(x))) 
        x$Lambda.struc = "unstructured"
      if (!("Ab.struct" %in% names(x))) 
        x$Ab.struct = "unstructured"
      if (!("diag.iter" %in% names(x))) 
        x$diag.iter = 5
      if (!("Ab.diag.iter" %in% names(x))) 
        x$Ab.diag.iter = 0
      if (!("Lambda.start" %in% names(x))) 
        x$Lambda.start = c(0.1, 0.1, 0.1)
      x
    }
    fill_control.start = function(x){
      if (!("starting.val" %in% names(x))) 
        x$starting.val = "res"
      if (!("n.init" %in% names(x))) 
        x$n.init = 1
      if (!("jitter.var" %in% names(x))) 
        x$jitter.var = 0
      if (!("start.fit" %in% names(x))) 
        x$start.fit = NULL
      if (!("start.lvs" %in% names(x))) 
        x$start.lvs = NULL
      if (!("randomX.start" %in% names(x))) 
        x$randomX.start = "res"
      x
    }
    control <- fill_control(control)
    control.va <- fill_control.va(control.va)
    control.start <- fill_control.start(control.start)
    
    reltol = control$reltol; TMB = control$TMB; optimizer = control$optimizer; max.iter = control$max.iter; maxit = control$maxit; trace = control$trace;
    Lambda.struc = control.va$Lambda.struc; Ab.struct = control.va$Ab.struct; diag.iter = control.va$diag.iter; Ab.diag.iter=control.va$Ab.diag.iter; Lambda.start = control.va$Lambda.start
    starting.val = control.start$starting.val; n.init = control.start$n.init; jitter.var = control.start$jitter.var; start.fit = control.start$start.fit; start.lvs = control.start$start.lvs; randomX.start = control.start$randomX.start
    
    if(!is.null(X)){
      if(!is.matrix(X) && !is.data.frame(X) ) 
        stop("X must be a matrix or data.frame.")
    }
    if(!is.null(TR)){
      if(!is.matrix(TR) && !is.data.frame(TR) ) 
        stop("TR must be a matrix or data.frame.")
    }

    if (!is.null(y)) {
      y <- as.matrix(y)
      if (is.null(X) && is.null(TR)) {
        datayx <- list(y)
        m1 <- model.frame(y ~ NULL, data = datayx)
        term <- terms(m1)
      } else if (is.null(TR)) {
        if (is.null(formula)) {
          ff <- formula(paste("~", "0", paste("+", colnames(X), collapse = "")))
          if (is.data.frame(X)) {
            datayx <- list(y = y, X = model.matrix(ff, X))
          } else {
            datayx <- list(y = y, X = X)
          }
          m1 <- model.frame(y ~ X, data = datayx)
          term <- terms(m1)
        } else {
          datayx <- data.frame(y, X)
          m1 <- model.frame(formula, data = datayx)
        }
        term <- terms(m1)
      } else {
        term <- NULL
      }
      p <- NCOL(y)
      n <- NROW(y)
      if (p == 1)
        y <- as.matrix(y)
    } else {
      if (!is.null(data)) {
        if (is.null(formula))
          stop("Define formula when 'data' attribute is used.")
        if ("id" %in% colnames(data)) {
          id <- data[, "id"]
          n <- max(id)
          p <- dim(data)[1] / n
        } else {
          n <- NROW(data)
          p <- 1
          id <- 1:n
        }
      }

      cl <- match.call()
      mf <- match.call(expand.dots = FALSE)
      m <- match(c("formula", "data", "na.action"), names(mf), 0)
      mf <- mf[c(1, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      term <- attr(mf, "terms")
      abundances <- model.response(mf, "numeric")
      if (any(is.na(abundances)))
        stop("There are NA values in the response.")
      y <- abundances
      #
      X <- model.matrix(term, mf)

      atr <- c(attr(X, "assign"))
      if (sum(atr) > 0) {
        X <- X[, (atr > 0) * 1:ncol(X)]
      } else{
        X <- NULL
      }

      if (NCOL(y) == 1 &&
          !is.null(data)) {
        y <- matrix(y, n, p)
        colnames(y) <- paste("y", 1:p, sep = "")
      }
      try(
      if (is.null(X)) {
        datayx <- data.frame(y = y)
      } else {
        datayx <- data.frame(y = y, X = X)
      }, silent = TRUE)

      if (!is.null(data)) {
        frame1 <- mf
        X <- TR <- NULL
        if (length(attr(term, "term.labels")) > 0) {
          datax <- frame1[, colnames(frame1)!="y"]
          colnames(datax) <- colnames(frame1)[colnames(frame1)!="y"]
          #datax <- frame1[, attr(term, "term.labels")[attr(term, "order") == 1]]
          # colnames(datax) <- attr(term, "term.labels")[attr(term, "order") == 1]

          for (k in 1:ncol(datax)) {
            lngth <- NULL
            namek <- colnames(datax)[k]
            for (i in 1:n) {
              lngth <- c(lngth, length(unique(datax[(id == i), k])))
            }
            if (max(lngth) == 1) {
              if (!is.null(X))
                X <- data.frame(X, datax[1:n, k])
              else
                X <- data.frame(datax[1:n, k])

              colnames(X)[ncol(X)] <- namek
            } else {
              if (!is.null(TR)){
                TR <- data.frame(TR, datax[id == 1, k])
              } else{
                TR <- data.frame(datax[id == 1, k])
                }
              colnames(TR)[ncol(TR)] <- namek
            }
          }
        }
      }
    }
    p <- NCOL(y)
    n <- NROW(y)
    if (p == 1)
      y <- as.matrix(y)

    if (class(family) == "family") {
      la.link.bin <- family$link
      family <- family$family
    }

    if(any(colSums(y)==0))
      warning("There are responses full of zeros. \n");

    if(row.eff %in% c("fixed", "random", TRUE)){
      if(p<2)
        stop("There must be at least two responses in order to include row effects. \n");
      if(any(rowSums(y)==0))
        warning("There are rows full of zeros in y. \n");
      }


    if (row.eff == "random" && family == "ordinal" && TMB==FALSE) {
      stop("Random row effect model is not implemented for ordinal family. \n")
    }
    if (method == "LA" && family == "ordinal") {
      cat("Laplace's method cannot yet handle ordinal data, so VA method is used instead. \n")
      method <- "VA"
    }

    if (method == "LA" && !TMB) {
      cat("Laplace's method is not implemented without TMB, so 'TMB = TRUE' is used instead. \n")
      TMB = TRUE
    }
    if (method == "VA" && (family == "tweedie" || family == "ZIP")) {
      cat("VA method cannot handle", family, " family, so LA method is used instead. \n")
      method <- "LA"
    }
    if (p < 3 && !is.null(TR)) {
      stop("Fourth corner model can not be fitted with less than three response variables.\n")
    }
    if (row.eff == "random" && !TMB) {
      cat("Random row effect model is not implemented without TMB, so 'TMB = TRUE' is used instead. \n")
      TMB <- TRUE
    }
    
    if (family == "gaussian" && !TMB) {
      TMB <- TRUE
      cat("Only TMB implementation available for ", family, " family, so 'TMB = TRUE' is used instead. \n")
    }


    if (!is.null(start.fit)) {
      if (class(start.fit) != "gllvm")
        stop("Only object of class 'gllvm' can be given as a starting parameters.")

      if (!(family %in% c("poisson", "negative.binomial", "ZIP")))
        stop("Starting parameters can be given only for count data.")

    }
    #  if(num.lv>=p){ stop("Number of latent variables (",num.lv,") must be less than number of response variables (",p,").");}


    if (is.null(offset))
      O <- matrix(0, nrow = n, ncol = p)
    else if (NCOL(offset) == 1)
      O <- matrix(rep(offset), nrow = n, ncol = p)
    else
      O <- as.matrix(offset)

    if (is.matrix(start.lvs)) {
      starting.val <- "random"
      if (ncol(start.lvs) != num.lv || nrow(start.lvs) != n)
        stop("Given starting value matrix for latent variables has a wrong dimension.")
      colnames(start.lvs) <-  paste("LV",1:num.lv, sep = "")
    }
    n.i <- 1

    out <- list( y = y, X = X, TR = TR, data = datayx, num.lv = num.lv,
        method = method, family = family, row.eff = row.eff, randomX = randomX, n.init = n.init,
        sd = FALSE, Lambda.struc = Lambda.struc, TMB = TMB, terms = term, beta0com = beta0com)

    if (family == "binomial") {
      if (method == "LA")
        out$link <- la.link.bin
      if (method == "VA")
        out$link <- "probit"
    }
    out$offset <- offset


    if (TMB) {
      trace = FALSE
      if (row.eff == TRUE)
        row.eff <- "fixed"
      if (!is.null(TR)) {
        fitg <- trait.TMB(
            y,
            X = X,
            TR = TR,
            formula = formula,
            num.lv = num.lv,
            family = family,
            Lambda.struc = Lambda.struc,
            row.eff = row.eff,
            reltol = reltol,
            seed = seed,
            maxit = maxit,
            start.lvs = start.lvs,
            offset = O,
            sd.errors = sd.errors,
            trace = trace,
            link = la.link.bin,
            n.init = n.init,
            start.params = start.fit,
            optimizer = optimizer,
            starting.val = starting.val,
            method = method,
            Power = Power,
            diag.iter = diag.iter,
            Ab.diag.iter = Ab.diag.iter,
            Ab.struct = Ab.struct,
            dependent.row = dependent.row,
            Lambda.start = Lambda.start,
            jitter.var = jitter.var,
            randomX = randomX,
            randomX.start = randomX.start,
            beta0com = beta0com, 
            scale.X = scale.X,
            zeta.struc = zeta.struc
        )
        out$X <- fitg$X
        out$TR <- fitg$TR

      } else {
        fitg <- gllvm.TMB(
            y,
            X = X,
            formula = formula,
            num.lv = num.lv,
            family = family,
            method = method,
            Lambda.struc = Lambda.struc,
            row.eff = row.eff,
            reltol = reltol,
            seed = seed,
            maxit = maxit,
            start.lvs = start.lvs,
            offset = O,
            sd.errors = sd.errors,
            trace = trace,
            link = la.link.bin,
            n.init = n.init,
            restrict = restrict,
            start.params = start.fit,
            optimizer = optimizer,
            starting.val = starting.val,
            Power = Power,
            diag.iter = diag.iter,
            dependent.row = dependent.row,
            Lambda.start = Lambda.start,
            jitter.var = jitter.var,
            zeta.struc = zeta.struc
          )
      }

      out$X.design <- fitg$X.design
      out$TMBfn = fitg$TMBfn
      out$logL <- fitg$logL
      
      if (num.lv > 0)
        out$lvs <- fitg$lvs
      out$X <- fitg$X
      
      

      out$params <- fitg$params
      if (sd.errors) {
        out$sd <- fitg$sd
      }
      if (family == "tweedie") {
        out$Power <- fitg$Power
      }
      if(family == "ordinal"){
        out$zeta.struc = zeta.struc
      }
      if (method == "VA") {
        out$A <- fitg$A
        out$Ar <- fitg$Ar
      }
      if (!is.null(randomX)) {
        out$corr <- fitg$corr
        out$Xrandom <- fitg$Xrandom
        out$Ab <- fitg$Ab
      }
      out$start <- fitg$start

    } else {
      if (row.eff == "fixed")
        row.eff <- TRUE

      fitg <- gllvm.VA(
          y,
          X = X,
          TR = TR,
          family = family,
          formula = formula,
          num.lv = num.lv,
          max.iter = max.iter,
          eps = reltol,
          row.eff = row.eff,
          Lambda.struc = Lambda.struc,
          trace = trace,
          plot = plot,
          sd.errors = sd.errors,
          start.lvs = start.lvs,
          offset = O,
          maxit = maxit,
          diag.iter = diag.iter,
          seed = seed,
          n.init = n.init,
          restrict = restrict,
          constrOpt = constrOpt,
          start.params = start.fit,
          starting.val = starting.val,
          Lambda.start = Lambda.start,
          jitter.var = jitter.var
        )
      out$logL <- fitg$logLik
      if (num.lv > 0)
        out$lvs <- fitg$lvs
      out$X <- fitg$X
      out$TR <- fitg$TR
      out$X.design <- fitg$X.design
      out$params <- fitg$coef
      if (sd.errors) {
        out$sd <- fitg$sd
      }
      out$Lambda.struc <- fitg$Lambda.struc
      out$A <- fitg$Lambda
      out$start <- fitg$start
    }
    if (family == "negative.binomial")
      out$params$inv.phi <- 1 / out$params$phi
    if (is.infinite(out$logL)){
      warning("Algorithm converged to infinity, try other starting values or different method.")
      cat("Algorithm converged to infinity, try other starting values or different method. \n")
      }
    out$formula <- fitg$formula
    if (is.null(out$terms))
      out$terms <- fitg$terms
    if (is.finite(out$logL) && !is.null(TR) && NCOL(out$TR)>0 && NCOL(out$X)>0) {
      out$fourth.corner <- try(getFourthCorner(out),silent = TRUE)
    }
    if (is.finite(out$logL) && row.eff == "random" && FALSE){
      if(method == "LA"){
        if(abs(out$params$sigma)<0.02)
          cat("Random row effects ended up to almost zero. Might be a false convergence or local maxima. You can try simpler model, less latent variables or change the optimizer. \n")
      } else{
        if(abs(out$params$sigma)<0.02 && max(abs(out$params$sigma-sqrt(out$Ar))) < 1e-3)
          cat("Random row effects ended up to almost zero. Might be a false convergence or local maxima. You can try simpler model, less latent variables or change the optimizer. \n")
      }
    }

    out$Hess = fitg$Hess
    out$prediction.errors = fitg$prediction.errors
    out$call <- match.call()
    class(out) <- "gllvm"
    return(out)
  }
