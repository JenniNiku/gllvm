#' @title Generalized Linear Latent Variable Models.
#' @description Fit generalized linear latent variable models for multivariate data. The models can be fitted using Laplace approximation method or Variational
#' approximation method.
#'
#' @param y (n x m) matrix of responses. If data is in long format, \code{y} is a data.frame including a row index vector, a column index vector and a vector of response variables in this order.
#' @param X Matrix of environmental covariates.
#' @param TR Matrix of trait covariates. Can be used only with VA method
#' @param data Matrix of environmental and trait covariates. Can be used only if responses y are in long format.
#' @param n.X Number of  environmental covariates in data matrix above.
#' @param n.TR Number of  trait covariates in data matrix above.
#' @param num.lv  Number of latent variables, d in GLLVM model. Non-negative integer. Defaults to 2.
#' @param family  Distribution function to be used in the model. Options are "poisson" (with log link), "negative.binomial" (with log link) and "ZIP" (with log link), "binomial" (with logit/cloglog link when method="LA" and probit link when method="VA"), Tweedie (with log link, only with "LA"-method), "ordinal" (only with "VA"-method).
#' @param method  Model can be fitted using Laplace Approximation method (\code{method="LA"}) or Variational Approximation method (\code{method="VA"}). Defaults to \code{"VA"}.
#' @param row.eff  Logical, indicating whether row effects are included in the model. Defaults to \code{FALSE}.
#' @param get.fourth Logical. If \code{TRUE} (default) fourth corner/interaction terms of environmental and trait covariates are included in the model.
#' @param get.trait Logical. If \code{TRUE} (default) main effects of trait covariates are included in the model.
#' @param start.lvs Initialize starting values for latent variables with (n x \code{num.lv}) matrix. Defaults to NULL when starting values are generated randomly.
#' @param start.params object of class 'gllvm' which can be given as starting parameters for count data (poisson, NB, or ZIP).
#' @param sd.errors  Logical. If \code{TRUE} (default) standard errors for parameter estimates are calculated.
#' @param n.init Number of initial runs. Uses multiple runs and picks the one giving highest log-likelihood value. This is recommendable because the method is quite sensitive to initial values of the latent variables. Defaults to 1.
#' @param offset Vector or matrix of offset terms
#' @param la.phi.upd  Defaults to "phi" when dispersion parameters with LA method are updated using parametrization \eqn{\phi}. If "inv.phi", dispersion parameters are updated using parametrization \eqn{1/\phi}.
#' @param la.link.bin Link function for binomial family if \code{method="LA"}. Options are "logit" and "cloglog".
#' @param constrOpt Logical, if \code{TRUE} parameters are estimated by constraining linear predictors in optimization. This is recommendable if algorithm is very unstable in different runs. Defaults to \code{FALSE}. Cannot be used when family is "tweedie" or "ZIP". Constraint for the absolute value of the linear predictor is defined by \code{restrict}.
#' @param fixed.power Logical, if \code{TRUE} (default) the power parameter in Tweedie model is fixed.
#' @param Power Fixed power parameter if \code{fixed.power=TRUE} or starting value for power parameter if \code{fixed.power=FALSE} in Tweedie model. Scalar from interval (1,2). Defaults to 1.5.
#' @param Lambda.struc  Covariance structure for latent variables when \code{method = "VA"}, "unstructured" or "diagonal".
#' @param diag.iter  Non-negative integer which is used to speed up the updating of variational parameters Lambda in VA method. Defaults to 5.
#' @param trace  Logical, if \code{TRUE} in each iteration step information on current step will be printed. Defaults to \code{FALSE}.
#' @param plot  Logical, if \code{TRUE} ordination plots will be printed on each iteration step. Defaults to \code{FALSE}.
#' @param eps  Convergence criteria for log likelihood, defaults to 1e-4.
#' @param max.iter Maximum number of iterations, defaults to 100.
#' @param maxit Maximum number of iterations within \code{optim} function, defaults to 1000.
#' @param seed a single seed value, defaults to \code{NULL}.
#'
#' @details
#' Fits generalized linear latent variable models.
#' Method can be used with two types of latent variable models depending on covariates. If only
#' site related environmental covariates are used, the expectation of response \eqn{Y_{ij}} is determined by
#'
#' \deqn{g(\mu_{ij}) = \eta_{ij} = \alpha_i + \beta_{0j} + x_i'\beta_j + u_i'\theta_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{u_i} are \eqn{d}-variate latent variables, \eqn{\alpha_i} is an optional row effect
#' at site \eqn{i}, \eqn{beta_{0j}} is intercept term for species \eqn{j}, \eqn{\beta_j} and \eqn{\theta_j} are column
#' specific coefficients related to covariates and the latent variables, respectively.
#' If also trait covariates are included, the expectation of response \eqn{Y_{ij}} is
#'
#'\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + x_i'\beta_x + TR_j'\beta_t + vec(B)*kronecker(TR_j,X_i) + u_i'\theta_j}
#'
#' where g(.), \eqn{u_i}, \eqn{\beta_{0j}} and \eqn{\theta_j} are defined as above. Vectors \eqn{beta_x} and \eqn{beta_t} are main effects
#' or coefficients related to environmental and trait covariates, respectively, matrix \eqn{B} includes interaction terms.
#' This is the so called 'fourth corner model', but now interaction/fourth corner terms are optional as
#' well as are the main effects of trait covariates. The fourth corner model is implemented only using
#' variational approximation method.
#'
#'
#' \subsection{Distributions}{
#'   Mean and variance for distributions are defined as follows.
#'
#'   For count data \code{family="poisson"}: Expectation \eqn{E[Y_{ij}]=\mu_{ij}}, variance \eqn{V(\mu_{ij})=\mu_{ij}}.
#'   OR \code{family="negative.binomial"}: Expectation \eqn{E[Y_{ij}]=\mu_{ij}}, variance \eqn{V(\mu_{ij})=\mu_{ij}+\phi_j*\mu_{ij}^2}
#'   OR \code{family="ZIP"}: Expectation \eqn{E[Y_{ij}]=(1-p)\mu_{ij}}, variance \eqn{V(\mu_{ij})=\mu_{ij}(1-p)(1+\mu_{ij}p)}
#'
#'   For binary data \code{family="binomial"}: Expectation \eqn{E[Y_{ij}]=\mu_{ij}}, variance \eqn{V(\mu_{ij})=\mu_{ij}(1-\mu_{ij})}.
#'
#'   For biomass data \code{family="tweedie"}: Expectation \eqn{E[Y_{ij}]=\mu_{ij}}, variance \eqn{V(\mu_{ij})=\phi_j*\mu_{ij}^\nu}, where \eqn{\nu} is a power parameter of Tweedie distribution. See details Dunn and Smyth (2005).
#'
#'   For ordinal data \code{family="ordinal"}: See Hui et.al. (2016).
#' }
#'
#' - Method is sensitive for the choices of initial values of the latent variables. Therefore it is
#' recommendable to use multiple runs and pick up the one giving the highest log-likelihood value.
#' However, if this is computationally too demanding you can use the result of some classical
#' ordination method (for example nMDS) as initial values for latent variables.
#'
#' - If algorithm is very unstable in different runs, it is recommendable to use \code{constrOpt=TRUE}. The parameters are then estimated by constraining linear predictors in optimization. This has been done using function \code{constrOptim}.
#' Especially with binomial family, algorithm may often converge to a poor local maximum and constrained optimization may help in such cases.
#'
#'
#'
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
#'    \item{$beta0 }{ column specific intercepts}
#'    \item{$Xcoef }{ coefficients related to environmental covariates X}
#'    \item{$Tcoef }{ coefficients related to trait covariates TR}
#'    \item{$fourth }{ interaction terms}
#'    \item{$row.params }{ row specific intercepts}
#'    \item{$phi }{ dispersion parameters \eqn{\phi} for negative binomial or Tweedie family}
#'    \item{$inv.phi }{ dispersion parameters \eqn{1/\phi} for negative binomial}
#'    \item{$p }{ Probability of zero inflation for ZIP family}
#'    }}
#'  \item{Power }{ power parameter \eqn{\nu} for Tweedie family}
#'  \item{sd }{list of standard errors of parameters}
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui, Sara Taskinen
#' @references
#' Dunn, P. K. and Smyth, G. K. (2005).  Series evaluation of tweedie exponential dispersion model densities. Statistics and Computing, 15:267-280.
#'
#' Hui, F. K. C., Taskinen, S., Pledger, S., Foster, S. D., and Warton, D. I. (2015).  Model-Based Approaches to Unconstrained Ordination. Methods in Ecology and Evolution, 6:399-411.
#'
#' Hui, F. K. C., Warton, D., Ormerod, J., Haapaniemi, V., and Taskinen, S. (2016).  Variational Approximations for Generalized Linear Latent Variable Models. Journal of Computational and Graphical Statistics. Journal of Computational and Graphical Statistics, 26:35-43.
#'
#' Niku, J., Warton,  D. I., Hui, F. K. C., and Taskinen, S. (2017). Generalized Linear Latent Variable Models for Multivariate Count and Biomass Data in Ecology. Under revision.
#'
#' Warton, D. I., Guillaume Blanchet, F., O’Hara, R. B., Ovaskainen, O., Taskinen, s., Walker, S. C. and Hui, F. K. C. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology & Evolution, 30:766-779.
#'
#'@seealso \code{\link{ordplot.gllvm}}, \code{\link{summary.gllvm}}, \code{\link{residuals.gllvm}}, \code{\link{confint.gllvm}}, \code{\link{coefplot.gllvm}}.
#' @examples
#' \dontrun{
#'library(mvabund) ## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- as.matrix(antTraits$env[,4:5])
#'TR <- antTraits$traits[,4:5]
#'
#'## Example 1: Fit model with two latent variables
#'# Using Laplace approximation:
#'fitl0 <- gllvm(y=y,family = "negative.binomial", method="LA",seed=1)
#'ordplot.gllvm(fitl0)
#'# Using Variational analysis:
#'fitv0 <- gllvm(y=y,family = "negative.binomial", method="VA",seed=1)
#'ordplot.gllvm(fitv0)
#'
#'# Poisson family:
#'fit.p <- gllvm(y=y,family = "poisson", method="LA")
#'ordplot.gllvm(fit.p)
#'# Use poisson model as a starting parameters for ZIP-model:
#'fit.z <- gllvm(y=y,family = "ZIP", method="LA",start.param = fit.p)
#'ordplot.gllvm(fit.z)
#'
#'
#'## Example 2: Fit model with two latent variables and environmental covariates,
#'fitlX <- gllvm(y=y, X = X,family = "negative.binomial", method="LA")
#'ordplot.gllvm(fitlX)
#'# Use 5 initial runs and pick the best one
#'fitlX_5 <- gllvm(y=y, X = X,family = "negative.binomial", method="LA")
#'ordplot.gllvm(fitlX_5)
#'
#'
#'# Example 3: Fit fourth corner model with two latent variables
#'# This can be fitted using only Variational analysis:
#'fitF <- gllvm(y=y, X = X, TR=TR,family = "negative.binomial", method="VA")
#'ordplot.gllvm(fitF)
#'
#'# Example 4: Fit Tweedie model for coral data
#'data(tikus)
#'ycoral <- tikus$abund
#'# Exclude species which have observed at less than 5 sites
#'ycoral <- ycoral[(colSums(ycoral>0)>4)]
#'fit.twe <- gllvm(y=ycoral,family = "tweedie", method="LA")
#'ordplot.gllvm(fit.twe)
#'}
#' @export
#'

gllvm<-function(y, X = NULL, TR = NULL,data=NULL, n.X=0,n.TR=0, num.lv = 2, family, method = "VA", row.eff = FALSE, get.fourth=TRUE, get.trait=TRUE, offset=NULL, start.lvs = NULL, start.params=NULL, sd.errors = TRUE, Lambda.struc = "unstructured", diag.iter = 5, trace = FALSE, plot = FALSE,la.phi.upd="phi",la.link.bin="logit",n.init=1,fixed.power=TRUE,Power=1.5,constrOpt=FALSE,restrict=30, eps = 1e-5, seed = NULL, max.iter = 200, maxit = 1000){
  if(is.data.frame(y) && ncol(y)==3){
    id1<-!duplicated(y[,1])
    id2<-!duplicated(y[,2])
    n <- sum(id1); p <- sum(id2);
    Y<-reshape(y, direction = "wide", idvar = names(y)[1], timevar = names(y)[2],v.names = names(y)[3])[,-1]
    if(nrow(Y)!=n || ncol(Y)!=p) {stop("Data format of y not allowed. See instructions from help.");} else{y <- Y}

    if(n.X>0){
      X<-data[id1,1:n.X]
    }
    if(n.TR>0){
      TR<-data[id2,(n.X+1):(n.X+n.TR)]
    }
  }

  if((family == "tweedie" || family == "ZIP") && !is.null(TR)) {
    stop(paste(family," family cannot yet handle trait covariates."));
  }
  if(method == "LA" && family == "ordinal") {
    cat("Laplace's method cannot yet handle ordinal data, so VA method is used instead. \n")
    method = "VA"
  }
  if(method == "VA" && (family == "tweedie" || family == "ZIP")) {
    cat("VA method cannot handle", family," family, so LA method is used instead. \n")
    method = "LA"
  }
  if(method=="LA" && !is.null(TR)) {
    cat("Laplace's method cannot yet handle trait covariates, so VA method is used instead. \n")
    method = "VA"
  }
  if(!is.null(start.params)){
    if(class(start.params)!="gllvm") stop("Only object of class 'gllvm' can be given as a starting parameters.");
    if(!(family %in% c("poisson","negative.binomial","ZIP"))) stop("Starting parameters can be given only for count data.");
    }

  n<-nrow(y); p<-ncol(y);
  if (is.null(offset))
    O <- matrix(0, nrow = n, ncol = p)
  else if (NCOL(offset) == 1)
    O <- matrix(rep(offset), nrow = n, ncol = p)
  else O <- as.matrix(offset)

  n.i<-1;
  out<-list(y = y, X = X, TR = TR, num.lv = num.lv, method = method, family=family, row.eff = row.eff,n.init=n.init,get.fourth=get.fourth, get.trait=get.trait,sd=FALSE)
  if(family=="binomial"){
    if(method=="LA") out$link=la.link.bin
    if(method=="VA") out$link="probit"
  }
  out$offset <- offset;

  if(method == "VA") {

    fitg <- gllvm.VA(y, X = X, T = TR, family = family, num.lv = num.lv, max.iter = max.iter, eps = eps, row.eff = row.eff, get.fourth=get.fourth, get.trait=get.trait, Lambda.struc = Lambda.struc, trace = trace, plot = plot, sd.errors = sd.errors, start.lvs = start.lvs, offset=O, maxit = maxit, diag.iter = diag.iter, seed=seed,n.init = n.init,restrict=restrict,constrOpt=constrOpt,start.params=start.params)
    out$logL <- fitg$logLik
    if(num.lv>0) out$lvs <- fitg$lvs
    out$X = fitg$X; out$TR <- fitg$T
    out$params <- fitg$coef
    if(sd.errors){ out$sd <- fitg$sd}
    out$Lambda.struc <- fitg$Lambda.struc
    out$Lambda <- fitg$Lambda
  }
  if(method == "LA"){
    if(family=="tweedie"){
      fitg <- gllvm.tweedie(y, X = X, num.lv = num.lv, row.eff = row.eff, max.iter = max.iter, eps = eps, seed = seed, maxit = maxit, start.lvs = start.lvs, offset=O, info = sd.errors, trace = trace, plot = plot, n.init = n.init,fixed.power=fixed.power,Power=Power)
      out$logL <- fitg$logL
      if(num.lv > 0) out$lvs <- fitg$lvs
      out$Power <- fitg$Power
      out$X <- fitg$X;
      out$params <- fitg$params
      if(sd.errors) {
        out$sd <- fitg$sd
        out$lvs.cov <- fitg$lvs.cov
      }
    } else if(family=="ZIP"){
      fitg <- gllvm.ZIP(y, X = X, num.lv = num.lv, row.eff = row.eff, max.iter = max.iter, eps = eps, seed = seed, maxit = maxit, start.lvs = start.lvs, offset=O, info = sd.errors, trace = trace, plot = plot, n.init = n.init,start.params=start.params)
      out$logL <- fitg$logL
      if(num.lv > 0) out$lvs <- fitg$lvs
      out$X <- fitg$X;
      out$params <- fitg$params
      if(sd.errors) {
        out$sd <- fitg$sd
        out$lvs.cov <- fitg$lvs.cov
      }
      out$p0=fitg$p0
    } else {
      fitg <- gllvm.LA(y, X = X, num.lv = num.lv, family = family, row.eff = row.eff, max.iter = max.iter, eps = eps, seed = seed, maxit = maxit, start.lvs = start.lvs, offset=O, info = sd.errors, trace = trace, plot = plot,link=la.link.bin, n.init = n.init,phi.upd=la.phi.upd,restrict=restrict,constrOpt=constrOpt,start.params=start.params)
      out$logL <- fitg$logL
      if(num.lv > 0) out$lvs <- fitg$lvs
      out$X <- fitg$X;
      out$params <- list(theta = fitg$lambdas, beta0 = fitg$beta0, Xcoef = fitg$betas, row.params = fitg$row.params)
      if(sd.errors) {
        out$sd <- list(theta = fitg$se.lambdas, beta0 = fitg$se.beta0, Xcoef = fitg$se.Xcoefs, row.params = fitg$se.row.params,inv.phi = fitg$se.phis)
        out$lvs.cov <- fitg$lvs.cov
      }
      if(family == "negative.binomial") {
        out$params$phi <- fitg$phis
      }
    }
  }
  if(family=="negative.binomial") out$params$inv.phi=1/out$params$phi
  out$call <- match.call()
  class(out) <- "gllvm"
  return(out)
}




