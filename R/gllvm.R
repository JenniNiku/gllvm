#' @title Generalized Linear Latent Variable Models
#' @description Fits generalized linear latent variable model for multivariate data. The model can be fitted using Laplace approximation method or variational
#' approximation method.
#'
#' @param y (n x m) matrix of responses.
#' @param X matrix or data.frame of environmental covariates.
#' @param TR matrix or data.frame of trait covariates.
#' @param data data in long format, that is, matrix of responses, environmental and trait covariates and row index named as "id". When used, model needs to be defined using formula. This is alternative data input for y, X and TR.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted (for fixed-effects predictors).
#' @param family  distribution function for responses. Options are \code{"negative.binomial"} (with log link), \code{poisson(link = "log")}, \code{binomial(link = "probit")} (and also with \code{link = "logit"} when \code{method = "LA"} or \code{method = "EVA"}), zero-inflated poisson (\code{"ZIP"}), zero-inflated negative-binomial (\code{"ZINB"}), \code{gaussian(link = "identity")}, Tweedie (\code{"tweedie"}) (with log link, for \code{"LA"} and \code{"EVA"}-method), \code{"gamma"} (with log link), \code{"exponential"} (with log link), beta (\code{"beta"}) (with logit and probit link, for \code{"LA"} and  \code{"EVA"}-method) and \code{"ordinal"} (only with \code{"VA"}-method).
#' @param num.lv  number of latent variables, d, in gllvm model. Non-negative integer, less than number of response variables (m). Defaults to 2, if \code{num.lv.c=0} and \code{num.RR=0}, otherwise 0.
#' @param num.lv.c  number of latent variables, d, in gllvm model to constrain, with residual term. Non-negative integer, less than number of response (m) and equal to, or less than, the number of predictor variables (k). Defaults to 0. Requires specification of "lv.formula" in combination with "X" or "datayx". Can be used in combination with num.lv and fixed-effects, but not with traits.
#' @param num.RR number of latent variables, d, in gllvm model to constrain, without residual term (reduced rank regression). Cannot yet be combined with traits.
#' @param lv.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted (for latent variables).
#' @param lvCor (Under development, not to be used at the moment!) correlation structure for latent variables, defaults to \code{NULL} Correlation structure for latent variables can be defined via formula, eg. \code{~struc(1|groups)}, where option to 'struc' are \code{corAR1} (AR(1) covariance), \code{corExp} (exponentially decaying, see argument '\code{dist}') and \code{corCS} (compound symmetry). The grouping variable needs to be included either in 'X' or 'studyDesign'. Works at the moment only with unconstrained ordination without quadratic term.
#' @param studyDesign variables related to eg. sampling/study design, used for defining correlation structure of the latent variables and row effects.
#' @param method  model can be fitted using Laplace approximation method (\code{method = "LA"}) or variational approximation method (\code{method = "VA"}), or with extended variational approximation method (\code{method = "EVA"}) when VA is not applicable. If particular model has not been implemented using the selected method, model is fitted using the alternative method as a default. Defaults to \code{"VA"}.
#' @param row.eff  \code{FALSE}, \code{fixed}, \code{"random"} or formula to define the structure for the row parameters. Indicating whether row effects are included in the model as a fixed or as a random effects. Defaults to \code{FALSE} when row effects are not included. Structured random row effects can be defined via formula, eg. \code{~(1|groups)}, when unique row effects are set for each group, not for all rows, the grouping variable needs to be included in \code{X}. Correlation structure between random group effects/intercepts can also be set using \code{~struc(1|groups)}, where option to 'struc' are \code{corAR1} (AR(1) covariance), \code{corExp} (exponentially decaying, see argument '\code{dist}') and \code{corCS} (compound symmetry). Correlation structure can be set between or within groups, see argument '\code{corWithin}'.
#' @param corWithin logical. Vector of length equal to the number of row effects. For structured row effects with correlation, If \code{TRUE}, correlation is set between row effects of the observation units within group. Correlation and groups can be defined using \code{row.eff}. Defaults to \code{FALSE}, when correlation is set for row parameters between groups.
#' @param corWithinLV logical. For LVs with correlation, If \code{TRUE}, correlation is set between rows of the observation units within group. Defaults to \code{FALSE}, when correlation is set for rows between groups.
#' @param dist list of length equal to the number of row effects with correlation structure \code{corExp} that holds the matrix of coordinates or time points.
#' @param distLV matrix of coordinates or time points used for LV correlation structure \code{corExp}.
#' @param colMat matrix of similarity for the column effects, such as from a Phylogeny. Must be positive definite.
#' @param quadratic either \code{FALSE}(default), \code{TRUE}, or \code{LV}. If \code{FALSE} models species responses as a linear function of the latent variables. If \code{TRUE} models species responses as a quadratic function of the latent variables. If \code{LV} assumes species all have the same quadratic coefficient per latent variable.
#' @param randomB either \code{FALSE}(default), "LV", "P", or "single". Fits concurrent or constrained ordination (i.e. models with num.lv.c or num.RR) with random slopes for the predictors. "LV" assumes LV-specific variance parameters, "P" predictor specific, and "single" the same across LVs and predictors.
#' @param sd.errors  logical. If \code{TRUE} (default) standard errors for parameter estimates are calculated.
#' @param offset vector or matrix of offset terms.
#' @param Ntrials number of trials for binomial family.
#' @param link link function for binomial family if \code{method = "LA"} and beta family. Options are "logit" and "probit.
#' @param Power fixed power parameter in Tweedie model. Scalar from interval (1,2). Defaults to 1.1. If set to NULL it is estimated (note: experimental). 
#' @param seed a single seed value, defaults to \code{NULL}.
#' @param plot  logical. If \code{TRUE} ordination plots will be printed in each iteration step when \code{TMB = FALSE}. Defaults to \code{FALSE}.
#' @param zeta.struc structure for cut-offs in the ordinal model. Either "common", for the same cut-offs for all species, or "species" for species-specific cut-offs. For the latter, classes are arbitrary per species, each category per species needs to have at least one observations. Defaults to "species".
#' @param randomX  formula for species specific random effects of environmental variables in fourth corner model. Defaults to \code{NULL}, when random slopes are not included.
#' @param dependent.row logical. Whether or not random row effects are correlated (dependent) with the latent variables. Defaults to \code{FALSE} when correlation terms are not included.
#' @param beta0com logical. If \code{FALSE} column-specific intercepts are assumed. If \code{TRUE}, a common intercept is used which is allowed only for fourth corner models.
#' @param scale.X logical. If \code{TRUE}, covariates are scaled when fourth corner model is fitted.
#' @param return.terms logical. If \code{TRUE} 'terms' object is returned.
#' @param gradient.check logical. If \code{TRUE} gradients are checked for large values (>0.01) even if the optimization algorithm did converge.
#' @param disp.formula formula, or alternatively a vector of indices, for the grouping of dispersion parameters (e.g. in a negative-binomial distribution). Defaults to NULL so that all species have their own dispersion parameter. Is only allowed to include categorical variables. If a formula, data should be included as named rows in y.
#' @param setMap UNDER DEVELOPMENT, DO NOT USE! list of a set of parameters to be fixed
#' @param control A list with the following arguments controlling the optimization:
#' \itemize{
#'  \item{\emph{reltol}: }{ convergence criteria for log-likelihood, defaults to 1e-8.}
#'  \item{\emph{reltol.c}: }{ convergence criteria for equality constraints in ordination with predictors, defaults to 1e-8.}  
#'  \item{\emph{TMB}: }{ logical, if \code{TRUE} model will be fitted using Template Model Builder (TMB). TMB is always used if \code{method = "LA"}.  Defaults to \code{TRUE}.}
#'  \item{\emph{optimizer}: }{ if \code{TMB=TRUE}, log-likelihood can be optimized using \code{"\link{optim}"} (default) or \code{"\link{nlminb}"}. For ordination with predictors (num.RR>0 or num.lv.c>0) this can additionally be one of \code{alabama}(default), \code{nloptr(agl)} or \code{nloptr(sqp)}.}
#'  \item{\emph{max.iter}: }{ maximum number of iterations when \code{TMB = FALSE} or for \code{optimizer = "nlminb"} when \code{TMB = TRUE}, defaults to 4000.}
#'  \item{\emph{maxit}: }{ maximum number of iterations for optimizer, defaults to 4000.}
#'  \item{\emph{trace}: }{ logical, if \code{TRUE} in each iteration step information on current step will be printed. Defaults to \code{FALSE}. Only with \code{TMB = FALSE}.}
#'  \item{\emph{optim.method}: }{ optimization method to be used if optimizer is \code{"\link{optim}"},\code{"alabama"}, or  \code{"\link{nloptr}"}, but the latter two are only available in combination with at least two latent variables (i.e., num.RR+num.lv.c>1). Defaults to \code{"BFGS"}, but to \code{"L-BFGS-B"} for Tweedie family due the limited-memory use. For optimizer='alabama' this can be any \code{"\link{optim}"} method, or  \code{"\link{nlminb}"}. If optimizer = 'nloptr(agl)' this can be one of: "NLOPT_LD_CCSAQ", "NLOPT_LD_SLSQP", "NLOPT_LD_TNEWTON_PRECOND" (default), "NLOPT_LD_TNEWTON", "NLOPT_LD_MMA".}
#' }
#' @param control.va A list with the following arguments controlling the variational approximation method:
#' \itemize{
#'  \item{\emph{Lambda.struc}: }{ covariance structure of VA distributions for latent variables when \code{method = "VA"}, "unstructured" or "diagonal".}
#'  \item{\emph{Ab.struct}: }{ covariance structure of VA distributions for random slopes when \code{method = "VA"}, "unstructured" or "diagonal".}
#'  \item{\emph{Ar.struc}: }{ covariance structure of VA distributions for random row effects when \code{method = "VA"}, "unstructured" or "diagonal".}
#'  \item{\emph{sp.Ar.struc}: }{ covariance structure of VA distributions for random column effects when \code{method = "VA"}, "unstructured" or "diagonal".}
#'  \item{\emph{diag.iter}: }{ non-negative integer which can sometimes be used to speed up the updating of variational (covariance) parameters in VA method. Can sometimes improve the accuracy. If \code{TMB = TRUE} either 0 or 1. Defaults to 1.}
#'  \item{\emph{Ab.diag.iter}: }{ As above, but for variational covariance of random slopes.}
#'  \item{\emph{Lambda.start}: }{ starting values for variances in VA distributions for latent variables, random row effects and random slopes in variational approximation method. Defaults to 0.3.}
#'  \item{\emph{NN}: }{ Number of nearest neighbors for NN variational covariance. Defaults to ...}
#' }
#' @param control.start A list with the following arguments controlling the starting values:
#' \itemize{
#'   \item{\emph{starting.val}: }{ starting values can be generated by fitting model without latent variables, and applying factorial analysis to residuals to get starting values for latent variables and their coefficients (\code{starting.val = "res"}). Another options are to use zeros as a starting values (\code{starting.val = "zero"}) or initialize starting values for latent variables with (n x num.lv) matrix. Defaults to \code{"res"}, which is recommended.}
#'   \item{\emph{n.init}: }{ number of initial runs. Uses multiple runs and picks up the one giving highest log-likelihood value. Defaults to 1.}
#'   \item{\emph{n.init.max}: }{ maximum number of refits try try for n.init without improvement, defaults to 10.}
#'   \item{\emph{start.fit}: }{ object of class 'gllvm' which can be given as starting parameters for count data (poisson, NB, or ZIP).}
#'   \item{\emph{start.lvs}: }{ initialize starting values for latent variables with (n x num.lv) matrix. Defaults to \code{NULL}.}
#'   \item{\emph{jitter.var}: }{ jitter variance for starting values of latent variables. Defaults to 0, meaning no jittering.}
#'   \item{\emph{randomX.start}: }{ starting value method for the random slopes. Options are \code{"zero"} and \code{"res"}. Defaults to \code{"zero"}.}
#'   \item{\emph{start.struc}: }{ starting value method for the quadratic term. Options are \code{"LV"} (default) and \code{"all"}.}
#'   \item{\emph{quad.start}: }{ starting values for quadratic coefficients. Defaults to 0.01.}
#'   \item{\emph{MaternKappa}: }{ Starting value for smoothness parameter kappa of Matern covariance function. Defaults to 3/2.}
#'   \item{\emph{scalmax}: }{ Sets starting value for the scale parameter for the coordinates. Defaults to 10, when the starting value for scale parameter scales the distances of coordinates between 0-10.}
#'   \item{\emph{rangeP}: }{ Sets starting value for the range parameter for the correlation structure.}
#' }
#' @param ... Not used.
#'
#' @details
#' Fits generalized linear latent variable models as in Hui et al. (2015 and 2017) and Niku et al. (2017).
#' Method can be used with two types of latent variable models depending on covariates. If only
#' site related environmental covariates are used, the expectation of response \eqn{Y_{ij}} is determined by
#'
#' \deqn{g(\mu_{ij}) = \eta_{ij} = \alpha_i + \beta_{0j} + x_i'\beta_j + u_i'\theta_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{u_i} are \eqn{d}-variate latent variables (\eqn{d}<<\eqn{m}), \eqn{\alpha_i} is an optional row effect
#' at site \eqn{i}, and it can be fixed or random effect (also other structures are possible, see below), \eqn{\beta_{0j}} is an intercept term for species \eqn{j}, \eqn{\beta_j} and \eqn{\theta_j} are column
#' specific coefficients related to covariates and the latent variables, respectively.
#'
#' \subsection{Quadratic model}{
#'Alternatively, a more complex version of the model can be fitted with \code{quadratic = TRUE}, where species are modeled as a quadratic function of the latent variables:
#' \deqn{g(\mu_{ij}) = \eta_{ij} = \alpha_i + \beta_{0j} + x_i'\beta_j + u_i'\theta_j - u_i' D_j u_i}.
#'Here, D_j is a diagonal matrix of positive only quadratic coefficients, so that the model generates concave shapes only. This implementation follows
#'the ecological theoretical model where species are generally recognized to exhibit non-linear response curves.
#'For a model with quadratic responses, quadratic coefficients are assumed to be the same for all species: \deqn{D_j = D}. This model requires less information
#'per species and can be expected to be more applicable to most datasets. The quadratic coefficients D can be used to calculate the length of 
#'ecological gradients.
#'For quadratic responses, it can be useful to provide the latent variables estimated with a GLLVM with linear responses, or estimated with (Detrended) Correspondence Analysis.
#'The latent variables can then be passed to the \code{start.lvs} argument inside the \code{control.start} list, which in many cases gives good results. 
#'}
#'
#' \subsection{Ordination with predictors}{
#'For GLLVMs with both linear and quadratic response model, a series of predictors \eqn{x_{lv}} can be included to explain the latent variables:
#'
#'\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + x_i'\beta_j + (B' x_{lv,i} + \epsilon_i)' \gamma_j - (B' x_{lv,i} + \epsilon_i)' D_j (B' x_{lv,i} + \epsilon_i) ,}
#'where \eqn{z_i = B' x_{lv,i} + \epsilon_i} are latent variables informed by the predictors, but not constrained compared to unconstrained ordination as in methods such as CCA or RDA.
#' Omitting the predictors results in an unconstrained ordination, and omitting \eqn{\epsilon_i} in the usual constrained ordination, which can also be fitted.
#'}
#'
#' \subsection{Fourth corner model}{
#' An alternative model is the fourth corner model (Brown et al., 2014, Warton et al., 2015) which will be fitted if also trait covariates
#' are included. The expectation of response \eqn{Y_{ij}} is
#'
#' \deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + x_i'(\beta_x + b_j) + TR_j'\beta_t + vec(B)*kronecker(TR_j,X_i) + u_i'\theta_j - u_i'D_ju_i}
#'
#' where g(.), \eqn{u_i}, \eqn{\beta_{0j}} and \eqn{\theta_j} are defined as above. Vectors \eqn{\beta_x} and \eqn{\beta_t} are the main effects
#' or coefficients related to environmental and trait covariates, respectively, matrix \eqn{B} includes interaction terms. Vectors \eqn{b_j} are 
#' optional species-specific random slopes for environmental covariates.
#' The interaction/fourth corner terms are optional as well as are the main effects of trait covariates.
#'}
#'
#' \subsection{Structured row effects}{
#' In addition to the site-specific random effects, \eqn{\alpha_i}, it is also possible to set arbitrary structure/design for the row effects. 
#' That is, assume that observations / rows \eqn{i=1,...,n} in the data matrix are from groups \eqn{t=1,...,T}, so that each row \eqn{i} belongs to one of the groups, denote \eqn{G(i) \in \{1,...,T\}}. Each group \eqn{t} has a number of observations \eqn{n_t}, so that \eqn{\sum_{t=1}^{T} n_t =n}.
#' Now we can set random intercept for each group \eqn{t}, (see argument '\code{row.eff}'):
#' 
#'  \deqn{g(\mu_{ij}) = \eta_{ij} = \alpha_{G(i)} + \beta_{0j} + x_i'\beta_j + u_i'\theta_j,}
#'  
#'  There is also a possibility to set correlation structure for the random intercepts between groups, so that \eqn{(\alpha_{1},...,\alpha_{T})^\top \sim N(0, \Sigma_r)}. That might be the case, for example, when the groups are spatially or temporally dependent.
#'  Another option is to set row specific random intercepts \eqn{\alpha_i}, but to set the correlation structure for the observations within groups, (see argument '\code{corWithin}'). That is, we can set \eqn{corr(\alpha_{i},\alpha_{i'}) = C(i,i') \neq 0} according to some correlation function \eqn{C}, when \eqn{G(i)=G(i')}.
#'  This model is restricted to the case, where each group has equal number of observations (rows), that is \eqn{n_t=n_{t'}} for all \eqn{t,t' \in \{1,...,T\}}.
#'  
#'  The correlation structures available in the package are 
#'\itemize{
#'   \item{\code{corAR1}} { autoregressive process of order 1.}
#'   \item{\code{corExp}} { exponentially decaying, see argument '\code{dist}'.}
#'   \item{\code{corCS}} { compound symmetry.}
#' }  
#' }
#'
#'\subsection{Starting values}{
#' The method is sensitive for the choices of initial values of the latent variables. Therefore it is
#' recommendable to use multiple runs and pick up the one giving the highest log-likelihood value (see argument '\code{n.init}').
#' However, sometimes this is computationally too demanding, and default option
#' \code{starting.val = "res"} is recommended. For more details on different starting value methods, see Niku et al., (2018).
#'}
#' Models are implemented using TMB (Kristensen et al., 2015) applied to variational approximation (Hui et al., 2017), extended variational approximation (Korhonen et al., 2021) and Laplace approximation (Niku et al., 2017).
#'
#' With ordinal family response classes must start from 0 or 1.
#'
#' \subsection{Distributions}{
#'
#'   Mean and variance for distributions are defined as follows.
#'\itemize{
#'   \item{For count data \code{family = poisson()}:} {Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}}, or}
#'   \item{ \code{family = "negative.binomial"}:}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}+\mu_{ij}^2\phi_j}, or}
#'   \item{ \code{family = "ZIP"}:}{ Expectation \eqn{E[Y_{ij}] = (1-p)\mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}(1-p_j)(1+\mu_{ij}p)}.}
#'   \item{ \code{family = "ZINB"}:}{ Expectation \eqn{E[Y_{ij}] = (1-p)\mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}(1-p_j)(1+\mu_{ij}(\phi_j+p_j))}.}
#'   \item{For binary data \code{family = binomial()}:}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}(1-\mu_{ij})}.}
#'   
#'   \item{For percent cover data \eqn{0 < Y_{ij} < 1} \code{family = "beta"}:}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}(1-\mu_{ij})//(1+\phi_j)}.}
#'
#'   \item{For positive continuous data \code{family = "gamma"}:}{Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}^2/\phi_j}, where \eqn{\phi_j} is species specific shape parameter.}
#'   
#'   \item{For non-negative  continuous data \code{family = "exponential"}:}{Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \mu_{ij}^2}.}
#'   
#'   \item{For non-negative continuous or biomass data \code{family = "tweedie"}}{ Expectation \eqn{E[Y_{ij}] = \mu_{ij}}, variance \eqn{V(\mu_{ij}) = \phi_j*\mu_{ij}^\nu}, where \eqn{\nu} is a power parameter of Tweedie distribution. See details Dunn and Smyth (2005).}
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
#'  \item{call }{ function call}
#'  \item{y}{ (n x m) matrix of responses.}
#'  \item{X}{ matrix or data.frame of environmental covariates.}
#'  \item{lv.X}{ matrix or data.frame of environmental covariates for latent variables.}
#'  \item{TR}{ Trait matrix}
#'  \item{formula}{ Formula for predictors}
#'  \item{lv.formula}{ Formula of latent variables in constrained and concurrent ordination}
#'  \item{randomX }{ Formula for species specific random effects in fourth corner model}
#'  \item{randomB }{ Boolean flag for random slopes in constrained and concurrent ordination}
#'  \item{num.lv}{ Number of unconstrained latent variables}
#'  \item{num.lv.c}{ Number of latent variables in concurrent ordination}
#'  \item{num.RR}{ Number of latent variables in constrained ordination}
#'  \item{Ntrials}{ Number of trials in a binomial model}
#'  \item{method}{ Method used for integration}
#'  \item{family}{ Response distribution}
#'  \item{row.eff}{ Type of row effect used}
#'  \item{n.init}{ Number of model runs for best fit}
#'  \item{disp.group}{ Groups for dispersion parameters}
#'  \item{sd }{ List of standard errors}
#'  \item{lvs }{ Latent variables}
#'  \item{params}{ List of parameters
#'  \itemize{
#'    \item{theta }{ latent variables' loadings relative to the diagonal entries of loading matrix}
#'    \item{sigma.lv }{ diagonal entries of latent variables' loading matrix}
#'    \item{LvXcoef }{ Predictor coefficients (or predictions for random slopes) related to latent variables, i.e. canonical coefficients}
#'    \item{beta0 }{ column specific intercepts}
#'    \item{Xcoef }{ coefficients related to environmental covariates X}
#'    \item{B }{ coefficients in fourth corner model}
#'    \item{betar}{ column random effects}
#'    \item{row.params }{ row-specific intercepts}
#'    \item{sigma }{ scale parameters for row-specific random effects}
#'    \item{phi }{ dispersion parameters \eqn{\phi} for negative binomial or Tweedie family, probability of zero inflation for ZIP family, standard deviation for gaussian family or shape parameter for gamma/beta family}
#'    \item{inv.phi }{ dispersion parameters \eqn{1/\phi} for negative binomial}
#'    }}
#'  \item{Power }{ power parameter \eqn{\nu} for Tweedie family}
#'  \item{sd }{ list of standard errors of parameters}
#'  \item{prediction.errors }{ list of prediction covariances for latent variables and variances for random row effects when method \code{"LA"} is used}
#'  \item{A, Ar, Ab_lv, spArs}{ covariance matrices for variational densities of latent variables, random row effects, random slopes, and colum effects respectively}
#'  \item{seed}{ Seed used for calculating starting values}
#'  \item{TMBfn}{ TMB objective and derivative functions}
#'  \item{logL }{ log likelihood}
#'  \item{convergence }{ convergence code of optimizer}
#'  \item{quadratic }{ flag for quadratic model}
#'  \item{Hess }{ List holding matrices of second derivatives}
#'  \item{beta0com }{ Flag for common intercept in fourth corner models}
#'  \item{cstruc }{ Correlation structure for row effects}
#'  \item{cstruclv }{ Correlation structure for LVs}
#'  \item{dist }{ Matrix of coordinates or time points used for row effects}
#'  \item{distLV }{ Matrix of coordinates or time points used for LVs}
#'  \item{col.eff }{ list of components for column random effects}
#'  \itemize{
#'  \item{col.eff }{flag indicating if column random effects are included}
#'  \item{colL }{ lower cholesky factor for covariance matrix of column effects}
#'  \item{sdrp }{ sparse design matrix}
#'  \item{nsp }{ vector of levels per random effect}
#'  
#'  }
#'  \item{terms }{ Terms object for main predictors}
#'  \item{start }{ starting values for model}
#'  \item{optim.method }{ Optimization method when using 'optim', 'alabama', or 'nloptr'}
#'  
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Wesley Brooks, Riki Herliansyah, Francis K.C. Hui, Pekka Korhonen, Sara Taskinen, Bert van der Veen, David I. Warton
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
#' Korhonen, P., Hui, F. K. C., Niku, J., and Taskinen, S. (2021). Fast, universal estimation of latent variable models using extended variational approximations. Stat Comput 33, 26 (2023).
#'
#' Niku, J., Warton,  D. I., Hui, F. K. C., and Taskinen, S. (2017). Generalized linear latent variable models for multivariate count and biomass data in ecology. Journal of Agricultural, Biological, and Environmental Statistics, 22:498-522.
#'
#' Niku, J., Brooks, W., Herliansyah, R., Hui, F. K. C., Taskinen, S., and Warton,  D. I. (2018). Efficient estimation of generalized linear latent variable models. PLoS One, 14(5):1-20.
#'
#' Warton, D. I., Guillaume Blanchet, F., O'Hara, R. B., Ovaskainen, O., Taskinen, S., Walker, S. C. and Hui, F. K. C. (2015). So many variables: Joint modeling in community ecology. Trends in Ecology & Evolution, 30:766-779.
#'
#'@seealso  \code{\link{coefplot.gllvm}}, \code{\link{confint.gllvm}}, \code{\link{ordiplot.gllvm}}, \code{\link{plot.gllvm}}, \code{\link{summary.gllvm}}.
#' @examples
#'# Extract subset of the microbial data to be used as an example
#'data(microbialdata)
#'X <- microbialdata$Xenv
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, X, formula = ~ pH + Phosp, family = poisson())
#'fit$logL
#'ordiplot(fit)
#'coefplot(fit)
#'
#' \donttest{
#'# Inclusion of structured random row effect
#'sDesign<-data.frame(Site = microbialdata$Xenv$Site)
#'fit <- gllvm(y, X, formula = ~ pH + Phosp, family = poisson(), 
#'             studyDesign=sDesign, row.eff=~(1|Site))
#'
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
#'## Example 1: Fit model with two unconstrained latent variables
#'# Using variational approximation:
#'fitv0 <- gllvm(y, family = "negative.binomial", method = "VA")
#'ordiplot(fitv0)
#'plot(fitv0, mfrow = c(2,2))
#'summary(fitv0)
#'confint(fitv0)
#'
#'## Example 1a: Fit concurrent ordination model with two latent variables and with 
#'# quadratic response model
#'# We scale and centre the  predictors to improve convergence
#'fity1 <- gllvm(y, X = scale(X), family = "negative.binomial", 
#'               num.lv.c=2, method="VA")
#'ordiplot(fity1, biplot = TRUE)
#'
#'#'## Example 1b: Fit constrained ordination model with two latent variables and with 
#'# random canonical coefficients
#'fity2 <- gllvm(y, X = scale(X), family = "negative.binomial", 
#'               num.RR=2, randomB="LV", method="VA")
#'               
#'# Using Laplace approximation: (this line may take about 30 sec to run)
#'fitl0 <- gllvm(y, family = "negative.binomial", method = "LA")
#'ordiplot(fitl0)
#'
#'# Poisson family:
#'fit.p <- gllvm(y, family = poisson(), method = "LA")
#'ordiplot(fit.p)
#'# Use poisson model as a starting parameters for ZIP-model, this line 
#'# may take few minutes to run
#'fit.z <- gllvm(y, family = "ZIP", method = "LA", 
#'               control.start = list(start.fit = fit.p))
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
#'@importFrom mvabund manyglm
#'@importFrom graphics abline axis par plot segments text points boxplot panel.smooth lines polygon arrows 
#'@importFrom grDevices rainbow hcl
#'@importFrom stats dnorm pnorm qnorm rnorm dbinom pbinom rbinom pnbinom rnbinom pbeta rbeta pexp rexp pgamma rgamma ppois rpois runif pchisq qchisq qqnorm lm AIC binomial constrOptim factanal glm model.extract model.frame model.matrix model.response nlminb optim optimHess reshape residuals terms BIC qqline sd formula ppoints quantile gaussian cov princomp as.formula residuals.lm coef printCoefmat nobs predict
#'@importFrom Matrix bdiag chol2inv diag
#'@importFrom MASS ginv polr mvrnorm
#'@importFrom mgcv gam predict.gam
#'@importFrom nloptr nloptr
#'@importFrom alabama auglag
#'@importFrom utils combn tail relist
#'@importFrom methods cbind2 rbind2
#'



gllvm <- function(y = NULL, X = NULL, TR = NULL, data = NULL, formula = NULL, family,
                  num.lv = NULL, num.lv.c = 0, num.RR = 0, lv.formula = NULL,
                  lvCor = NULL, studyDesign=NULL,dist = list(matrix(0)), distLV = matrix(0), colMat = NULL, corWithin = FALSE, corWithinLV = FALSE,
                  quadratic = FALSE, row.eff = FALSE, sd.errors = TRUE, offset = NULL, method = "VA", randomB = FALSE,
                  randomX = NULL, dependent.row = FALSE, beta0com = FALSE, zeta.struc = "species",
                  plot = FALSE, link = "probit", Ntrials = 1,
                  Power = 1.1, seed = NULL, scale.X = TRUE, return.terms = TRUE, 
                  gradient.check = FALSE, disp.formula = NULL,
                  control = list(reltol = 1e-10, reltol.c = 1e-8, TMB = TRUE, optimizer = ifelse((num.RR+num.lv.c)==0 | randomB!=FALSE,"optim","alabama"), max.iter = 6000, maxit = 6000, trace = FALSE, optim.method = NULL), 
                  control.va = list(Lambda.struc = "unstructured", Ab.struct = "unstructured", Ar.struc="diagonal", sp.Ar.struc  = "diagonal", diag.iter = 1, Ab.diag.iter=0, Lambda.start = c(0.3, 0.3, 0.3), NN = 3),
                  control.start = list(starting.val = "res", n.init = 1, n.init.max = 10, jitter.var = 0, start.fit = NULL, start.lvs = NULL, randomX.start = "zero", quad.start=0.01, start.struc = "LV", scalmax = 10, MaternKappa=1.5, rangeP=NULL), setMap=NULL, ...
                  ) {
  # Dthreshold=0,
  if(!method%in%c("LA","VA","EVA"))stop("Selected method is not supported.")
  if(method=="EVA" && !isFALSE(quadratic))stop("The EVA method is not available for a quadratic GLLVM.")
  if(!is.null(lvCor) & corWithin) warning("'lvCor' with 'corWithin = TRUE' is under development, so all properties may not work properly.")
    #change default behavior of num.lv.
    #if num.lv.c>0, num.lv defaults to 0 if it is 0. Otherwise, it defaults to 2
  if(randomB!=FALSE&&quadratic!=FALSE&&(num.lv.c+num.RR)>0&&method=="LA"){
    stop("Model with quadratic responses and random slopes not allowed with method 'LA'")
  }
  if((num.RR+num.lv.c)==0){
    randomB <- FALSE
  }

  if(!randomB%in%c(FALSE,"single","P","LV")){
    stop("RandomB should be one of FALSE, 'single', 'P', or 'LV'")
  }
  
  if(is.null(num.lv)&num.lv.c==0&num.RR==0){
    num.lv <- 2
  }else if(is.null(num.lv)){num.lv<-0}

    constrOpt <- FALSE
    restrict <- 30
    term <- NULL
    term2 <- NULL
    datayx <- NULL
    if(is.null(X)|!is.null(X)&(num.lv.c+num.RR)==0)lv.X <- NULL
    pp.pars <- list(...)
    
    if (inherits(family,"family")) {
      link <- family$link
      family <- family$family
    }
    
    fill_control = function(x){
      if (!("reltol" %in% names(x))) 
        x$reltol = 1e-8
      if (!("reltol.c" %in% names(x))) 
        x$reltol.c = 1e-8
      if (!("TMB" %in% names(x))) 
        x$TMB = TRUE
      if (!("optimizer" %in% names(x))) 
        x$optimizer = ifelse((num.RR+num.lv.c)<=1 | randomB!=FALSE,"optim","alabama")
      if((num.lv.c+num.RR)>1 && family =="tweedie") x$optimizer = "alabama"
      if (!("optim.method" %in% names(x)) | is.null(x$optim.method)) {
        if(family=="tweedie") x$optim.method = "L-BFGS-B" else x$optim.method = "BFGS"
        if((num.RR+num.lv.c)>1 && randomB == FALSE && family!="tweedie" && x$optimizer%in%c("nloptr(agl)","nloptr(sqp)")) x$optim.method = "NLOPT_LD_TNEWTON_PRECOND"
      }
      if (!("max.iter" %in% names(x))) 
        x$max.iter = 6000
      if (!("maxit" %in% names(x))) 
        x$maxit = 6000
      if (!("trace" %in% names(x))) 
        x$trace = FALSE
      x
    }
    
    fill_control.va = function(x){
      if (!("Lambda.struc" %in% names(x)) & is.null(lvCor)) 
        x$Lambda.struc = "unstructured"
      if (!("Lambda.struc" %in% names(x)) & !is.null(lvCor)) 
        x$Lambda.struc = "diagonal"
      if (!("Ab.struct" %in% names(x))) 
        x$Ab.struct = "unstructured"
      if (!("Ar.struc" %in% names(x))) 
        x$Ar.struc = "diagonal"
      if (!("sp.Ar.struc" %in% names(x))) 
        x$sp.Ar.struc = "diagonal"
      if (!("diag.iter" %in% names(x))) 
        x$diag.iter = 5
      if (!("Ab.diag.iter" %in% names(x))) 
        x$Ab.diag.iter = 0
      if (!("Lambda.start" %in% names(x))) 
        x$Lambda.start = c(0.3, 0.3, 0.3)
      if (!("NN" %in% names(x)))
        x$NN = 3   
      x
    }
    fill_control.start = function(x){
      if (!("starting.val" %in% names(x))) 
        x$starting.val = "res"
      if (!("n.init" %in% names(x))) 
        x$n.init = 1
      if (!("n.init.max" %in% names(x))) 
        x$n.init.max = 10
      if (!("jitter.var" %in% names(x))) 
        x$jitter.var = 0
      if (!("start.fit" %in% names(x))) 
        x$start.fit = NULL
      if (!("start.lvs" %in% names(x))) 
        x$start.lvs = NULL
      if (!("randomX.start" %in% names(x))) 
        x$randomX.start = "res"
      if (!("quad.start" %in% names(x))) 
        x$quad.start = 0.01
      if (!("start.struc" %in% names(x))) 
        x$start.struc = "LV"
      if (!("scalmax" %in% names(x))) 
        x$scalmax = 10
      if (!("rangeP" %in% names(x))) 
        x$rangeP = NULL
      if (!("MaternKappa" %in% names(x))) 
        x$MaternKappa = 1.5
      x
    }

    if (inherits(family,"family")) {
      link <- family$link
      family <- family$family
    }  
    
    control <- fill_control(c(pp.pars, control))
    control.va <- fill_control.va(c(pp.pars, control.va))
    control.start <- fill_control.start(c(pp.pars, control.start))
    
  #some checks for optimizer

  # Cannot use nloptr or alabama with randomB
  if(randomB!=FALSE && control$optimizer %in% c("alabama","nloptr(sqp)","nloptr(agl)")){
    warning("Random slope models should use 'nlminb' or 'optim' as optimizer. Changing to 'optim'.")
    control$optimizer <- 'optim'
    if(family != "tweedie") {control$optim.method <- 'BFGS'}else{control$optim.method <- 'L-BFGS-B'}
    
  }
  
  # Define valid optimization routines
    if(control$optimizer=="optim" && !control$optim.method%in%c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")){
      stop("Invalid optim.method '", control$optim.method,"'")
    }
  if(!control$optimizer%in%c("optim","nlminb","alabama","nloptr(sqp)","nloptr(agl)")){
    stop("Optimizer must be one of 'optim', 'nlminb', 'alabama', 'nloptr(sqp)' or 'nloptr(agl)'.")
  }else if(control$optimizer%in%c("nloptr(sqp)","nloptr(agl)")){
    # Change to NLOPT algorithm names
    if(control$optimizer=="nloptr(sqp)"){control$optimizer <- "NLOPT_LD_SLSQP"}else if(control$optimizer=="nloptr(agl)"){control$optimizer <- "NLOPT_LD_AUGLAG_EQ"}
  }
  # cannot use alabama or nloptr without num.lv.c or num.RR for now
    if((num.lv.c+num.RR)<=1 && control$optimizer %in% c("alabama","nloptr(sqp)","nloptr(agl)")){
      warning("Selected optimizer not available for this model. Using optim instead.")
      control$optimizer <- "optim"
      if(family!="tweedie")control$optim.metod <- "BFGS"
      if(family=="tweedie")control$optim.method <- "L-BFGS-B"
    }
    
  if((num.RR+num.lv.c)>1 && control$optimizer%in%c("optim","nlminb") && randomB == FALSE){
    warning("Cannot fit ordination with predictors using 'optim' or 'nlminb', using 'nloptr(agl)' instead.")
    control$optimizer <- "nloptr(agl)"
  }
  if(family=="tweedie" && (num.lv.c+num.RR)>1 && control$optimizer != "alabama"){
    warning("Due to memory issues only optimizer 'alabama' with 'optim.method='L-BFGS-B' can be used with Tweedie.")
    control$optimizer <- "alabama"
    control$optim.method <- "L-BFGS-B"
  }
    
    # Check if local solver for nloptr augmented lagranian algorithm is one of the defined options
    if((num.RR+num.lv.c)>1 && randomB == FALSE && control$optimizer == "nloptr(agl)"){
      if(!control$optim.method%in%c("NLOPT_LD_CCSAQ", "NLOPT_LD_SLSQP", "NLOPT_LD_TNEWTON_PRECOND", "NLOPT_LD_TNEWTON", "NLOPT_LD_MMA"))control$optim.method <- "NLOPT_LD_TNEWTON_PRECOND"
    }
    
    if(randomB!=FALSE&!control$TMB){
      stop("Random slopes in ordination only allows with TMB = TRUE.")
    }
    if(family == "binomial" && max(Ntrials) == 1 && !is.null(y) && max(y, na.rm=TRUE)>1){
    stop("Using the binomial distribution requires setting the `Ntrials` argument.")
    }
    if(family == "binomial" && method == "EVA" && max(Ntrials) >1){
      stop("Binomial distribution not yet supported with the EVA method.")
    }
    # if(num.RR>0&quadratic>0&(num.lv+num.lv.c)==0){
    #   control.start$start.struc <- "all"
    # }
    reltol = control$reltol; reltol.c = control$reltol.c; TMB = control$TMB; optimizer = control$optimizer; max.iter = control$max.iter; maxit = control$maxit; trace = control$trace; optim.method = control$optim.method
    Lambda.struc = control.va$Lambda.struc; Ab.struct = control.va$Ab.struct; Ar.struc = control.va$Ar.struc; sp.Ar.struc = control.va$sp.Ar.struc; diag.iter = control.va$diag.iter; Ab.diag.iter=control.va$Ab.diag.iter; Lambda.start = control.va$Lambda.start; NN = control.va$NN;
    starting.val = control.start$starting.val; n.init = control.start$n.init; n.init.max = control.start$n.init.max; jitter.var = control.start$jitter.var; start.fit = control.start$start.fit; start.lvs = control.start$start.lvs; randomX.start = control.start$randomX.start
    start.struc = control.start$start.struc;quad.start=control.start$quad.start; scalmax=control.start$scalmax; rangeP=control.start$rangeP; MaternKappa=control.start$MaternKappa
    
    if(!is.null(TR)&num.lv.c>0|!is.null(TR)&num.RR>0){
      stop("Cannot fit model with traits and reduced rank predictors. \n")
    }
    
    if(!is.null(start.fit)){
    if(start.fit$num.lv.c!=num.lv.c&start.fit$num.lv!=start.fit$num.lv){
      stop("Cannot use gllvm with different num.lv and num.lv.c as starting values.")
    }
    if(!inherits(start.fit,"gllvm.quadratic")&quadratic!=FALSE){
      stop("Cannot use gllvm with linear responses as starting fit for gllvm with quadratic responses.")
    }
    }
    
    if((num.lv.c+num.RR)>0&method=="VA"&TMB==FALSE){
      warning("Concurrent and constrained ordination only implemented with TMB. Setting TMB to TRUE.\n")
      control$TMB <- TRUE
    }

    # if (inherits(family,"family")) {
    #   link <- family$link
    #   family <- family$family
    # }  

    if(is.null(optim.method) && optimizer == "optim") optim.method <- ifelse(family == "tweedie", "L-BFGS-B", "BFGS")

    if(!is.null(X)){
      if(!is.matrix(X) && !is.data.frame(X) ) 
        stop("X must be a matrix or data.frame.")
      if(any(is.na(X)) )
        stop("NAs are not allowed in 'X'.")
    }
    if(!is.null(TR)){
      if(!is.matrix(TR) && !is.data.frame(TR) )
        stop("TR must be a matrix or data.frame.")
      if(any(is.na(TR)) )
        stop("NAs are not allowed in 'TR'.")
    }
    #is.null(X)&is.null(data)&num.lv.c>0|
    if((num.RR+num.lv.c)>0&is.null(X)&is.null(data)){
      stop("Cannot constrain latent variables without predictors. Please provide X, or set num.lv.c=0 or num.RR=0. \n")
    }
    
    if(!is.null(disp.formula)&!TMB){
      stop("Grouped dispersion parameters not allowed with TMB = FALSE.")
    }
    
    if(!is.null(disp.formula)){
      if(!is.vector(disp.formula)){
      if(!is.null(y)){
        if(all(all.vars(disp.formula)%in%row.names(y))){
          disp.group <- as.factor(do.call(paste,list(c(t(y)[,all.vars(disp.formula)]))))
          y <- y[!row.names(y)%in%all.vars(disp.formula),]
          levels(disp.group) <- 1:length(levels(disp.group))
          #check if row numbers are still sequential if so renumber
          if(all(diff(as.numeric(row.names(y)))==1)){
            row.names(y) <- 1:nrow(y)
          }
        }else{
          stop("Grouping variable for dispersion needs to be included as named rows in 'Y'")
        }
      }
      }
    }
    if((num.lv.c+num.RR)>0&!is.null(formula)&is.null(lv.formula)){
      stop("'lv.formula' should be provided when 'formula' is used with concurrent or constrained ordination.")
    }
    
    # separate species random effects
    if(length(X)>0 & length(studyDesign)>0){
      X.col.eff <- cbind(X,studyDesign)
    }else if(length(X)>0){
      X.col.eff <- X
    }else if(length(studyDesign)>0){
      X.col.eff <- studyDesign
    }else{
      X.col.eff <- NULL
    }
    col.eff <- FALSE;col.eff.formula = ~0
    if(anyBars(formula)){
      col.eff <- "random"
      col.eff.formula <- reformulate(sprintf("(%s)", sapply(findbars1(formula), deparse1)))# take out fixed effects
      formula <- nobars1_(formula) # take out random effects
      if(is.null(formula) && is.null(lv.formula)){
        X <- NULL
      }
    }
    if(anyBars(col.eff.formula) && is.null(X.col.eff)){
      stop("Covariates for species random effects must be provided.")
    }
    
    if (!is.null(y)) {
      y <- as.matrix(y)
      if (is.null(X) && is.null(TR)) {
        datayx <- list(y)
        m1 <- model.frame(y ~ NULL, data = datayx)
        term <- terms(m1)
      } else if (is.null(TR)) {
        if (is.null(formula)&is.null(lv.formula)&(num.lv.c+num.RR)==0) {
          ff <- formula(paste("~", "0", paste("+", colnames(X), collapse = "")))
          if (is.data.frame(X)) {
            datayx <- list(y = y, X = model.matrix(ff, X))
          } else {
            datayx <- list(y = y, X = X)
          }
          m1 <- model.frame(y ~ X, data = datayx)
          term <- terms(m1)
        } else if(is.null(formula)&is.null(lv.formula)&(num.lv.c+num.RR)>0){

          if(inherits(row.eff,"formula")){
            if(any(colnames(X)==all.vars(row.eff))){
              lv.formula <- formula(paste("~", paste(colnames(X[,-which(colnames(X)==all.vars(row.eff))]), collapse = "+")))              
            }else{
              lv.formula <- formula(paste("~", paste(colnames(X), collapse = "+")))
            }

            if (is.data.frame(X)) {
              datayx <- list(X = model.matrix(lv.formula, X)[,-1,drop=F])
            } else {
              datayx <- list(X = X)
            }
            lv.X <- as.matrix(model.frame(~ X, data = datayx))
            if(!is.null(row.names(lv.X)))row.names(lv.X)<-row.names(X)
            if(any(colnames(X)==all.vars(row.eff))){
              X <-  X[,all.vars(row.eff),drop=F]
            }else{
             X <- NULL 
            }
          }else{
            lv.formula <- formula(paste("~", paste(colnames(X), collapse = "+")))
            if (is.data.frame(X)) {
              datayx <- list(X = model.matrix(lv.formula, X)[,-1],drop=F)
            } else {
              datayx <- list(X = X)
            }
            lv.X <- as.matrix(model.frame(~ X, data = datayx))
            if(!is.null(row.names(lv.X)))row.names(lv.X)<-row.names(X)
            colnames(lv.X)  <- gsub("X.","",colnames(lv.X))
            lv.formula <- formula(paste("~", paste(colnames(lv.X), collapse = "+")))
            X <- NULL
          }
          
          m1 <- model.frame(y ~ NULL, data = datayx)
          term <- terms(m1)
        }else if(is.null(lv.formula)&!is.null(formula)){
          datayx <- data.frame(y, X)
          m1 <- model.frame(formula, data = datayx)
          term <- terms(m1)
          lv.X <- NULL
          lv.formula <- ~ 1
        } else if(is.null(formula)&!is.null(lv.formula)){
          if(inherits(row.eff,"formula")){
            if(any(colnames(X)==all.vars(row.eff))){
              datayx <- data.frame(y, X[,-which(colnames(X)==all.vars(row.eff)),drop=F])
              X <-  X[,all.vars(row.eff),drop=F]
            } else {
              datayx <- data.frame(y, X)
              X <- NULL
            }
          }else{
            datayx <- data.frame(y, X)
            X <- NULL
          }
          m1 <- model.frame(y ~ NULL, data = datayx)
          labterm <- labels(terms(lv.formula))
          if(any(labterm==1)|any(labterm==0)){
            labterm<-labterm[labterm!=1&labterm!=0]
          }

          lv.formula <- formula(paste("~", paste(labterm, collapse = "+")))
          lv.X<- model.matrix(lv.formula,data=datayx)[,-1,drop=F]
        }else if(!is.null(formula)&!is.null(lv.formula)){
          datayx <- data.frame(y, X)
          m1 <- model.frame(formula, data = datayx)
          labterm <- labels(terms(lv.formula))
          if(any(labterm==1)|any(labterm==0)){
            labterm<-labterm[labterm!=1&labterm!=0]
          }
          lv.formula <- formula(paste("~", paste(labterm, collapse = "+")))
          lv.X<- model.matrix(lv.formula,data=datayx)[,-1,drop=F]
          term <- terms(m1)
        }

      } 
      if(!is.null(X)&(num.lv.c+num.RR)>0|!is.null(data)&(num.lv.c+num.RR)>0){
      if(!is.null(formula)&!is.null(lv.formula)){
        if(any(attr(term,"term.labels")==labterm)){
          stop("Cannot include the same variables for fixed-effects and for constraining the latent variables.")
        }
      }
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
    #If not empty but a vector..
    if(!is.null(disp.formula)){
      if(is.vector(disp.formula)){
        #Defensive coding
          if(length(disp.formula)!=p){
            stop("disp.formula must be a vector of same length as the number of species.")
          } 
        if(any(diff(unique(sort(disp.formula)))!=1)){
          stop("disp.formula indices must form a sequence without gaps.")
        }
        if(min(disp.formula)!=1&max(disp.formula)!=length(unique(disp.formula))){
          stop("disp.formula must start at 1 and end at length(unique(disp.formula)).")
        }
        disp.group <- disp.formula
      }
    }else{
      #if empty we default to the number of species
      disp.group <- 1:NCOL(y)
    }
    
    #check for redundant predictors
    
    if(!is.null(lv.X)){
      if((num.RR+num.lv.c)>ncol(lv.X)){
        stop("Cannot have more reduced dimensions than the number of predictor variables. Please reduce num.RR or num.lv.c \n")
      }
      if((num.RR+num.lv.c)>=p){
        stop("num.RR and num.lv.c should be less than the number of species.")
      }
      #check for redundant predictors
      QR<-qr(lv.X)
      if(QR$rank<ncol(lv.X)){
        warning("Redundant predictors detected, some have been omitted as they explain similar information. \n")
        if(num.lv.c>ncol(lv.X)&num.RR==0){
          num.lv.c <- QR$rank
          warning("Setting num.lv.c to number of non-redunant predictors")
        }else if(num.RR>ncol(lv.X)&num.lv.c==0){
          num.RR <- QR$rank
          warning("Setting num.RR to number of non-redunant predictors")
        }
        if(num.RR>ncol(lv.X)|num.lv.c>ncol(lv.X)){
          stop("Please reduce num.RR and/or num.lv.c, to at maximum the number of non-redundant predictor variables.")
        }
        lv.X.red <- colnames(lv.X)[QR$pivot[-c(1:QR$rank)]]
        lv.X<-lv.X[,QR$pivot[1:QR$rank],drop=F]
        
        #remove redundant terms from formulas
        if(!is.null(lv.formula)){
          lv.formula <- formula(paste("~",paste(attr(terms(lv.formula),"term.labels")[!attr(terms(lv.formula),"term.labels")%in%lv.X.red],collapse="+")))
        }
      }
    }
    
    p <- NCOL(y)
    n <- NROW(y)
    if (p == 1)
      y <- as.matrix(y)

    # Species random effects    
    spdr = NULL; cstruc = "diag"
    if(col.eff == "random") {
      bar.f <- findbars1(col.eff.formula) # list with 3 terms
      if(!is.null(bar.f)) {
        mf <- model.frame(subbars1(col.eff.formula),data=X.col.eff)
        RElist <- mkReTrms1(bar.f,mf)
        spdr <- Matrix::t(RElist$Zt)
        cs <- RElist$cs
        # colnames(spdr) <- rep(names(RElist$nms),RElist$nms)
      }
    }

    if(!inherits(row.eff, "formula") && row.eff == "random"){
      row.eff <- ~(1|site)
      if(is.null(studyDesign)){
        studyDesign <- data.frame(site = 1:n)
      }else{
        studyDesign <- cbind(studyDesign, data.frame(site = 1:n)) 
      }
    }
    
# Structured row parameters
    dr = NULL; cstruc = "diag"
    if(inherits(row.eff,"formula")) {
      bar.f <- findbars1(row.eff) # list with 3 terms
      grps <- unique(unlist(lapply(bar.f, all.vars)))
      if(is.null(studyDesign)){
        if(!is.null(data)) {
        if(any(colnames(data) %in% grps)){
          xgrps<-as.data.frame(data[1:n,(colnames(data) %in% grps)])
          colnames(xgrps) <- grps
          studyDesign<-cbind(studyDesign, xgrps)
        }
          
      } else if(all(grps %in% colnames(X))) {
        if (!is.null(studyDesign)){ 
          studyDesign=cbind(studyDesign, X)
        } else {
          studyDesign=X
        }
        xnames<-colnames(X)[!(colnames(X) %in% grps)]
        X <- as.data.frame(X[,!(colnames(X) %in% grps)])
        colnames(X)<-xnames
        if(ncol(X)==0) X<-NULL
      } else if(is.null(studyDesign)){
        stop("Grouping variable needs to be included in 'studyDesign'")
      }
      } else if(!is.null(studyDesign) && any(colnames(studyDesign) %in% colnames(X))){
        X <- X[,-which(colnames(X)%in%colnames(studyDesign)),drop=F]
        if(ncol(X)==0)X<-NULL
      }
      
      if(is.null(bar.f)) {
        stop("Incorrect definition for structured random effects. Define the structure this way: 'row.eff = ~(1|group)'")
        # } else if(!all(grps %in% colnames(X))) {
        # stop("Grouping variable need to be included in 'X'")
      }
      # else if(!all(apply(studyDesign[,(colnames(studyDesign) %in% grps)],2,order)==c(1:n)) && (corWithin)) {
      #   stop("Data (response matrix Y and covariates X) needs to be grouped according the grouping variable: '",grps,"'")
      # } 
      # if there are nested components, the formula needs to be expanded to ensure 
      # a correct number of entries in corstruc
      form.parts <- strsplit(deparse1(row.eff),split="\\+")[[1]]
      nested.parts <- grepl("/",form.parts)
      if(any(nested.parts)){
      form.parts <- strsplit(deparse1(row.eff),split="\\+")[[1]]
      
      for(i in which(nested.parts))
        form.parts[i] <- paste0("(", findbars1(formula(paste0("~",form.parts[i]))),")",collapse="+")
      
      corstruc.form <- as.formula(paste0("~", paste0(form.parts,collapse="+")))
      }else{
      corstruc.form <- row.eff
      }
      cstruc <- corstruc(corstruc.form)
      corWithin <- ifelse(cstruc == "diag", FALSE, corWithin)
      
      if(!is.null(bar.f)) {
        mf <- model.frame(subbars1(row.eff),data=studyDesign)
        # adjust correlated terms for "corWithin = TRUE"; site-specific random effects with group-specific structure
        # consequence: we can use Zt everywhere
        mf.new <- mf
        if(any(corWithin)){
          mf.new[, corWithin] <- apply(mf[, corWithin, drop=F],2,function(x)order(order(x)))
        }
        colnames(mf.new) <- colnames(mf)
        RElist <- mkReTrms1(bar.f,mf.new)
        dr <- Matrix::t(RElist$Zt)
        colnames(dr) <- rep(names(RElist$grps),RElist$grps)
        
        # add unique column names with corWithin so that we can identify them as separate random effects later
      if(any(corWithin)){
        corWithinNew <- corWithin
        cstrucNew <- cstruc
        for(re.nm in colnames(mf)[corWithin]){
          cstrucNew <- rep(cstrucNew, times = ifelse(colnames(mf)==re.nm, length(unique(mf[, re.nm])), 1))
          corWithinNew <- rep(corWithinNew, times = ifelse(colnames(mf)==re.nm, length(unique(mf[, re.nm])), 1))
          colnames(dr)[colnames(dr) == re.nm] <- paste0(re.nm, sort(mf[, re.nm]))
        }
        corWithin <- corWithinNew
        cstruc <- cstrucNew
      }
      }
      row.eff <- "random"
    }
    
    dLV = NULL;cstruclv = "diag"
    if(inherits(lvCor,"formula")) {
      bar.lv <- findbars1(lvCor) # list with 3 terms
      grps <- unlist(lapply(bar.lv,function(x) as.character(x[[3]])))
      if(is.null(studyDesign)){
        if(all(grps %in% colnames(X))) {
        studyDesign=X
        xnames<-colnames(X)[!(colnames(X) %in% grps)]
        X <- as.data.frame(X[,!(colnames(X) %in% grps)])
        colnames(X)<-xnames
        if(ncol(X)==0) X<-NULL
      } else if(is.null(studyDesign)){
        stop("Grouping variable needs to be included in 'studyDesign'")
      }
      }
      if(!is.null(data)) {
        if(any(colnames(data) %in% grps)){
          xgrps<-as.data.frame(data[1:n,(colnames(data) %in% grps)])
          colnames(xgrps) <- grps
          X<-cbind(X,xgrps)
        }
      }
      if(is.null(bar.lv)) {
        stop("Incorrect definition for structured random effects. Define the structure this way: 'row.eff = ~(1|group)'")
        # } else if(!all(grps %in% colnames(X))) {
        # stop("Grouping variable need to be included in 'X'")
      } else if(!all(order(studyDesign[,(colnames(studyDesign) %in% grps)])==c(1:n)) && (corWithinLV)) {
        stop("Data (response matrix Y and covariates X) needs to be grouped according the grouping variable: '",grps,"'")
      } else {
        # if(quadratic != FALSE) {warning("Structured row effects model may not work properly with the quadratic model yet.")}
        mf <- model.frame(subbars1(lvCor),data=studyDesign)
        dLV <- Matrix::t(mkReTrms1(bar.lv,mf)$Zt)
      }
      cstruclv = corstruc(lvCor)[1]
      num.lv.cor = num.lv
      # Have to set this to diagonal to avoid too many parameters:
      # Lambda.struc = "diagonal"
    } else {
      num.lv.cor = 0
    }
    if(Lambda.struc %in% c("bdNN","UNN") & num.lv.cor>0){
      NN<-t(apply(as.matrix(dist(distLv, upper = TRUE, diag = TRUE)),1, order)[1+(1:NN),])
      i1<-rep(1:nrow(NN), each=ncol(NN))
      i2<-c(t(NN))
      indM<-cbind(i1,i2)
      indM[i1<i2,1]<- i2[i1<i2]
      indM[i1<i2,2]<- i1[i1<i2]
      # indM[,1]>indM[,2]
      indM<-indM[order(indM[,2]),]
      indM<-indM[order(indM[,1]),]
      dupl<-c(TRUE, rowSums(abs(indM[-1,]-indM[1:(nrow(indM)-1),]), na.rm=TRUE)!=0)
      NN<-indM[dupl,]
    } else if(Lambda.struc %in% c("LR") & num.lv.cor>0){
      NN <- as.matrix(NN)
    } else {NN=matrix(0)}
    
    if(num.lv==0&num.lv.c==0&num.RR==0)quadratic <- FALSE

    if(any(colSums(y, na.rm = TRUE) == 0))
      warning("There are responses full of zeros. \n");

    # if(row.eff %in% c("fixed", "random", TRUE) ){
    #   if((p<2) & (any(rstruc == 0)))
    #     stop("There must be at least two responses in order to include unstructured row effects. \n");
    # }
      if(any(rowSums(y, na.rm = TRUE)==0))
        warning("There are rows full of zeros in y. \n");
    # if(row.eff == "random" && quadratic != FALSE && Lambda.struc == "unstructured"){
    #   stop("Dependent row-effects can only be used with quadratic != FALSE if Lambda.struc == 'diagonal'' '. \n")
    #   #This can potentially be relaxed for the gaussian, binomial and ordinal distributions because the linear and quadratic approximation terms can be separated.
    # }
    if( any(!is.finite(y[!is.na(y)])) ) stop("Infinite values are not allowed in 'y'")
    if(any(is.na(y)))y[is.na(y)]<-NA_real_
    if (row.eff == "random" && family == "ordinal" && TMB==FALSE) {
      stop("Random row effect model is not implemented for ordinal family without TMB. \n")
    }
    
    if ((method %in% c("LA", "EVA")) && family == "ordinal") {
      cat("Laplace's (or EVA) method cannot yet handle ordinal data, so VA method is used instead. \n")
      method <- "VA"
    }
    if (method == "LA" && quadratic != FALSE && (num.lv+num.lv.c)>0){
      cat("Laplace's method cannot model species responses as a quadratic function of the latent variables, so attempting VA instead. \n")
      method <- "VA"
    }
    if (method == "VA" && quadratic != FALSE && TMB == FALSE){
      cat("The quadratic model is not implemented without TMB, so 'TMB = TRUE' is used instead. \n")
      TMB <- TRUE
    }
    if (method == "LA" && !TMB) {
      cat("Laplace's method is not implemented without TMB, so 'TMB = TRUE' is used instead. \n")
      TMB = TRUE
    }
    # if ((method %in% c("VA", "EVA")) && (family == "ZIP")) {
    #   cat("VA method cannot handle", family, " family, so LA method is used instead. \n")
    #   method <- "LA"
    # }
    if (quadratic != FALSE && family %in% c("tweedie", "beta")){
      stop("The quadratic model is not implemented for ", family, " family yet. \n")
    }
    if (method == "VA" && family %in% c("tweedie", "beta")){
      cat("Note that, the", family, "family is implemented using the extended variational approximation method. \n")
    }
    # if (method == "EVA"){
    #   method = "VA"
    # }
    if (p < 2 && !is.null(TR)) {
      stop("Fourth corner model can not be fitted with less than tvo response variables.\n")
    }
    if (row.eff == "random" && !TMB) {
      cat("Random row effect model is not implemented without TMB, so 'TMB = TRUE' is used instead. \n")
      TMB <- TRUE
    }
    
    if (family %in% c("gaussian","ZIP","ZINB","beta","Tweedie","gamma","exponential") && !TMB) {
      TMB <- TRUE
      cat("Only TMB implementation available for ", family, " family, so 'TMB = TRUE' is used instead. \n")
    }
    
    # if(family == "ordinal" && num.lv ==0 && zeta.struc == "common"){
    #   stop("Ordinal model with species-common cut-offs without latent variables not yet implemented. Use `TMB = FALSE` and `zeta.struc = `species` instead.")
    # }
    # if(family == "ordinal" && TMB && num.lv==0){
    #   stop("Ordinal model without latent variables not yet implemented using TMB.")
    # }
    
    if(!TMB && zeta.struc == "common"){
      stop("Ordinal model with species-common cut-offs not implemented without TMB.")
    }

    if (!is.null(start.fit)) {
      if (!inherits(start.fit,"gllvm"))
        stop("Only object of class 'gllvm' can be given as a starting parameters.")

      # if (!(family %in% c("poisson", "negative.binomial", "ZIP")))
      #   stop("Starting parameters can be given only for count data.")

    }
    # RRR with num.lv is a special (weird) case
    # where unconstrained LVs are not uncorrelated with predictors
    # better inform the user this might be a bad idea
    if(!is.null(lv.X)){
      if(num.lv.c>0&(num.lv>0)&((num.RR+num.lv.c)<ncol(lv.X))){
        warning("Are you sure you want to fit this model? It might be better to increase num.RR or num.lv.c until the number of predictors is reached, before adding unconstrained LVs. \n")
      }
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
      if (ncol(start.lvs) != (num.lv+num.lv.c) || nrow(start.lvs) != n)
        stop("Given starting value matrix for latent variables has a wrong dimension.")
      if(num.lv>0&num.lv.c==0)colnames(start.lvs) <-  paste("LV",1:num.lv, sep = "")
      if(num.lv==0&num.lv.c>0)colnames(start.lvs) <-  paste("CLV",1:num.lv.c, sep = "")
      if(num.lv>0&num.lv.c>0)colnames(start.lvs) <-  c(paste("CLV",1:num.lv.c, sep = ""),paste("LV",1:num.lv, sep = ""))
    }
 
    n.i <- 1


    out <- list( y = y, X = X, lv.X = lv.X, TR = TR, data = datayx, num.lv = num.lv, num.lv.c = num.lv.c, num.RR = num.RR, num.lvcor =num.lv.cor, lv.formula = lv.formula, lvCor = lvCor, formula = formula,
        method = method, family = family, row.eff = row.eff, col.eff = list(col.eff = col.eff, col.eff.formula = col.eff.formula, spdr = spdr, sp.Ar.struc = sp.Ar.struc),  corP=list(cstruc = cstruc, cstruclv = cstruclv, corWithin = corWithin, corWithinLV = corWithinLV, Astruc=0), dist=dist, distLV = distLV, randomX = randomX, n.init = n.init,
        sd = FALSE, Lambda.struc = Lambda.struc, TMB = TMB, beta0com = beta0com, optim.method=optim.method, disp.group = disp.group, NN=NN, Ntrials = Ntrials, quadratic = quadratic, randomB = randomB)
    if(return.terms) {out$terms = term} #else {terms <- }

    if("la.link.bin" %in% names(pp.pars)){link = pp.pars$la.link.bin}
    if (family == "binomial") {
      if (method %in% c("LA", "EVA"))
        out$link <- link
      if (method == "VA")
        out$link <- "probit"
    }
    if (family %in% c("beta","betaH","orderedBeta")) {
        out$link <- link
    }

    if (family == "orderedbeta") {
      out$link <- "probit"
    }
    out$offset <- offset
    if(quadratic=="LV")start.struc <- "LV"
    if (TMB) {
      if (family == "betaH") {
        if(is.null(colnames(y))) colnames(y)= paste("y",1:NCOL(y))
        y01 = y #(y>0)*1; 
        if(is.null(colnames(y01))) colnames(y01)= paste("y",1:NCOL(y01))
        colnames(y01) = paste("H01",colnames(y01), sep = "_")
        y=cbind(y,y01)
        if(!is.null(TR)){
          TR=rbind(TR,TR)
        }
        if(!is.null(disp.group)){
          disp.group=c(disp.group,disp.group)
        }
        O = cbind(O,O)
      }
      
      # trace = FALSE
      if (row.eff == TRUE)
        row.eff <- "fixed"
      if (!is.null(TR)) {
        fitg <- trait.TMB(
            y,
            X = X,
            # lv.X = lv.X,
            TR = TR,
            formula = formula,
            # lv.formula = lv.formula,
            num.lv = num.lv,
            # num.lv.c = num.lv.c,
            num.lv.cor=num.lv.cor,
            family = family,
            Lambda.struc = Lambda.struc,
            row.eff = row.eff,
            col.eff = col.eff,
            reltol = reltol,
            # reltol.c = reltol.c,
            seed = seed,
            maxit = maxit,
            max.iter=max.iter,
            start.struc = start.struc,
            quad.start = quad.start,
            start.lvs = start.lvs,
            offset = O,
            sd.errors = sd.errors,
            trace = trace,
            link = link,
            n.init = n.init,
            n.init.max = n.init.max,
            start.params = start.fit,
            optimizer = optimizer,
            starting.val = starting.val,
            method = method,
            Power = Power,
            diag.iter = diag.iter,
            Ab.diag.iter = Ab.diag.iter,
            Ab.struct = Ab.struct, Ar.struc = Ar.struc, sp.Ar.struc = sp.Ar.struc,
            dependent.row = dependent.row,
            Lambda.start = Lambda.start,
            jitter.var = jitter.var,
            randomX = randomX,
            randomX.start = randomX.start,
            beta0com = beta0com, 
            scale.X = scale.X,
            zeta.struc = zeta.struc,
            quadratic = quadratic,
            optim.method=optim.method, 
            dr=dr,dLV=dLV, cstruc = cstruc, cstruclv = cstruclv, dist =dist, distLV = distLV, corWithinLV = corWithinLV, NN=NN, 
            scalmax = scalmax, MaternKappa = MaternKappa, rangeP = rangeP,
            setMap = setMap, #Dthreshold=Dthreshold,
            disp.group = disp.group,
            Ntrials = Ntrials, 
            spdr = spdr
        )
        out$X <- fitg$X
        out$TR <- fitg$TR
        out$formula <- fitg$formula
        
      } else {
        fitg <- gllvm.TMB(
            y,
            X = X,
            lv.X = lv.X,
            formula = formula,
            lv.formula = lv.formula,
            num.lv = num.lv,
            num.lv.c = num.lv.c,
            num.RR = num.RR,
            num.lv.cor=num.lv.cor,
            family = family,
            method = method,
            Lambda.struc = Lambda.struc, Ar.struc = Ar.struc, sp.Ar.struc = sp.Ar.struc,
            row.eff = row.eff,
            col.eff = col.eff, colMat = colMat,
            reltol = reltol,
            reltol.c = reltol.c,
            seed = seed,
            maxit = maxit,
            max.iter=max.iter,
            start.struc = start.struc,
            quad.start = quad.start,
            start.lvs = start.lvs,
            offset = O,
            sd.errors = sd.errors,
            trace = trace,
            link = link,
            Ntrials = Ntrials,
            n.init = n.init,
            n.init.max = n.init.max,
            restrict = restrict,
            start.params = start.fit,
            optimizer = optimizer,
            starting.val = starting.val,
            Power = Power,
            diag.iter = diag.iter,
            dependent.row = dependent.row,
            Lambda.start = Lambda.start,
            jitter.var = jitter.var,
            zeta.struc = zeta.struc,
            quadratic = quadratic,
            randomB = randomB,
            optim.method=optim.method, 
            dr=dr, dLV=dLV, cstruc = cstruc, cstruclv = cstruclv, dist =dist, distLV = distLV, corWithinLV = corWithinLV, NN=NN, 
            scalmax = scalmax, MaternKappa = MaternKappa, rangeP = rangeP,
            setMap=setMap, #Dthreshold=Dthreshold,
            disp.group = disp.group,
            spdr = spdr,
            cs = cs
        )
        if(is.null(formula)) {
          out$formula <- fitg$formula
          out$X <- fitg$X
        }
      if(col.eff == "random"){
        if(!is.null(colMat)){
          out$col.eff$colL <- fitg$colL
        }
        out$col.eff$nsp <- fitg$nsp
      }
}
      out$disp.group <- disp.group
      out$seed <- fitg$seed
      out$X.design <- fitg$X.design
      out$TMBfn = fitg$TMBfn
      out$logL <- fitg$logL
      if (num.lv|num.lv.c > 0)
        out$lvs <- fitg$lvs
      # out$X <- fitg$X
      
      out$params <- fitg$params
      
      #### Try to calculate sd errors
      if (!is.infinite(out$logL) && sd.errors) {
        trsd <- try({
          ses <- se.gllvm(out)
          out$sd <- ses$sd
          out$Hess <- ses$Hess
          out$prediction.errors <- ses$prediction.errors
        if(!is.null(out$sd)&(num.lv.c+num.lv)>0|!is.null(out$sd)&row.eff=="random"){
          if(!is.finite(determinant(out$Hess$cov.mat.mod)$modulus)){
            warning("Determinant of the variance-covariance matix is zero. Please double check your model for e.g. overfitting or lack of convergence. \n")
          }
        }
        }, silent = TRUE)
        if(inherits(trsd, "try-error")) { cat("Standard errors for parameters could not be calculated, due to singular fit.\n") }
      }
      
      if (family == "tweedie") {
        out$Power <- fitg$Power
      }
      
      if(family == "ordinal"){
        out$zeta.struc = zeta.struc
      }
      if ((method %in% c("VA", "EVA"))) {
        out$A <- fitg$A
        out$Ar <- fitg$Ar
        out$AQ <- fitg$AQ
        if(randomB!=FALSE){
          out$Ab.lv <- fitg$Ab.lv
        }
        if(col.eff == "random"){
          out$spArs <- fitg$spArs
        }
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
      if ((num.lv+num.lv.c) > 0)
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
    
    #post-hoc processing for sigma.lv
    # if((num.lv.c+num.lv)>0){
    # if(num.lv.c>1){
    #   if(num.lv.c>1){out$params$sigma.lv <- abs(diag(out$params$theta[,1:num.lv.c]));diag(out$params$theta[,1:num.lv.c])<-1}else{out$params$sigma.lv <- abs(out$params$theta[1,1]);out$params$theta[1,1]<-1}
    #   if(!is.null(out$sd)){
    #     if(num.lv.c>1){out$sd$sigma.lv <- diag(out$sd$theta[,1:num.lv.c]);diag(out$sd$theta[,1:num.lv.c])<-0}else{out$sd$sigma.lv <- out$sd$theta[1];out$sd$theta[1]<-0}
    #   }
    # }
    # if(num.lv>0){
    #   if(num.lv.c==0){
    #     if(num.lv>1){out$params$sigma.lv <- abs(diag(out$params$theta[,1:num.lv]));diag(out$params$theta[,1:num.lv])<-1}else{out$params$sigma.lv <- abs(out$params$theta[1,1]);out$params$theta[1,1]<-1}
    #     if(!is.null(out$sd)){
    #       if(num.lv>1){out$sd$sigma.lv <- diag(out$sd$theta[,1:num.lv]);diag(out$sd$theta[,1:num.lv])<-0}else{out$sd$sigma.lv <- out$sd$theta[1];out$sd$theta[1]<-0}
    #     }
    #   }else{
    #     if(num.lv>1){out$params$sigma.lv <- c(out$params$sigma.lv, abs(diag(out$params$theta[,-c(1:num.lv.c),drop=F])));diag(out$params$theta[,-c(1:num.lv.c)])<-1}else{out$params$sigma.lv <- c(out$params$sigma.lv,out$params$theta[1,-c(1:num.lv.c)]);out$params$theta[1,-c(1:num.lv.c)]<-1}
    #     if(!is.null(out$sd)){
    #       if(num.lv>1){out$sd$sigma.lv <- c(out$sd$sigma.lv, diag(out$sd$theta[,-c(1:num.lv.c),drop=F]));diag(out$sd$theta[,-c(1:num.lv.c)])<-0}else{out$sd$sigma.lv <- c(out$sd$sigma.lv,out$sd$theta[1,-c(1:num.lv.c)]);out$sd$theta[1,-c(1:num.lv.c)]<-0}
    #     }
    #   }
    # }
    #   names(out$params$sigma.lv) <- names(out$sd$sigma.lv) <- colnames(out$params$theta[,1:(num.lv+num.lv.c)])
    #   }
    if (family == "negative.binomial")
      out$params$inv.phi <- 1 / out$params$phi
    if (is.infinite(out$logL)){
      warning("Algorithm converged to infinity, try other starting values or different method.")
      cat("Algorithm converged to infinity, try other starting values or different method. \n")
      if(num.lv.c>0|num.RR>0){
        cat("Try scaling and centering your predictors before entering them into the model, if you haven't. \n")
      }
    }
    if(family == "binomial")out$Ntrials = fitg$Ntrials
    if (is.null(out$terms) && return.terms)
      out$terms <- fitg$terms
    if (is.finite(out$logL) && !is.null(TR) && NCOL(out$TR)>0 && NCOL(out$X)>0) {
      out$fourth.corner <- try(getFourthCorner(out),silent = TRUE)
    }
    if(row.eff == "random"){
      out$dr = fitg$dr
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
    
    out$convergence = fitg$convergence
    if(is.finite(out$logL)&TMB){
    if(!out$convergence) {
      warning("The algorithm did not converge, the maximum number of iterations might have been reached.")
      } else if(gradient.check && TMB){
        if(any(abs(c(out$TMBfn$gr(out$TMBfn$par)))> 0.05)) warning("Algorithm converged with large gradients (>0.05). Stricter convergence criterion (reltol) might help. \n")
      }
      # check if constraints on Bs hold, otherwise indicates lack of convergence
      # out$convergence should usually catch this (status code 5 in nloptr, maxeval was reached), but sometimes doesn't.
      if((num.RR+num.lv.c)>1 && out$convergence && randomB==FALSE){
        BB <- t(out$params$LvXcoef)%*%out$params$LvXcoef
        if(any(abs(unique(BB[col(BB)!=row(BB)]))>=1e-2)) warning("Canonical coefficients are not orthogonal, refit the model with a different set of starting values, fit with a different optimizer, or change the optimization criteria for e.g. 'reltol.c'.")

      }
    }


    if(is.null(out$sd)){
      out$sd <- FALSE
    }
    if(TMB && !isFALSE(quadratic)){
    # check if quadratic coefs have converged or have stuck to "LV"
    if(isTRUE(quadratic)){
      if(length(unique(round(out$params$theta[,-c(1:(num.RR+num.lv.c+num.lv)),drop=F],6)))==(num.RR+num.lv.c+num.lv)){
        warning("Quadratic model seems to have converged to species-common tolerances. Try refitting with `start.struc='all'`, or with different starting values.\n")
      out$quadratic <- "LV"        
      }else if(length(unique(out$params$theta[,-c(1:(num.RR+num.lv.c+num.lv)),drop=F]))==1 && starting.val == "zero"){
        warning("It looks like the optimizer failed to move the quadratic coefficients away from the starting values. Please change the starting values. \n")
        out$quadratic <- quadratic
      }else{
        out$quadratic <- TRUE
      }
    }else if(quadratic == "LV"){
      if(length(unique(out$params$theta[,-c(1:(num.RR+num.lv.c+num.lv)),drop=F]))==1 && starting.val == "zero"){
        warning("It looks like the optimizer failed to move the quadratic coefficients away from the starting values. Please change the starting values. \n")
      }
      out$quadratic <- quadratic
    }
    }else{
      out$quadratic <- FALSE
    }
    if(!TMB&family=="ordinal"){
      out$zeta.struc <- "species"
    }else if(TMB & family == "ordinal"){
      out$zeta.struc = fitg$zeta.struc
    }

    out$call <- match.call()
    if(quadratic == FALSE){
      class(out) <- "gllvm"  
    }else{
      class(out) <- c("gllvm","gllvm.quadratic")
    }
    
    return(out)
  }
