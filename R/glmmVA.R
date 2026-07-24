#' Fitting Generalized Linear Mixed-Effects Models with Variational Approximation
#'
#' Wraps the \code{"\link{gllvm}"} function to fit a (univariate) generalized linear mixed-effects model (GLMM) in a more familiar syntax. 
#' Both fixed and random effects are specified in \code{lme4}-style via the formula argument.
#'
#' @param formula a formula object describing both the fixed- and random-effects part of the model. A response should be present on the left-hand side of the operator, and otherwise passed as a named 'y' argument to the function. Random effects are written in \code{lme4}-style. Some structured random effects are supported, see \code{"\link{gllvm}"} for more information.
#' @param data an optional data frame containing the variables named in formula.
#' @param family a family as supported by the \code{"\link{gllvm}"} function. For mixed response type models, a character vector of length equal to the number of observations.
#' @param response.group optional factor of length equal to the number of observations. Observations sharing a level are mapped to the same distribution.
#' @param control A list with the following arguments controlling the optimization:
#' \describe{
#'  \item{\emph{reltol}: }{ convergence criteria for log-likelihood, defaults to 1e-10.}
#'  \item{\emph{optimizer}: }{the log-likelihood can be optimized using \code{"\link{optim}"} (default) or \code{"\link{nlminb}"}.}
#'  \item{\emph{max.iter}: }{ maximum number of iterations for \code{optimizer = "nlminb"}, defaults to 6000.}
#'  \item{\emph{maxit}: }{ maximum number of iterations for optimizer, defaults to 6000.}
#'  \item{\emph{optim.method}: }{ optimization method to be used if optimizer is \code{"\link{optim}"}. Defaults to \code{"BFGS"}, but to \code{"L-BFGS-B"} for Tweedie family due the limited-memory use.}
#' }
#' @param control.va A list with the following arguments controlling the variational approximation method:
#' \describe{
#'  \item{\emph{Ar.struc}: }{ covariance structure of VA distributions, "unstructured" or "diagonal". Defaults to "unstructured". "Unstructured" meaning block diagonal for ordinary random effects, a kronecker product for propto structures with correlaton parameters, and fully unstructured for structured random effects (such as "corExp").}
#'  \item{\emph{diag.iter}: }{ non-negative integer which can sometimes be used to speed up the updating of variational (covariance) parameters, whichh can sometimes improve the fit. Either 0 or 1. Defaults to 0.}
#'  \item{\emph{Lambda.start}: }{ starting value for variances in VA distributions. Defaults to 0.3.}
#' }
#' @param control.start A list with the following arguments controlling the starting values:
#' \describe{
#'   \item{\emph{starting.val}: }{ defaults to \code{"zero"}. See \code{"\link{gllvm}"} for details.}
#'   \item{\emph{n.init}: }{ number of initial runs. Uses multiple runs and picks up the one giving highest log-likelihood value. Defaults to 1.}
#'   \item{\emph{n.init.max}: }{ maximum number of refits try try for n.init without improvement, defaults to 10.}
#'   \item{\emph{start.fit}: }{ object that inherits of class 'gllvm' which can be given as starting parameters.}
#'   \item{\emph{MaternKappa}: }{ Starting value for smoothness parameter of Matern covariance function. Defaults to 3/2.}
#'   \item{\emph{scalmax}: }{ Sets starting value for the scale parameter for the coordinates. Defaults to 10, when the starting value for scale parameter scales the distances of coordinates between 0-10.}
#'   \item{\emph{rangeP}: }{ Sets starting value for the range parameter for the correlation structure.}
#'   \item{\emph{zetacutoff}: }{ A vector of length 2. Sets starting value for the cutoff parameters of the ordered beta model.}
#'   \item{\emph{start.optimizer}: }{ optimizer for starting value generation, see "optimizer" for more information.}
#'   \item{\emph{start.optim.method}: }{ optimizer method for starting value generation, see "optim.method" for more information.}
#' }
#' @param ... other arguments passed onto the \code{"\link{gllvm}"} function.
#' 
#' @return An object of class "glmmVA" that inherits from the "gllvm" class.
#' 
#' @seealso \code{\link{gllvm}}
#' 
#' @author Bert van der Veen
#'
#' @examples
#' data(eSpider)
#' data <- cbind(data.frame(y = c(eSpider$abund[eSpider$nonNA,]),
#'            species = factor(rep(1:ncol(eSpider$abund), each = length(eSpider$nonNA))),
#'            site = factor(rep(1:length(eSpider$nonNA), ncol(eSpider$abund)))),
#'            do.call(rbind, replicate(ncol(eSpider$abund), scale(eSpider$X[eSpider$nonNA,]), 
#'            simplify = FALSE)))
#' # Example 1: crossed random slope effects                          
#' model <- glmmVA(y~species + ConWate + ConHumu + (0+ConHumu|species) + (0+ConWate|species),
#'                 family = "poisson", data = data)
#'
#' # Example 2: correlated random slopes
#' model1 <- glmmVA(y~species + ConWate + ConHumu + (0+ConWate+ConHumu|species),
#'                  family = "poisson", data = data)
#'
#' # Example 3: joint model for mixed Poisson and and negative-binomial responses.
#' fam <- c(rep("poisson", 28), rep("negative.binomial", nrow(data)-28))
#' model2 <- update(model1, family=fam, sd.errors=FALSE)
#'
#' @export
glmmVA <- function(formula, data, family, response.group = NULL,
                   control = list(reltol = 1e-10, optimizer = "optim", max.iter = 6000, maxit = 6000, optim.method = "BFGS"),
                   control.va = list(Ar.struc="unstructured", diag.iter = 0, Lambda.start = 0.3),
                   control.start = list(starting.val = "zero", n.init = 1, n.init.max = 10, start.fit = NULL, scalmax = 10, MaternKappa=1.5, rangeP=NULL, zetacutoff = NULL, start.optimizer = "nlminb", start.optim.method = "BFGS"),
                   ...) {

  args  <- list(...)

  # Auto-infer response.group when family is a per-observation vector
  mf_tmp <- model.frame(subbars1(formula), data = data)
  y_tmp  <- model.response(mf_tmp)
  n_tmp  <- if (is.null(y_tmp)) length(args[["y"]]) else NROW(y_tmp)
  if (is.null(response.group) && is.character(family) && length(family) == n_tmp) {
    response.group <- as.factor(family)
  }
  rm(mf_tmp, y_tmp, n_tmp)

  mixed_mode <- !is.null(response.group)
  
  mf <- model.frame(subbars1(formula), data = data)
  y <- as.matrix(model.response(mf))
  offset <- model.offset(mf)
  if(!is.null(offset)) offset <- matrix(offset, nrow = NROW(y), ncol = NCOL(y))
  row.eff.formula = remove.offset(formula)
  
  if("Lambda.start" %in% names(control.va)){
    control.va$Lambda.start <- c(0.3, control.va$Lambda.start[1], 0.3)
  }
  
  if(is.null(y) && !"y" %in% names(args)){
    stop("Response variable not found")
  }else if("y" %in% names(args)){
    y <- args[[names(args)]]
    args[[-which("y"%in%names(args))]]
  }else{
    row.eff.formula <- row.eff.formula[-2]
  }

  # Mixed response type: reshape long-format y into a wide matrix (n_obs x n_groups)
  # with NA in cells where an observation does not belong to that response group.
  if (mixed_mode) {
    response.group <- as.factor(response.group)
    if (length(response.group) != nrow(y))
      stop("`response.group` must have the same length as the number of observations.")

    groups <- levels(response.group)
    n_obs  <- nrow(y)
    n_grp  <- length(groups)

    # Resolve family to a per-column character vector
    if (length(family) == n_obs) {
      family_col <- character(n_grp)
      for (g in seq_along(groups)) {
        idx  <- which(response.group == groups[g])
        fams <- unique(family[idx])
        if (length(fams) != 1L)
          stop(sprintf("All observations in group '%s' must share the same family.", groups[g]))
        family_col[g] <- fams
      }
      family <- family_col
    } else if (length(family) == n_grp) {
      # already per-column — nothing to do
    } else if (length(family) == 1L) {
      family <- rep(family, n_grp)
    } else {
      stop("`family` must have length 1, nrow(data), or nlevels(response.group) in mixed-response mode.")
    }

    # Build wide matrix: n_obs rows x n_grp columns, NA where obs is not in that group
    y_wide <- matrix(NA_real_, nrow = n_obs, ncol = n_grp,
                     dimnames = list(NULL, paste0("y", seq_len(n_grp))))
    for (g in seq_along(groups)) {
      idx <- which(response.group == groups[g])
      y_wide[idx, g] <- y[idx, 1L]
    }
    y <- y_wide
  }

  X = model.frame(subbars1(row.eff.formula), data)
  
  if (inherits(family,"family") || is.function(family) && inherits(family(), "family")) {
    if(inherits(family, "family")){
    link <- family$link
    family <- family$family
    }else if(is.function(family) && inherits(family(), "family")){
      link <- family()$link
      family <- family()$family
    }
  }else if("link" %in% names(args)){
    link <- args[["link"]]
    if(length(args)>1){
    args <- args[[names(args)!="link"]]
    }else{
      args <- NULL
    }
  }else{
    link <- "probit"
  }
  
  Ntrials = matrix(1)

  # Facilitate y to be passed as a 2-column matrix for binomial responses (univariate only)
  if (!mixed_mode) {
    if((ncol(y) == 2) && family %in% c("binomial", "ZIB", "ZNIB", "beta.binomial") && !("Ntrials" %in% names(args))){
      Ntrials <- matrix(rowSums(y))
      y <- y[, 1,drop=FALSE]
    }else if((ncol(y)==2) && family %in% c("binomial", "ZIB", "ZNIB", "beta.binomial") && ("Ntrials" %in% names(args))){
      stop("Cannot both pass 'y' as a matrix and provide the 'Ntrials' argument.")
    }
  }
  
  # Detect no-intercept formula (0+ or -1 in fixed part) before deciding beta0com.
  fixed_part <- nobars1_(row.eff.formula)
  no_intercept <- inherits(fixed_part, "formula") &&
    !attr(terms(fixed_part), "intercept")

  # In mixed mode, collapse the p per-column intercepts to one shared intercept
  # so the count matches the univariate case (one beta0, not p).
  # Skip when no_intercept: fixing all p to 0 is equivalent to sharing one zero.
  if (mixed_mode && !no_intercept && !"beta0com" %in% names(args)) {
    args$beta0com <- TRUE
  }

  # If the user specified no intercept (0+ or -1), fix b to 0 via setMap.
  # Use ncol(y) so the map length matches the number of intercept parameters
  # (p=1 for univariate, p=n_grp for mixed), avoiding the beta0com else-branch conflict.
  if (no_intercept) {
    map_b <- factor(rep(NA_integer_, ncol(y)))
    if(!"setMap" %in% names(args)) {
      args$setMap <- list(b = map_b)
    } else if(!"b" %in% names(args$setMap)) {
      args$setMap[["b"]] <- map_b
    }
  }

  row.eff_gllvm <- if(no_intercept) update(row.eff.formula, ~ . + 1) else row.eff.formula

  gllvm_args <- c(list(y = y, studyDesign = X, family = family, num.lv = 0,
                       row.eff = row.eff_gllvm, offset = offset, link = link,
                       Ntrials = Ntrials, control = control,
                       control.va = control.va, control.start = control.start),
                  args)
  if (mixed_mode) {
    object <- withCallingHandlers(
      do.call(gllvm, gllvm_args),
      warning = function(w) {
        if (grepl("rows full of zeros", conditionMessage(w)))
          invokeRestart("muffleWarning")
      }
    )
  } else {
    object <- do.call(gllvm, gllvm_args)
  }

  # b is fixed at 0 via setMap; results in NA
  if(no_intercept) object$params$beta0 <- rep(0, length(object$params$beta0))

  object$call <-  match.call()
  class(object) <-  c("glmmVA", "gllvm")

  return(object)
}

#' @title Predict method for glmmVA objects
#'
#' @description Wrapper around \code{\link{predict.gllvm}} that, for
#'   mixed-response models (those fitted with a per-observation \code{family}
#'   vector), collapses the internal wide-format prediction matrix back to a
#'   long-format numeric vector aligned with the original observations.
#'
#' @param object an object of class \code{"glmmVA"}.
#' @param ... further arguments passed to \code{\link{predict.gllvm}}.
#'
#' @return For univariate models: the same object returned by
#'   \code{\link{predict.gllvm}}. For mixed-response models:
#'   a numeric vector of length \code{nobs(object)}.
#'
#' @seealso \code{\link{predict.gllvm}}, \code{\link{glmmVA}}
#' @author Bert van der Veen
#' @export
predict.glmmVA <- function(object, ...) {
  out <- predict.gllvm(object, ...)

  # For mixed-response models the wide matrix (n_obs x n_grp) has the same
  # number of rows as the original long-format data; collapse to a vector by
  # taking the non-NA column for each observation.
  if (any(is.na(object$y))) {
    col_idx <- apply(object$y, 1L, function(row) which(!is.na(row))[1L])
    idx2    <- cbind(seq_len(nrow(object$y)), col_idx)
    if (is.matrix(out)) {
      out <- out[idx2]
    } else if (is.list(out) && is.matrix(out$fit)) {
      out$fit   <- out$fit[idx2]
      out$lower <- out$lower[idx2]
      out$upper <- out$upper[idx2]
    }
  }
  out
}

#' @title Generic function for extracting random effects
#' @description S3 generic dispatched on the model class.
#' @param object a fitted model object of class glmmVA.
#' @param ... additional arguments passed to the method.
#' @export
ranef <- function(object, ...) UseMethod("ranef")

#' @title Extract random effects from a glmmVA object
#'
#' @description Returns the random effects estimate for a model fitted
#'   with \code{\link{glmmVA}}, optionally with their conditional variances.
#'
#' @param object an object of class \code{"glmmVA"}.
#' @param condVar logical; if \code{TRUE}, attach conditional variances as an
#'   attribute \code{"condVar"} on the returned object. Defaults to \code{FALSE}.
#' @param ... not used.
#'
#' @return A vector of random-effect estimates. When \code{condVar = TRUE}, the attribute \code{"condVar"}
#'   is a vector of conditional variances.
#'   returned by \code{\link{getPredictErr}}.
#'
#' @seealso \code{\link{glmmVA}}, \code{\link{getPredictErr}}
#'
#' @author Bert van der Veen
#'
#' @export
ranef.glmmVA <- function(object, condVar = FALSE, ...) {
  est <- coef(object, "row.params.random")
  if (condVar) {
    cv <- getPredictErr(object, cov = TRUE, ...)$row.effects
    attr(est, "condVar") <- cv
  }
  est
}
