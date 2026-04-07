#' @title Predict species richness from a gllvm object
#' @description Obtains predictions from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param spp index vector, defaults to NULL. If omitted, returns a prediction with all species in the data.
#' @param expected character, defaults to "mean", in which case the returned measure of central tendency for species richness is the Poisson-Binomial expectation. Alternatively can be "mode".
#' @param se.fit integer or logical, defaults to 1000. Number of simulations for confidence interval. No confidence interval is returned when set to \code{FALSE}.
#' @param level specification for how to predict. Level one (\code{level = 1}) attempts to use the predicted site scores from variational approximations or laplace approximation or given site scores in \code{newLV}. Level 0 sets the latent variable to zero. Defaults to 1.
#' @param ci character vector, defaults to "expected" but can also be "pmf" or both. If "pmf" provides CI on the whole pmf prediction.
#' @param alpha numeric between 0 and 1, defaults to 0.95. Confidence level of the prediction.
#' @param seed numeric, defaults to 42. Seed for the simulation of the confidence interval.
#' @param return.pred logical, defaults to \code{FALSE}. Returns the point-estimate results from \code{"\link{predict.gllvm}"}.
#' @param batch integer or \code{NULL}, defaults to 200. If provided, species predictions needed to calculate richness are processed in chunks of this size per simulation draw, reducing peak memory at the cost of more \code{predict} calls. \code{NULL} processes all (selected) species in one call.
#' @param ... arguments passed to \code{"\link{predict.gllvm}"}.
#'
#' @details This function returns probabilities of richness for species richness from gllvm objects, by first calling \code{"\link{predict.gllvm}"}, and following by predicting from a Poisson-Binomial distribution for richness. The distribution for species richness follows from a sum of binary responses (i.e., occurrence), and is naturally extended to non-normal data types as the probability to get a non-zero observation. Especially with many species and rows, the calculation may take a while.
#' @note The point estimate ("fit") can sometimes be outside of the simulated intervals when the model fit is poor, the standard errors of parameter estimates are on very different scales, or the number of simulations is too low. Try either 1) scaling and centering covariates in the model, 2) repeatedly refitting the model (see the n.init argument in \code{"\link{gllvm}"}) to find a better set of starting values and improve the fit, or  3) increasing the number of simulations via 'n.init'.
#'
#' @return list with entries "predicted" and "expected". "predicted" includes the prediction from the Poisson-Binomial distribution, returning a matrix of size sites by length(SR), with statistical uncertianty if se.fit = TRUE. "expected" includes the expected species' richness, i.e., a vector of size n, with statistical uncertianty if se.fit = TRUE.
#' 
#' @seealso \code{\link{predict.gllvm}}, \code{\link{residuals.predictSR.gllvm}}
#' 
#' @author Bert van der Veen
#'
#' @examples
#' \donttest{
#'# Load a dataset from the mvabund package
#'data(spider, package = "mvabund")
#'y <- as.matrix(spider$abund)
#'X <- scale(spider$x)
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, formula = ~soil.dry, family = poisson())
#'# Fitted values
#'newX = data.frame(soil.dry = seq(min(X[,"soil.dry"]),max(X[,"soil.dry"]),length.out=100))
#'predSR <- predictSR(fit, newX = newX, level = 0)
#'
#'# Visualize the results
#' par(mfrow = c(4,4))
#' for(i in 0:ncol(fit$y)){
#' plot(x = newX$soil.dry, y = predSR$predicted$fit[,i+1], xlab = "Soil dry matter content", ylab = bquote(p(SR == .(i))), type = "l", ylim = c(0,1))
#' }
#'
#'}
#'
#'@author Bert van der Veen
#'
#'@aliases predictSR predictSR.gllvm
#'@method predictSR gllvm
#'@export
#'@export predictSR.gllvm

predictSR.gllvm <- function(object, spp = NULL, expected = "mean", se.fit = 1000, level = 1, ci = "expected", alpha = 0.95, seed = 42, return.pred = FALSE, batch = NULL, ...){

  fit <- predict(object, type = "response", se.fit = FALSE, ordinal.cat = 1L, spp = spp, level = level, ...)

  probs <- gllvm.presence.prob(fit, object, spp = spp)

  if(is.null(spp)) spp <- seq_len(ncol(probs))
  SR  <- 0:ncol(probs)
  n   <- nrow(probs)

  eSR <- rowSums(probs)
  if(expected == "mode") eSR <- sapply(eSR, pbMode, n = length(spp))

  TMBcores    <- TMB::openmp(DLL = "gllvm")[1]
  PoisBinCores <- get_omp_threads()
  if(PoisBinCores == 1 && TMBcores > 1) set_omp_threads(TMBcores)

  predSR <- poisbinom(probs)
  colnames(predSR) <- paste0("SR_", SR)

  out <- list(predicted = list(fit = predSR), expected = list(fit = eSR))

  if(is.numeric(se.fit) || isTRUE(se.fit)){
    R <- if(is.numeric(se.fit)) se.fit else 1000L

    params <- simulate.params.gllvm(object, R, seed, level = level, n = n)

    type <- if(expected == "mode") 1L else 7L

    do_expected <- "expected" %in% ci
    do_pmf      <- "pmf" %in% ci

    if(do_expected) eSRmat     <- matrix(0.0, nrow = R, ncol = n)
    if(do_pmf)      predSR.sim <- array(0, dim = c(R, n, length(SR)))

    batches <- if(do_pmf || is.null(batch)) list(spp) else
                 split(spp, ceiling(seq_along(spp) / batch))

    # pre-compute constant objects to avoid repeating work inside the loop
    skeleton     <- object$TMBfn$env$parList()
    obj_template <- object
    obj_template$params <- lapply(object$params, function(x){ if(is.numeric(x)) x[] <- NA; x })
    if(!is.null(obj_template$lvs)) obj_template$lvs[] <- NA

    for(r in seq_len(R)){
      newobj <- perturb.gllvm(object, params, r, type = "response",
                              skeleton = skeleton, template = obj_template)

      if(do_pmf){
        pred_r  <- predict(newobj, type = "response", se.fit = FALSE,
                           ordinal.cat = 1L, spp = spp, ...)
        probs_r <- gllvm.presence.prob(pred_r, object, spp = spp)
        predSR.sim[r, , ] <- poisbinom(probs_r)
        if(do_expected){
          eSRmat[r, ] <- rowSums(probs_r)
          if(expected == "mode")
            eSRmat[r, ] <- sapply(eSRmat[r, ], pbMode, n = length(spp))
        }
      } else if(do_expected){
        for(b in batches){
          pred_r <- predict(newobj, type = "response", se.fit = FALSE,
                            ordinal.cat = 1L, spp = b, ...)
          eSRmat[r, ] <- eSRmat[r, ] + rowSums(gllvm.presence.prob(pred_r, object, spp = b))
        }
        if(expected == "mode")
          eSRmat[r, ] <- sapply(eSRmat[r, ], pbMode, n = length(spp))
      }
    }

    if(do_expected){
      eSR.CI <- apply(eSRmat, 2, quantile,
                      prob = c((1 - alpha) / 2, 1 - (1 - alpha) / 2),
                      type = type)
      out$expected$lower <- eSR.CI[1, ]
      out$expected$upper <- eSR.CI[2, ]
    }

    if(do_pmf){
      SR.ci <- array(dim = c(2, n, length(SR)))
      for(i in seq_len(n)){
        center    <- pmin(pmax(predSR[i, ], .Machine$double.eps), 1 - .Machine$double.eps)
        sim_i     <- pmin(pmax(predSR.sim[, i, ], .Machine$double.eps), 1 - .Machine$double.eps)
        dists     <- hilbert_to_provided_center(sim_i, center)
        threshold <- quantile(dists, alpha)
        SR.ci[, i, ] <- apply(predSR.sim[dists <= threshold, i, , drop = FALSE], 3, range)
      }
      out$predicted$lower <- SR.ci[1, , ]
      out$predicted$upper <- SR.ci[2, , ]
    }
  }
  if(return.pred){
    out$predict.gllvm <- list(fit = fit)
  }

  out$spp <- spp
  
  class(out) <- "predictSR.gllvm"
  return(out)
  
}

#'@export
predictSR <- function(object, ...) {
  UseMethod("predictSR")
}

hilbert_to_provided_center <- function(mat, center) {
  # Ensure inputs are valid
  if (any(mat[!is.na(mat)] <= 0) || any(center[!is.na(center)] <= 0)) {
    stop("Hilbert metric requires strictly positive values (interior of the simplex).")
  }
  
  if (ncol(mat) != length(center)) {
    stop("Dimension mismatch: 'center' must have the same length as the number of columns in 'mat'.")
  }
  
  ratios <- sweep(mat, 2, center, "/")
  
  row_max <- apply(ratios, 1, max, na.rm = TRUE)
  row_min <- apply(ratios, 1, min, na.rm = TRUE)
  
  return(log(row_max / row_min))
}

#' @title Predict co-occurrence from a gllvm object
#' @description Obtains the probability of co-occurrnece for species from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param spp a vector of length 2, defaults to NULL. If NULL, returns the probability of co-occurrence for all species in the data.
#' @param se.fit logical, defaults to \code{FALSE}. If set to \code{TRUE}, returns a simulated confidence interval of the prediction.
#' @param alpha numeric between 0 and 1, defaults to 0.95. Confidence level of the prediction.
#' @param seed numeric, defaults to 42. Seed for the simulation of the confidence interval.
#' @param ... arguments passed to \code{"\link{predict.gllvm}"}.
#'
#' @details
#' The probability of co-occurrence is simply given by the joint probability of getting two presences for any pair of species. As such,
#' it is simply calculated as the product of the probability of occurrence for any pair of species, so that this function mostly serves
#' as convenience wrapper around \code{"\link{predict.gllvm}"}.
#' 
#' @return A matrix of size np(p-1)/2 by 2.
#' 
#' @seealso \code{\link{predict.gllvm}}, \code{\link{predictSR.gllvm}}
#' 
#' @author Bert van der Veen
#'
#' @examples
#' 
#' \donttest{
#'# Load a dataset from the mvabund package
#'data(spider, package = "mvabund")
#'y <- as.matrix(spider$abund)
#'X <- scale(spider$x)
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, formula = ~soil.dry, family = poisson())
#'# fitted values
#'newX = data.frame(soil.dry = seq(min(X[,"soil.dry"]),max(X[,"soil.dry"]),length.out=100))
#'predPair <- predictPairwise(fit, spp = c(1,2), newX = newX, level = 0)
#'plot(x = newX$soil.dry, y = predPair$prob, type = "l", xlab = "Soil dry moisture content", ylab = "Probability of co-occurrence", ylim = c(0,1))
#'
#'}
#'
#'@author Bert van der Veen
#'
#'@aliases predictPairwise predictPairwise.gllvm
#'@method predictPairwise gllvm
#'@export
#'@export predictPairwise.gllvm
predictPairwise.gllvm <- function(object, spp = NULL, se.fit = FALSE, alpha = 0.95, seed = 42, ...){
  
  preds <- predict(object, type = "response", se.fit = se.fit, alpha = NULL, ...)
  fit <- preds
  
  if(is.null(spp)){spp <- 1:p}else{spp <- sort(spp)}
  spp.pairs <- t(combn(spp, 2))
  
  if(is.list(preds)){
    fit <- preds$fit
    
    probs.fit <- gllvm.presence.prob(fit, object, spp)
    n <- nrow(probs.fit)
    p <- ncol(probs.fit)
    
    colnames(probs.fit) <- spp
    idx <- matrix(match(spp.pairs, colnames(probs.fit)), ncol = 2)
    
    probs.fit <- cbind(
      probs.fit[, idx[,1 ], drop = FALSE],
      probs.fit[, idx[,2 ], drop = FALSE]
    )
    probs.fit <- matrix(probs.fit, ncol = 2)
    
    predPair <- matrix(exp(rowSums(log(probs.fit))), ncol = 1) # rowProd via exp-sum-log
    
    predPair.sim <- matrix(nrow=dim(preds$ci.sim)[1], ncol =nrow(predPair))
    for(r in 1:dim(preds$ci.sim)[1]){
      probs.fit.sim <- gllvm.presence.prob(preds$ci.sim[r,,], object, spp)
      colnames(probs.fit.sim) <- spp
      probs.fit.sim <- cbind(
        probs.fit.sim[, idx[,1 ], drop = FALSE],
        probs.fit.sim[, idx[,2 ], drop = FALSE]
      )
      probs.fit.sim <- matrix(probs.fit.sim, ncol = 2)
      predPair.sim[r,] <-  exp(rowSums(log(probs.fit.sim))) # rowProd via exp-sum-log
      
    }
    predPair.ci <- apply(predPair.sim, 2, quantile, prob = c((1-alpha)/2, 1-(1-alpha)/2))
    out <- data.frame(prob = predPair,  site = rep(seq_len(n), times = nrow(spp.pairs)), spp1 = rep(spp.pairs[,1], each = n), spp2 = rep(spp.pairs[,2], each =n), lower = predPair.ci[1,], upper = predPair.ci[2,])
  }else{
    probs.fit <- gllvm.presence.prob(fit, object, spp)
    n <- nrow(probs.fit)
    
    colnames(probs.fit) <- spp
    idx <- matrix(match(spp.pairs, colnames(probs.fit)), ncol = 2)
    
    probs.fit <- cbind(
      probs.fit[, idx[,1 ], drop = FALSE],
      probs.fit[, idx[,2 ], drop = FALSE]
    )
    probs.fit <- matrix(probs.fit, ncol = 2)
    
    predPair <- matrix(exp(rowSums(log(probs.fit))), ncol = 1)
    
    out <- data.frame(prob = predPair,  site = rep(seq_len(n), times = nrow(spp.pairs)), spp1 = rep(spp.pairs[,1], each = n),spp2 = rep(spp.pairs[,2], each =n))
  }
  return(out)
}

#'@export
predictPairwise <- function(object, ...) {
  UseMethod("predictPairwise")
}

# prediction of binomial-based distributions can be used as-is (including orderedBeta)
# not implemented because not based (indirectly) on binary data: gaussian, beta, betaH, gamma, exponential
# families implemented: poisson, NB, NB1, ZIP, ZINB, Tweedie, ordinal
gllvm.presence.prob <- function(fit, object, spp = NULL) {

  family <- object$family
  n <- nrow(fit)
  if(is.null(spp)) spp <- seq_len(ncol(object$y))
  probs <- matrix(0, nrow = n, ncol = length(spp))

  # fit may already be pre-subsetted to length(spp) columns (e.g. from predict(..., spp=))
  # in that case use sequential column indices into fit, but spp[iter] for object lookups
  presubsetted <- ncol(fit) == length(spp) && ncol(fit) < ncol(object$y)

  fam_spp <- family[spp]

  for(fam in unique(fam_spp)){
    pos   <- which(fam_spp == fam)              # positions in spp / probs
    j     <- spp[pos]                           # original species indices for object lookups
    fcols <- if(presubsetted) pos else j        # column indices into fit
    mu    <- fit[, fcols, drop = FALSE]

    if(fam == "binomial"){
      Ntrials_mat <- matrix(object$Ntrials[, j, drop = FALSE], nrow = n, ncol = length(pos))
      probs[, pos] <- 1 - pbinom(0, size = Ntrials_mat, prob = mu)
    } else if(fam == "ZIB"){
      # params$phi for ZIB stores the zero-inflation probability directly (see gllvm.TMB.R line ~1669)
      sigma_vec <- rep(object$params$phi[j], each = n)
      Ntrials_vec <- c(object$Ntrials[, j, drop = FALSE])
      probs[, pos] <- matrix(1 - pzib(0, mu = as.vector(mu), sigma = sigma_vec, Ntrials = Ntrials_vec), nrow = n)
    } else if(fam == "ZNIB"){
      # params$phi stores p0 and ZINB.phi stores pN (already on probability scale, see residuals.gllvm)
      phis0 <- (object$params$phi[j]      / (1 + object$params$phi[j] + object$params$ZINB.phi[j]))
      phisN <- (object$params$ZINB.phi[j] / (1 + object$params$phi[j] + object$params$ZINB.phi[j]))
      p0_vec  <- rep(phis0, each = n)
      pN_vec  <- rep(phisN, each = n)
      Ntrials_vec <- c(object$Ntrials[, j, drop = FALSE])
      probs[, pos] <- matrix(1 - pznib(0, mu = as.vector(mu), p0 = p0_vec, pN = pN_vec, Ntrials = Ntrials_vec), nrow = n)
    } else if(fam == "poisson"){
      probs[, pos] <- ppois(0, lambda = mu, lower.tail = FALSE)
    } else if(fam == "negative.binomial"){
      sizes <- matrix(1 / object$params$phi[j], nrow = n, ncol = length(pos), byrow = TRUE)
      probs[, pos] <- pnbinom(0, mu = mu, size = sizes, lower.tail = FALSE)
    } else if(fam == "negative.binomial1"){
      sizes <- mu * matrix(object$params$phi[j], nrow = n, ncol = length(pos), byrow = TRUE)
      probs[, pos] <- pnbinom(0, mu = mu, size = sizes, lower.tail = FALSE)
    } else if(fam == "betaH"){
      probs[, pos] <- 1 - fit[, -seq_len(ncol(object$y)), drop = FALSE][, j, drop = FALSE]
    } else if(fam == "orderedBeta"){
      linkinv <- binomial(link = object$link)$linkinv
      linkfun <- binomial(link = object$link)$linkfun
      zeta1 <- if(object$zeta.struc == "species")
        matrix(object$params$zeta[j, 1], nrow = n, ncol = length(pos), byrow = TRUE)
      else
        object$params$zeta[1]
      probs[, pos] <- 1 - linkinv(zeta1 - linkfun(mu))
    } else if(fam == "ZIP"){
      phi_vec <- rep(object$params$phi[j], each = n)
      probs[, pos] <- matrix(1 - pzip(0, mu = as.vector(mu), sigma = phi_vec), nrow = n)
    } else if(fam == "ZINB"){
      phi_vec  <- rep(object$params$phi[j],      each = n)
      zphi_vec <- rep(object$params$ZINB.phi[j], each = n)
      probs[, pos] <- matrix(1 - pzinb(0, mu = as.vector(mu), p = phi_vec, sigma = zphi_vec), nrow = n)
    } else if(fam == "tweedie"){
      phi_vec <- rep(object$params$phi[j], each = n)
      probs[, pos] <- matrix(1 - fishMod::pTweedie(0, mu = as.vector(mu),
                                                    phi = phi_vec, p = object$Power), nrow = n)
    } else if(fam == "ordinal"){
      probs[, pos] <- 1 - mu # mu = p(y;k=1), i.e., 1- probability of absence
    } else if(fam == "beta.binomial"){
      phi_vec <- rep(object$params$phi[j], each = n)
      Ntrials_vec <- c(object$Ntrials[, j, drop = FALSE])
      probs[, pos] <- matrix(1 - pbetabinom(0, mu = as.vector(mu), phi = phi_vec, Ntrials = Ntrials_vec), nrow = n)
    }
  }

  return(probs)
}

get_omp_threads <- function() {
  .Call("C_get_omp_threads", PACKAGE = "gllvm")
}

set_omp_threads <- function(t) {
  t <- as.integer(t)
  .Call("C_set_omp_threads", t, PACKAGE = "gllvm")
}

poisbinom <- function(prob) {
  .Call("C_poisbinom", prob, PACKAGE = "gllvm")
}


#' @title Dunn-Smyth residuals for species richness from a gllvm object
#' @description Calculates Dunn-Smyth residuals for species richness.
#'
#' @param object an object of class 'gllvm'.
#' @param predSR object returned by \code{"\link{predictSR.gllvm}"}
#' @param ... not used.
#'
#' @details
#' See \code{"\link{residuals.gllvm}"} for details.
#' 
#' @return A list of length 2 with 1) the fitted and 2) the dunn-smyth residuals.
#' 
#' @seealso \code{\link{predict.gllvm}}, \code{\link{predictSR.gllvm}}, \code{\link{residuals.gllvm}}
#' 
#' @author Bert van der Veen
#'
#'@author Bert van der Veen
#'
#'@aliases residuals.predictSR.gllvm
#'@method residuals predictSR.gllvm
#'@export
#'@export residuals.predictSR.gllvm
residuals.predictSR.gllvm <- function(predSR, object, ...){
  
  # ensure ordinal starts at 0
  if(any(object$family == "ordinal"))object$y[,object$family=="ordinal"] <- object$y[,object$family=="ordinal"] - apply(object$y[,object$family=="ordinal"],2,min, na.rm = TRUE)
  
  y = rowSums(ifelse(object$y[, predSR$spp]==0,0,1),na.rm=TRUE)
  n = length(y)
  pmf = predSR$predicted$fit
  
  K <- ncol(pmf) - 1
  
  cdf <- t(apply(pmf, 1, cumsum))
  
  b <- cdf[cbind(seq_len(n), pmin(y, K) + 1)]
  
  a <- numeric(n)
  nz <- y > 0 & !is.na(y)
  a[nz] <- cdf[cbind(which(nz), y[nz])]
  
  u <- a + (b - a) * runif(n)
  
  if (any(u>(1-1e-16), na.rm = TRUE))
    u[u>(1-1e-16)] <- 1 - 1e-16
  
  if (any(u == 0, na.rm = TRUE))
    u[u == 0] <- 1e-16
  
  list(fitted = predSR$expected$fit, residuals = qnorm(u))
}

#'@export
plot.predictSR.gllvm <- function(predSR, object, which = c(1,2), ...){
  res <- residuals.predictSR.gllvm(predSR, object)
  
  par(mfrow=c(1,length(which)))
  
  if(1%in%which)plot(xlab = "Expected species richness", ylab = "Dunn-Smyth residuals", x = res$fitted, y = res$residuals)
  if(2%in%which)qqnorm(res$residuals);qqline(res$residuals)
}

pbMode <- function(mu, n) {
  k <- floor(mu)
  
  lower <- k + 1/(k + 2)
  upper <- (k + 1) - 1/(n - k + 1)
  
  if (mu < lower) {
    return(k)
  } else if (mu > upper) {
    return(k + 1)
  } else {
    return(k)#c(k, k + 1))  # bimodal region
  }
}
