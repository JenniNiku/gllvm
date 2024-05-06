#' @title Predict Method for gllvm Fits
#' @description Obtains predictions from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param type the type of prediction required. The default (\code{"link"}) is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable. that is, the predictions for the binomial model are predicted probabilities. In case of ordinal data, \code{type = "response"} gives predicted probabilities for each level of ordinal variable.
#' @param newX A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.
#' @param newTR A new data frame of traits for each response taxon. If omitted, the original matrix of traits is used.
#' @param newLV A new matrix of latent variables.  If omitted, the original matrix of latent variables is used. Note that number of rows/sites must be the same for \code{newX} (if X covariates are included in the model).
#' @param level specification for how to predict. Level one (\code{level = 1}) attempts to use the predicted site scores from variational approximations or laplace approximation or given site scores in \code{newLV}. Level 0 sets the latent variable to zero. Defaults to 1.
#' @param offset specification whether of not offset values are included to the predictions in case they are in the model, defaults to \code{TRUE} when offset values that are used to fit the model are included to the predictions. Alternatives are matrix/vector (number of rows must match with the \code{newX}) of new offset values or \code{FALSE}, when offsets are ignored.
#' @param ... not used.
#'
#' @details
#' If \code{newX}, \code{newTR} and \code{newLV} are omitted the predictions are based on the data used for fitting the model. Notice that \code{newTR} need to match with the number of species in the original data.
#' Instead, new sites can be specified in \code{newX}. If predictors \code{newX} (and \code{newTR}) are given, and \code{newLV} is not, latent variables are not used in the predictions.
#' 
#' @return A matrix containing requested predictor types.
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>,  David Warton
#'
#' @examples
#' \donttest{
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- scale(antTraits$env[, 1:3])
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, family = poisson())
#'# fitted values
#'predfit <- predict(fit, type = "response")
#'
#'# linear predictors
#'predlin <- predict(fit)
#'# Predict new sites:
#'# Generate matrix of environmental variables for 10 new sites
#'xnew <- cbind(rnorm(10), rnorm(10), rnorm(10))
#'colnames(xnew) <- colnames(X)
#'predfit <- predict(fit, newX = xnew, type = "response", level = 0)
#'
#'TR <- (antTraits$tr[, 1:3])
#'fitt <- gllvm(y = y, X, TR, family = poisson())
#'# linear predictors
#'predlin <- predict(fitt)
#'# Predict new sites:
#'# Generate matrix of environmental variables for 10 new sites
#'xnew <- cbind(rnorm(10), rnorm(10), rnorm(10))
#'colnames(xnew) <- colnames(X)
#'# Generate matrix of traits for species
#'trnew <- data.frame(Femur.length = rnorm(41), No.spines = rnorm(41),
#'  Pilosity = factor(sample(0:3, 41, replace = TRUE)))
#'predfit <- predict(fitt, newX = xnew, newTR = trnew, type = "response", level = 0)
#'}
#'@aliases predict predict.gllvm
#'@method predict gllvm
#'
#'@export
#'@export predict.gllvm
predict.gllvm <- function(object, newX = NULL, newTR = NULL, newLV = NULL, type ="link", level = 1, offset = TRUE, ...){

  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  
  r0 <- NULL
  newdata <- newX
  p <- ncol(object$y)
  if(object$family == "betaH")p <- p*2
  n <- max(nrow(object$y), nrow(newdata), nrow(newLV))
  if (!is.null(newdata)) 
    n <- nrow(newdata)
  if (is.null(newdata) && !is.null(object$X) && !is.null(newLV)) {
    if (nrow(newLV) != nrow(object$y)) {
      stop("Number of rows in newLV must equal to the number of rows in the response matrix, if environmental variables are included in the model and newX is not included.")
    }
  }
  if (is.null(newdata) & !is.null(newLV)) {
    n <- nrow(newLV)
    if(!is.null(object$X) || !is.null(object$lv.X)){
      if((nrow(object$y)!=n)) 
        stop("Number of rows in 'newLV' must be the same as the number of rows in the data used for model fitting, if new X values for the corresponding units (in 'newLV') are not given.")
    }
  }
  if (!is.null(object$X)) {
    formula <- formula(terms(object))
  }
  else {
    formula <- NULL
  }
  if (object$row.eff != FALSE) {
    if(!is.null(newX) & length(all.vars(object$call$row.eff))>0){
      if(any(all.vars(object$call$row.eff) %in% colnames(newX))) {
        warning("Using row effects for predicting new sites does not work yet.")
      }
    }
  }
  b0 <- object$params$beta0
  eta <- matrix(b0, n, p, byrow = TRUE)
  if (!is.null(newTR)) 
    if (nrow(newTR) != p) 
      stop("Number of rows in newTR must match to the number of responses in the original data matrix.")
  if (is.null(colnames(object$y))) {
    colnames(object$y) <- paste("y", 1:p, sep = "")
  }
  if (!is.null(object$X) && is.null(object$TR)) {
    B <- object$params$Xcoef
    if (is.null(newdata)) {
      X.d <- Xnew <- object$X
    }
    else {
      X.d <- Xnew <- newdata
    }
    if (!is.null(newdata)) {
      if (is.null(object$call$formula)) {
        n1 <- colnames(newdata)
        formula1 <- paste("~ ", n1[1], sep = "")
        if (length(n1) > 1) {
          for (i1 in 2:length(n1)) {
            formula1 <- paste(formula1, n1[i1], sep = "+")
          }
        }
        formula1 <- paste(formula1, "1", sep = "-")
        formula1 <- formula(formula1)
        Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      }
      formula <- as.formula(nobars1_(object$call$formula)) # due to potential REs
      xb <- as.matrix(model.matrix(formula, data = data.frame(Xnew)))
      X.d <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)")), 
                          drop = F])
      colnames(X.d) <- colnames(xb)[!(colnames(xb) %in% 
                                        c("(Intercept)"))]
    }
    else {
      X.d <- object$X.design
    }
    eta <- eta + X.d %*% t(B)
  }
  if (!is.null(object$X) && !is.null(object$TR)) {
    B <- object$params$B
    X.d <- object$X.design
    if (!is.null(newdata) && !is.null(newTR)) {
      Xnew <- newdata
      TRnew <- newTR
      y1 <- object$y[sample(1:nrow(object$y), nrow(Xnew), 
                            replace = TRUE), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", 
                    varying = colnames(data.frame(y1)), v.names = "y", 
                    timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if (!is.null(newdata) && is.null(newTR)) {
      formula <- formula(object$formula)
      n1 <- colnames(newdata)
      formula1 <- paste("~", n1[1], sep = "")
      if (length(n1) > 1) {
        for (i1 in 2:length(n1)) {
          formula1 <- paste(formula1, n1[i1], sep = "+")
        }
      }
      formula1 <- paste(formula1, "1", sep = "-")
      formula1 <- formula(formula1)
      Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      TRnew <- object$TR
      y1 <- object$y[sample(1:nrow(object$y), nrow(Xnew), 
                            replace = TRUE), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", 
                    varying = colnames(data.frame(y1)), v.names = "y", 
                    timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if (is.null(newdata) && !is.null(newTR)) {
      if (nrow(newTR) != p) 
        stop("Number of rows in TRnew must equal to the number of response variables.")
      formula <- formula(object$formula)
      TRnew <- newTR
      n1 <- colnames(newTR)
      formula1 <- paste("~", n1[1], sep = "")
      if (length(n1) > 1) {
        for (i1 in 2:length(n1)) {
          formula1 <- paste(formula1, n1[i1], sep = "+")
        }
      }
      formula1 <- paste(formula1, "1", sep = "-")
      formula1 <- formula(formula1)
      TRnew <- as.matrix(model.matrix(formula1, data = data.frame(newTR)))
      Xnew <- object$X
      y1 <- object$y[sample(1:nrow(object$y), nrow(Xnew), 
                            replace = TRUE), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", 
                    varying = colnames(data.frame(y1)), v.names = "y", 
                    timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    X.d <- as.matrix(X.d[, colnames(X.d) != "(Intercept)"])
    eta <- eta + matrix(X.d %*% B, n, p)
    if (!is.null(object$randomX)) {
      if (!is.null(newdata)) {
        tr <- try(X.xr <- as.matrix(model.matrix(object$randomX, 
                                                 data = yXT)), silent = TRUE)
        if (inherits(tr, "try-error") | (object$randomX=="~.")) {
          X.xr <- as.matrix(yXT[, colnames(object$Xrandom)])
          colnames(X.xr) <- colnames(object$Xrandom)
        }
        rnam <- colnames(X.xr)[!(colnames(X.xr) %in% 
                                   c("(Intercept)"))]
        xr <- as.matrix(X.xr[1:n, !(colnames(X.xr) %in% 
                                      c("(Intercept)"))])
        if (NCOL(xr) == 1) 
          colnames(xr) <- rnam
      } else {
        xr <- object$Xrandom
      }
      eta <- eta + matrix(xr %*% object$params$Br, n, p)
    }
  }
  if (level == 1) {
    if (is.null(newLV) && !is.null(newdata) & (object$num.lv + 
                                               object$num.lv.c) > 0) {
      if (nrow(newdata) != nrow(object$y)) {
        stop("Level 1 predictions cannot be calculated for new X values if new latent variable values are not given. Change to 'level = 0' predictions.")
      }
    }
    if (object$num.lv > 0 | (object$num.lv.c + object$num.RR) > 
        0) {
      
      if(!is.null(object$lvs) && inherits(object$lvCor,"formula")){
        if(nrow(object$lvs)!=n) object$lvs = object$TMBfn$env$data$dLV%*%object$lvs # !!!
      }
      
      if (!is.null(newLV)) {
        if (ncol(newLV) != (object$num.lv + object$num.lv.c)) 
          stop("Number of latent variables in input doesn't equal to the number of latent variables in the model.")
        if (!is.null(newdata)) {
          if (nrow(newLV) != nrow(newdata)) 
            stop("Number of rows in newLV must equal to the number of rows in newX, if newX is included, otherwise same as number of rows in the response matrix.")
        }
        
        if (object$num.RR == 0) {
          lvs <- t(t(newLV) * object$params$sigma.lv)
        } else {
          if (object$num.lv.c > 0) {
            lvs <- cbind(t(t(newLV[, 1:object$num.lv.c]) * 
                             object$params$sigma.lv[1:object$num.lv.c]), 
                         matrix(0, ncol = object$num.RR, nrow = n), 
                         t(t(newLV[, -c(1:object$num.lv.c)]) * 
                             object$params$sigma.lv[1:object$num.lv]))
          }
          else if (object$num.lv > 0 & object$num.lv.c == 0) {
            lvs <- cbind(matrix(0, ncol = object$num.RR, 
                                nrow = n), t(t(newLV) * object$params$sigma.lv))
          } else {
            lvs <- matrix(0, ncol = object$num.RR, nrow = n)
          }
        }
        
      } else {
        if (object$num.RR == 0) {
          lvs <- t(t(object$lvs) * object$params$sigma.lv)
        } else {
          if (object$num.lv.c > 0) {
            lvs <- cbind(t(t(object$lvs[, 1:object$num.lv.c]) * 
                             object$params$sigma.lv[1:object$num.lv.c]), 
                         matrix(0, ncol = object$num.RR, nrow = n), 
                         t(t(object$lvs[, -c(1:object$num.lv.c)]) * 
                             object$params$sigma.lv[1:object$num.lv]))
          }
          else if (object$num.lv > 0 & object$num.lv.c == 
                   0) {
            lvs <- cbind(matrix(0, ncol = object$num.RR, 
                                nrow = n), t(t(object$lvs) * object$params$sigma.lv))
          }
          else {
            lvs <- matrix(0, ncol = object$num.RR, nrow = n)
          }
        }
      }
      if ((object$num.lv.c + object$num.RR) > 0 & !is.null(newdata)) {
        lv.X <- model.matrix(object$lv.formula, as.data.frame(newdata))[,-1, drop = F]
      } else {
        lv.X <- object$lv.X.design
      }
      theta <- (object$params$theta[, 1:(object$num.lv + 
                                          (object$num.lv.c + object$num.RR)), drop = F])
      eta <- eta + lvs %*% t(theta)
      if ((object$num.lv.c + object$num.RR) > 0) {
        eta <- eta + lv.X %*% object$params$LvXcoef %*% 
          t((theta[, 1:(object$num.lv.c + object$num.RR), drop = F]))
      }
      if (object$quadratic != FALSE) {
        if(object$num.lv>0){
          theta2 <- (object$params$theta[, -c(1:(object$num.lv.c + object$num.RR+object$num.lv)), drop = F])
          theta2 <- (theta2[, (object$num.lv.c + object$num.RR+1):ncol(theta2), drop = F])
          eta <- eta + (lvs[,(ncol(lvs)-object$num.lv+1):ncol(lvs), drop = F])^2 %*% t(theta2)
        }
        if ((object$num.lv.c + object$num.RR) > 0) {
          theta2 <- object$params$theta[,-c(1:(object$num.lv+object$num.lv.c+object$num.RR))]
          theta2C <- abs(theta2[, 1:(object$num.lv.c + 
                                       object$num.RR), drop = F])
          lvs <- lvs[,1:(object$num.lv.c+object$num.RR)] + lv.X%*%object$params$LvXcoef
          for (j in 1:p) {
            eta[, j] <- eta[, j] - lvs^2%*%theta2C[j, ]
          }
        }
      }
    }
  }
  
  
  if ((object$row.eff %in% c("random", "fixed", "TRUE")) && is.null(r0) & is.null(newX)) {
    if(!is.null(object$params$row.params) & !(object$row.eff %in% c("fixed", "TRUE"))){
     object$params$row.params = object$TMBfn$env$data$dr0%*%object$params$row.params # !!!
    } else {
      object$params$row.params = as.matrix(object$params$row.params)
    }
    r0 <- object$params$row.params
    if((object$row.eff %in% "random") && (level==0)) r0 = r0*0
    eta <- eta + r0%*%rep(1,ncol(object$y))
  }
  
  if (object$col.eff$col.eff == "random" && is.null(newX)) {
    eta <- eta + as.matrix(object$col.eff$spdr%*%object$params$Br)
    eta <- eta + as.matrix(object$col.eff$spdr[,names(object$params$B),drop=FALSE]%*%matrix(object$params$B, ncol = ncol(object$y), nrow = length(object$params$B)))
  }else if(object$col.eff$col.eff == "random" && !is.null(newX) && level == 1){
    stop("Prediction with column effects not yet implemented.")
      # bar.f <- findbars1(object$col.eff$col.eff.formula) # list with 3 terms
      # mf <- model.frame(subbars1(object$col.eff$col.eff.formula),data=data.frame(newX))
      # RElist <- mkReTrms1(bar.f,mf)
      # newspdr <- Matrix::t(RElist$Zt)
      # eta <- eta + as.matrix(newspdr%*%object$params$Br)
      # eta <- eta + model.matrix(as.formula(paste0("~", paste0(all.vars(object$col.eff$col.eff.formula), collapse = "+"))), newX)[,-1, drop = FALSE]%*%t(object$params$Xcoef[,grepl("RE_mean_",colnames(object$params$Xcoef)),drop=FALSE])
  }

  if(!is.null(object$offset)){
    if(offset!=FALSE){
      if(is.matrix(offset)){
        if((NROW(offset) == NROW(eta))){
          eta <- eta+object$offset
        } else {stop(paste("Incorrect dimension for the 'offset', number of rows should now be ", NROW(eta)))}
      } else if((NROW(object$offset) == NROW(eta))){
        eta <- eta+object$offset
        } else {warning(paste("Could not include offset values as 'object$offset' has incorrect dimension, set 'offset = FALSE' or include new offset values"))}
    }
  }
  
  if(object$family %in% c("poisson", "negative.binomial", "tweedie", "gamma", "exponential"))
    ilinkfun <- exp
  if (object$family == "binomial" || (object$family == "beta") || (object$family == "betaH") || (object$family == "orderedBeta")) 
    ilinkfun <- binomial(link = object$link)$linkinv
  if (object$family == "ordinal") 
    ilinkfun <- pnorm
  if (object$family %in% c("ZIP","ZINB")) 
    ilinkfun <- function(eta) exp(eta) * (1 - matrix(object$params$phi, 
                                                     n, p, byrow = TRUE))
  if (object$family == "gaussian") 
    ilinkfun <- gaussian()$linkinv
  out <- NULL
  preds <- NULL
  if ("link" %in% type) 
    out <- eta
  if ("response" %in% type) 
    out <- ilinkfun(eta)
  if (is.null(newdata) && is.null(newTR) && is.null(newLV) && 
      "logL" %in% type) 
    out <- object$logL
  if (object$family == "ordinal" && type == "response") {
    if (object$zeta.struc == "species") {
      k.max <- apply(object$params$zeta, 1, function(x) length(x[!is.na(x)])) + 
        1
      preds <- array(NA, dim = c(max(k.max), nrow(eta), 
                                 p), dimnames = list(paste("level", 1:max(k.max), 
                                                           sep = ""), NULL, NULL))
      for (i in 1:n) {
        for (j in 1:p) {
          probK <- NULL
          probK[1] <- pnorm(object$params$zeta[j, 1] - 
                              eta[i, j], log.p = FALSE)
          probK[k.max[j]] <- 1 - pnorm(object$params$zeta[j, 
                                                          k.max[j] - 1] - eta[i, j])
          if (k.max[j] > 2) {
            j.levels <- 2:(k.max[j] - 1)
            for (k in j.levels) {
              probK[k] <- pnorm(object$params$zeta[j, 
                                                   k] - eta[i, j]) - pnorm(object$params$zeta[j, 
                                                                                              k - 1] - eta[i, j])
            }
          }
          preds[, i, j] <- c(probK, rep(NA, max(k.max) - 
                                          k.max[j]))
        }
      }
      out <- preds
    }
    else {
      k.max <- length(object$params$zeta) + 1
      preds <- array(NA, dim = c(k.max, nrow(eta), p), 
                     dimnames = list(paste("level", 1:max(k.max), 
                                           sep = ""), NULL, NULL))
      for (i in 1:n) {
        for (j in 1:p) {
          probK <- NULL
          probK[1] <- pnorm(object$params$zeta[1] - eta[i, 
                                                        j], log.p = FALSE)
          probK[k.max] <- 1 - pnorm(object$params$zeta[k.max - 
                                                         1] - eta[i, j])
          levels <- 2:(k.max - 1)
          for (k in levels) {
            probK[k] <- pnorm(object$params$zeta[k] - 
                                eta[i, j]) - pnorm(object$params$zeta[k - 
                                                                        1] - eta[i, j])
          }
          preds[, i, j] <- c(probK)
        }
      }
      out <- preds
    }
    dimnames(preds)[[3]] <- colnames(object$y)
  }
  try(rownames(out) <- 1:NROW(out), silent = TRUE)
  if(any(class(out) %in% "dgeMatrix")) out <- as.matrix(out)
  return(out)
}
