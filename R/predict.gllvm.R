#' @title Predict Method for gllvm Fits
#' @description Obtains predictions from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param type the type of prediction required. The default (\code{"link"}) is on the scale of the linear predictors; the alternatives are \code{"response"} and \code{"class"} (predicted classes, works only for binomial model with 1 trial or ordinal model). \code{"response"} is on the scale of the response variable. That is, the predictions for the binomial model are predicted probabilities. In case of ordinal data, \code{type = "response"} gives predicted probabilities for each level of ordinal variable.
#' @param newX A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.
#' @param newTR A new data frame of traits for each response taxon. If omitted, the original matrix of traits is used.
#' @param newLV A new matrix of latent variables.  If omitted, the original matrix of latent variables is used. Note that number of rows/sites must be the same for \code{newX} (if X covariates are included in the model).
#' @param level specification for how to predict. Level one (\code{level = 1}) attempts to use the predicted site scores from variational approximations or laplace approximation or given site scores in \code{newLV}. Level 0 sets the latent variable to zero. Defaults to 1.
#' @param offset specification whether of not offset values are included to the predictions in case they are in the model, defaults to \code{TRUE} when offset values that are used to fit the model are included to the predictions. Alternatives are matrix/vector (number of rows must match with the \code{newX}) of new offset values or \code{FALSE}, when offsets are ignored.
#' @param se.fit logical. If \code{TRUE}, performs 1000 simulations from the asymptotic covariance matrix for fixed effects, and from the CMSEP covariance matrix for the random effects. If an integer is provided, it is used as the number of simulations instead. 
#' @param alpha numerical between 0 and 1, defaults to 0.95. The confidence level for se.fit.
#' @param seed numeric, defaults to 42. Seed used for simulation in se.fit.
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
#'data(antTraits, package = "mvabund")
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
predict.gllvm <- function(object, newX = NULL, newTR = NULL, newLV = NULL, type ="link", level = 1, offset = TRUE, se.fit = FALSE, alpha = 0.95, seed = 42, ...){

  # if((is.numeric(se.fit)||se.fit) && any(object$family == "ordinal"))warning("Confidence intervals for ordinal models can be sensitive to the number of simulations.")
  if(type=="class" & all(!(object$family %in% c("binomial", "ordinal")))) {
    stop("type='class' can be calculated only for ordinal or binary data.")
  }
  # backward compatibility
  if(is.null(object$params$row.params.random) && !inherits(object$row.eff, "formula") && object$row.eff == "random"){
    object$params$row.params.random <- object$params$row.params
    object$row.eff <- ~(1|site)
  }
  if(is.null(object$params$row.params.fixed) && !inherits(object$row.eff, "formula") && object$row.eff == "fixed"){
    object$params$row.params.fixed <- object$params$row.params
    object$xr <- diag(nrow(object$y))
  }
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  
  # end backward compatibility
  
  newdata <- newX
  p <- ncol(object$y)
  if(length(object$family) != p) object$family = rep(object$family,p) [1:p]
  
  if(any(object$family == "betaH"))p <- p+ sum(object$family == "betaH")
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
        # formula1 <- paste(formula1, "1", sep = "-")
        formula1 <- formula(formula1)
        Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
        X.d <- as.matrix(Xnew[, !(colnames(Xnew) %in% c("(Intercept)")), 
                            drop = F])
      } else {
        formula <- as.formula(nobars1_(object$call$formula)) # due to potential REs
        xb <- as.matrix(model.matrix(formula, data = data.frame(Xnew)))
        X.d <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)")), 
                            drop = F])
        colnames(X.d) <- colnames(xb)[!(colnames(xb) %in% 
                                          c("(Intercept)"))]
      }
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
        object$lvs = as.matrix(object$TMBfn$env$data$dLV%*%object$lvs) # !!!
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
                             object$params$sigma.lv[-(1:object$num.lv.c)]))
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
                             object$params$sigma.lv[-c(1:object$num.lv.c)]))
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
        if(anyBars(object$lv.formula)){
          bar.f <- findbars1(object$lv.formula) # list with 3 terms
          lv.X <- model.frame(subbars1(reformulate(sprintf("(%s)", sapply(findbars1(object$lv.formula), deparse1)))),data=as.data.frame(newdata))
          RElistLV <- mkReTrms1(bar.f,lv.X, nocorr=corstruc(expandDoubleVerts2(object$lv.formula))) #still add find double bars
          lv.X = t(as.matrix(RElistLV$Zt))
          print(lv.X)
        }else{
          lv.X <- model.matrix(object$lv.formula, as.data.frame(newdata))[,-1, drop = F]
        }
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
          # set theta2 as matrix to avoid dimensionality error
          theta2 <- as.matrix(object$params$theta[, -c(1:(object$num.lv.c + object$num.RR+object$num.lv)), drop = F])
          theta2 <- (theta2[, (object$num.lv.c + object$num.RR+1):ncol(theta2), drop = F])
          eta <- eta + (lvs[,(ncol(lvs)-object$num.lv+1):ncol(lvs), drop = F])^2 %*% t(theta2)
        }
        if ((object$num.lv.c + object$num.RR) > 0) {
          # set theta2 as matrix to avoid dimensionality error
          theta2 <- as.matrix(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c+object$num.RR))])
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
  
  r0 <- NULL
  if (inherits(object$row.eff, "formula") & is.null(newX)) {
    if(!is.null(object$params$row.params.random)){
     object$params$row.params.random = object$TMBfn$env$data$dr0%*%object$params$row.params.random # !!!
     if(level==0)object$params$row.params.random = object$params$row.params.random*0
     r0 <- cbind(r0, as.matrix(object$params$row.params.random))
    } 
    if(!is.null(object$params$row.params.fixed)){
        object$params$row.params.fixed = object$TMBfn$env$data$xr%*%as.matrix(object$params$row.params.fixed)
        r0 <- cbind(r0, as.matrix(object$params$row.params.fixed))
    }
    
    eta <- eta + as.matrix(rowSums(r0))%*%rep(1,p)
  }else if(inherits(object$row.eff, "formula") & !is.null(newX)){
      row.eff <- object$row.eff

      # code is an undressed version of the code in gllvm.R
      # first, random effects part
      if(anyBars(row.eff)){
        row.form <- allbars(row.eff)
        bar.f <- findbars1(row.form) # list with 3 terms
        grps <- unique(unlist(lapply(bar.f, all.vars)))
        form.parts <- strsplit(deparse1(row.form),split="\\+")[[1]]
        nested.parts <- grepl("/",form.parts)
        if(any(nested.parts)){
          form.parts <- strsplit(deparse1(row.form),split="\\+")[[1]]
          
          for(i in which(nested.parts))
            form.parts[i] <- paste0("(", findbars1(formula(paste0("~",form.parts[i]))),")",collapse="+")
          
          corstruc.form <- as.formula(paste0("~", paste0(form.parts,collapse="+")))
        }else{
          corstruc.form <- row.form
        }
        cstruc <- corstruc(corstruc.form)
        corWithin <- ifelse(cstruc %in% c("diag","ustruc"), FALSE, corWithin)
        
        if(!is.null(bar.f)) {
          mf <- model.frame(subbars1(row.form),data=newX)
          # adjust correlated terms for "corWithin = TRUE"; site-specific random effects with group-specific structure
          # consequence: we can use Zt everywhere
          mf.new <- mf
          if(any(corWithin)){
            mf.new[, corWithin] <- apply(mf[, corWithin, drop=F],2,function(x)order(order(x)))
          }
          colnames(mf.new) <- colnames(mf)
          RElistRow <- mkReTrms1(bar.f, mf.new, nocorr=cstruc)
          dr <- Matrix::t(RElistRow$Zt)
        }
        row.eff <- nobars1_(row.eff)
      }
      # second, fixed effects part
      if(inherits(row.eff, "formula") && length(all.vars(terms(row.eff)))>0){
        xr <- model.matrix(row.eff, newX)[,-1,drop=FALSE]
      }

    r0 <- NULL
    if (inherits(object$row.eff, "formula") & is.null(newX)) {
      if(!is.null(object$params$row.params.random)){
        object$params$row.params.random = dr%*%object$params$row.params.random # !!!
        if(level==0)object$params$row.params.random = object$params$row.params.random*0
        r0 <- cbind(r0, as.matrix(object$params$row.params.random))
      } 
      if(!is.null(object$params$row.params.fixed)){
        object$params$row.params.fixed = xr%*%as.matrix(object$params$row.params.fixed)
        r0 <- cbind(r0, as.matrix(object$params$row.params.fixed))
      }
      
      eta <- eta + as.matrix(rowSums(r0))%*%rep(1,p)
    }
  }
  if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE # backward compatibility
  
  if (object$col.eff$col.eff == "random" && is.null(newX) && is.null(object$TR)) {
    eta <- eta + as.matrix(object$col.eff$spdr%*%object$params$Br)
    if(!is.null(object$params[["B"]]))eta <- eta + as.matrix(object$col.eff$spdr[,names(object$params$B),drop=FALSE]%*%matrix(object$params$B, ncol = ncol(object$y), nrow = length(object$params$B)))
  }else if(object$col.eff$col.eff == "random" && !is.null(newX) && level == 1 && is.null(object$TR)){
    bar.f <- findbars1(object$col.eff$col.eff.formula) # list with 3 terms
    mf <- model.frame(subbars1(reformulate(sprintf("(%s)", sapply(findbars1(object$col.eff$col.eff.formula), deparse1)))),data=data.frame(newdata))
    
    if(length(bar.f)==1 & bar.f[[1]]==bquote(1|1)){
      X.col.eff <- mf <- data.frame(Intercept=rep(1,nrow(object$y)))
    }
    
    RElistSP<- mkReTrms1(bar.f, mf, nocorr=corstruc(expandDoubleVerts2(object$col.eff$col.eff.formula)))
    spdr <- Matrix::t(RElistSP$Zt)
    eta <- eta + as.matrix(spdr%*%object$params$Br)
    if(!is.null(object$params[["B"]]))eta <- eta + as.matrix(spdr[,names(object$params$B),drop=FALSE]%*%matrix(object$params$B, ncol = ncol(object$y), nrow = length(object$params$B)))
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
  
  ilinkfun <- list()
  pointer <- NULL
  if(any(object$family %in% c("poisson", "negative.binomial","negative.binomial1", "tweedie", "gamma", "exponential"))){
    ilinkfun <- c(ilinkfun, exp)
    pointer[object$family %in% c("poisson", "negative.binomial","negative.binomial1", "tweedie", "gamma", "exponential")] <- length(ilinkfun)
  }
  if (any(object$family %in% c("binomial", "beta", "betaH", "orderedBeta", "ordinal", "ZIB", "ZNIB"))){
    ilinkfun <- c(ilinkfun, binomial(link = object$link)$linkinv)
    pointer[object$family %in% c("binomial", "beta", "betaH", "orderedBeta", "ordinal", "ZIB", "ZNIB")] <- length(ilinkfun)
  }
  if (any(object$family %in% c("ZIP","ZINB"))) {
    # ilinkfun <- c(ilinkfun, function(eta) exp(eta) * (1 - matrix(object$params$phi, n, p, byrow = TRUE)))
    ilinkfun <- c(ilinkfun, function(eta) exp(eta))
    pointer[object$family %in% c("ZIP","ZINB")] <- length(ilinkfun)
  }
  if (any(object$family == "gaussian")) {
    ilinkfun <- c(ilinkfun, gaussian()$linkinv)
    pointer[object$family %in% c("gaussian")] <- length(ilinkfun)
  }
  out <- NULL
  preds <- NULL
  if ("link" %in% type) 
    out <- eta
  if ("response" %in% type) {
    out <- sapply(1:p, function(j) ilinkfun[[pointer[j]]](eta[,j]))
    if(any(object$family %in% c("ZIP","ZINB"))) 
      out[,object$family %in% c("ZIP","ZINB")] <- out[,object$family %in% c("ZIP","ZINB")]*(1 - matrix(object$params$phi[object$family %in% c("ZIP","ZINB")], n, sum(object$family %in% c("ZIP","ZINB")), byrow = TRUE))
  }
  if ("class" %in% type & any(object$family == "binomial")) {
    if(is.null(out)){ out <- eta; out[] <- NA}
    out[,object$family == "binomial"] <- round(sapply((1:p)[object$family == "binomial"], function(j) ilinkfun[[pointer[j]]](eta[,j])))
  }
  if (is.null(newdata) && is.null(newTR) && is.null(newLV) && 
      "logL" %in% type) 
    out <- object$logL
  if (any(object$family == "ordinal") && (type == "response" | type == "class")) {
    ordi_ind <- c(1:p)[object$family == "ordinal"]
    if (object$zeta.struc == "species") {
      k.max <- apply(object$params$zeta, 1, function(x) length(x[!is.na(x)])) + 1
      preds <- array(NA, dim = c(max(k.max, na.rm = TRUE), nrow(eta), 
                                 p), dimnames = list(paste("level", 1:max(k.max, na.rm = TRUE), 
                                                           sep = ""), NULL, NULL))
      if(!all(object$family == "orderedBeta")) preds[1,,] <- out
        for (j in ordi_ind) {
          probK <- matrix(nrow=k.max[j], ncol = n)
          probK[1,] <- ilinkfun[[pointer[j]]](object$params$zeta[j, 1] -  eta[, j])
          probK[k.max[j],] <- 1 - ilinkfun[[pointer[j]]](object$params$zeta[j, k.max[j] - 1] - eta[, j])
          if (k.max[j] > 2) {
            j.levels <- 2:(k.max[j] - 1)
            for (k in j.levels) {
              probK[k,] <- ilinkfun[[pointer[j]]](object$params$zeta[j, k] - eta[, j]) - ilinkfun[[pointer[j]]](object$params$zeta[j, k - 1] - eta[, j])
            }
          }
          preds[1:k.max[j], , j] <- probK
        }
    } else {
      kz <- any(object$family == "orderedBeta")*2
      k.max <- length(object$params$zeta) + 1 - kz
      preds <- array(NA, dim = c(k.max, nrow(eta), p), 
                     dimnames = list(paste("level", 1:max(k.max), 
                                           sep = ""), NULL, NULL))
      if(!all(object$family == "orderedBeta")) preds[1,,] <- out
        for (j in 1:p) {
          probK <- matrix(nrow=k.max,ncol=n)
          probK[1,] <- ilinkfun[[pointer[j]]](object$params$zeta[1+kz] - eta[, j])
          probK[k.max,] <- 1 - ilinkfun[[pointer[j]]](object$params$zeta[k.max - 1 + kz] - eta[, j])
          levels <- 2:(k.max - 1)
          for (k in levels) {
            probK[k,] <- ilinkfun[[pointer[j]]](object$params$zeta[k + kz] - eta[, j]) - ilinkfun[[pointer[j]]](object$params$zeta[k - 1 + kz] - eta[, j])
          }
          preds[1:k.max,, j] <- probK
        }
      dimnames(preds)[[3]] <- colnames(object$y)
    }
    if(type == "response") {
      out <- preds
    }
    if(type == "class") {
      pred_class <- matrix(NA, dim(preds)[2], dim(preds)[3])
      for (j in ordi_ind) {
        for (i in 1:nrow(pred_class)) {
          pred_class[i,j] <- (order(preds[,i,j], decreasing = TRUE)[1]-1)
        }
      }
      colnames(pred_class) <- colnames(object$y)
      if(is.null(out)){ out <- eta; out[] <- NA}
      out[,object$family == "ordinal"] <- pred_class[,object$family == "ordinal"]
    }
  }
  try(rownames(out) <- 1:NROW(out), silent = TRUE)
  if(any(class(out) %in% "dgeMatrix")) out <- as.matrix(out)
  
  if((is.numeric(se.fit) || se.fit) && isFALSE(object$sd)){
    stop("Cannot calculate standard errors for the prediction without standard errors of the parameter estimates. Please refit the model with 'sd.erors = TRUE', or use 'se.gllvm'.")
  }else if(is.numeric(se.fit) || se.fit){
    num.lv <- object$num.lv
    num.RR <- object$num.RR
    num.lv.c <- object$num.lv.c
    num.lv.cor <- object$num.lvcor
    quadratic = object$quadratic
    
    incl <- object$Hess$incl
    incla <- object$Hess$incla
    if(object$method %in% c("VA","EVA")){
      incla <- rep(FALSE, length(object$TMBfn$par))
      if((num.lv.c+num.lv+num.lv.cor)>0)incla[names(object$TMBfn$par)=="u"] <- TRUE
      if("Br" %in% names(object$TMBfn$par))incla[names(object$TMBfn$par)=="Br"] <- TRUE
      if(num.RR>0 && !isFALSE(object$randomB))incla[names(object$TMBfn$par)=="b_lv"] <- TRUE
      if("r0r" %in% names(object$TMBfn$par))incla[names(object$TMBfn$par)=="r0r"] <- TRUE
    }
    
    if(is.null(incla)) incla <- 0
    
    R <- 1e3
    if(is.numeric(se.fit))R <- se.fit
    
    set.seed(seed)
    
    if(sum(incl)>0){
      Vf <- vcov(object)
      if(min(diag(Vf))/max(diag(Vf))<.Machine$double.eps)warning("Seeing some odd things in the fixed effects covariance matrix (variances of effects are not on the same scale), uncertainties may be inaccurate.\n")
      ffs <- try(MASS::mvrnorm(R, object$TMBfn$par[incl],Vf), silent = TRUE)
      if(inherits(ffs, "try-error"))stop("Covariance matrix of fixed effects is not semi positive-definite.")
      colnames(ffs) <- names(object$TMBfn$par)[incl]
      
      # Same code for organsing map as in se.gllvm
      if(any(colnames(ffs)%in%names(object$TMBfn$env$map))){
        map <- object$TMBfn$env$map[names(object$TMBfn$env$map)%in%colnames(ffs)]
        # rebuild se matrix if mapped parameters
        ffs.new <- NULL
        for(nm in unique(colnames(ffs))){
          if(!nm%in%names(map)){
            ffs.new <- cbind(ffs.new,ffs[,colnames(ffs)==nm, drop=FALSE])
          }else{
            ffs.new <- cbind(ffs.new, ffs[,colnames(ffs)==nm, drop=FALSE][,map[[nm]],drop=FALSE])
          }
          
        }
        ffs <- ffs.new
        rm(ffs.new)
      }
    
    }
    
    if(sum(incla)>0 && level>0){
      if(object$method == "LA"){
        Vr <- sdrandom(object$TMBfn, Vf, object$Hess$incl, return.covb = TRUE)
      }else{
        Vr <- CMSEPf(object, return.covb = TRUE)
      }
      if(min(diag(Vr))/max(diag(Vr))<.Machine$double.eps)warning("Seeing some odd things in the random effects (asymptotic) covariance matrix (variances of effects are not on the same scale), uncertainties may be inaccurate.\n")
      rfs <- try(MASS::mvrnorm(R, object$TMBfn$env$last.par.best[incla], Vr), silent = TRUE)
      if(inherits(rfs, "try-error"))stop("Covariance matrix of random effects is not semi positive-definite.")
    }
    
    # set all parameters to NA to avoid potential issues if we miss something
    object$params <- lapply(object$params, function(x) {
      if (is.numeric(x)) {
        x[] <- NA 
      }
      x
    })
    object$lvs[] <- NA
    
    predSims <- array(0.0, dim = c(R, dim(out)))
    newobject <- object
    
    # parameter organisation, partly from gllvm.TMB and traitTMB.R, should really store this in a separate function so we don't keep copying the same code...
    for(r in 1:R){
      # fixed effects first
      # fixed effects needed for prediction; beta0, Xcoef, sigmaLV, lambda, phi, zeta, 4th corner, r0f, b_lv, 
      
      if(sum(incl)>0){
      newpars <- relist_gllvm(ffs[r,], object$TMBfn$env$parList())
      newobject$params$beta0 <- newpars$b[1,]
      if(nrow(newpars$b)>1)newobject$params$Xcoef <- t(newpars$b[-1,,drop=FALSE])
      
      # LV stuff
      if((num.lv+num.lv.c)>0)newobject$params$sigma.lv <- abs(newpars$sigmaLV)
      
      # Organising theta is a pain, copied almost completely from gllvm.TMB.R
      if((num.lv+num.lv.c+num.RR)){
        theta <- matrix(0,p,num.lv+num.lv.c+num.RR)  
        if((num.lv.c+num.RR)>1){diag(theta[,1:(num.lv.c+num.RR)])<-1}else if((num.lv.c+num.RR)==1){theta[1,1]<-1}
        if(num.lv>1){diag(theta[,((num.lv.c+num.RR)+1):((num.lv.c+num.RR)+num.lv)])<-1}else if(num.lv==1){theta[1,((num.lv.c+num.RR)+1):((num.lv.c+num.RR)+num.lv)]<-1}
        if(num.lv>0&(num.lv.c+num.RR)==0){
          
          if(p>1) {
            theta[lower.tri(theta[,1:num.lv,drop=F],diag=FALSE)] <- newpars$lambda;
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(newpars$lambda2),ncol=num.lv,nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(newpars$lambda2))}  
          }
        }else if(num.lv==0&(num.lv.c+num.RR)>0){
          if(p>1) {
            theta[lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE)] <- newpars$lambda;
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(newpars$lambda2),ncol=(num.lv.c+num.RR),nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(newpars$lambda2))}  
          }
        }else if(num.lv>0&(num.lv.c+num.RR)>0){
          if(p>1) {
            theta[,1:(num.lv.c+num.RR)][lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE)] <- newpars$lambda[1:sum(lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE))];
            theta[,((num.lv.c+num.RR)+1):ncol(theta)][lower.tri(theta[,((num.lv.c+num.RR)+1):ncol(theta),drop=F],diag=FALSE)] <- newpars$lambda[(sum(lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE))+1):length(newpars$lambda)];
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(newpars$lambda2),ncol=num.lv+(num.lv.c+num.RR),nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(newpars$lambda2))}  
          }
        }
        newobject$params$theta <- theta
        
      }
    
      # Organising phi. Only necessary when predicing ZIP/ZINB stuff
      if(any(object$family %in% c("ZIP","ZINB"))){
        disp.group <- object$disp.group
        lp0 <- newpars$lg_phi[disp.group]
        object$params$phi <- (exp(lp0)/(1+exp(lp0)))[object$family %in% c("ZIP","ZINB")];
      }
      
      if(any(object$family %in% c("orderedBeta","ordinal")) && type == "response"){
        
        zetaO = NULL
          if(any(object$family%in%c("ordinal"))){
            K = max(object$TMBfn$env$data$y)-min(object$TMBfn$env$data$y)
          } else {K=2}
        
          if(object$zeta.struc =="common") {
            if(any(object$family%in%c("orderedBeta"))){
              zetaO <- c(zetaO, rep(TRUE,2))
             }
            if(any(object$family%in%c("ordinal"))){
              zetaO <- c(zetaO, rep(FALSE,(K-1)))
            }
          } else if(object$zeta.struc =="species") {
            o_ind <- c(1:ncol(object$y))[object$family%in%c("ordinal", "orderedBeta")]
            for (j in o_ind) {
              if(object$family[j]=="ordinal"){
                zetaO <- c(zetaO, rep(FALSE,length(na.omit(object$params$zeta[j,-1]))))
              } else {
                zetaO <- c(zetaO, rep(TRUE,2))
              }
            }
          }
 
        zetas <- newpars$zeta

        if(object$zeta.struc=="species"){
          zetanew <- matrix(NA,nrow=p,ncol=K)
          zetanew[,1] <- 0 
          idx<-0
          o_ind <- c(1:p)[object$family %in%c("ordinal","orderedBeta")]
          for(j in o_ind){
            if(object$family[j] == "ordinal"){
            k<-max(object$y[,j])-2
            if(k>0){
              for(l in 1:k){
                zetanew[j,l+1]<-zetas[idx+l]
              } 
            }
            idx<-idx+k
            zetanew[j,] <- cumsum(abs(zetanew[j,]))
          }else{
            zetanew[j,] <- c(zetas[idx +1], exp(zetas[idx +2]))
            idx<-idx+2
          }
          }
        }else{
          zetanew <- NULL
          if(any(object$family == "orderedBeta")){
            zetanew <- c(zetanew, zetas[1], exp(zetas[2]))
            names(zetanew) <- c("cutoff0","cutoff1")
          }
          if(any(object$family%in%c("ordinal"))){
            zetanew <- c(zetanew, 0,cumsum(abs(zetas[!zetaO])))
          }
        }
        newobject$params$zeta <- zetanew
      }
      # if(object$family == "orderedBeta"){
      #   zetas <- matrix(newpars$zeta,p,2)
      #   if(any(is.na(object$TMBfn$env$map$zeta))) zetas[is.na(object$env$TMBfn$map$zeta)] = attr(object$TMBfn$env$parameters$zeta, "shape")[is.na(object$TMBfn$env$map$zeta)]
      #   zetas[,2] = exp(zetas[,2])
      #   newobject$params$zeta <- zetas
      # }
      
      if(!is.null(newpars$B)){
        newobject$params$B <- newpars$B
        newobject$params$B <- newobject$params$B[!is.na(newpars$B)]
        names(newobject$params$B) <- names(object$params$B)
      }
      
      if(!is.null(newpars$params$r0f))
        newobject$params$row.params.fixed <- newpars$r0f
      
      if(num.RR>0 && isFALSE(object$randomB))
        newobject$params$LvXcoef <- newpars$b_lv
      }
      
      # random effects next
      # random effects needed for prediction; Br, r0r, u, b_lv, 
      if(sum(incla)>0 && level>0){
        newrfs <- relist_gllvm(rfs[r,], object$TMBfn$env$parList())
        
        if(!is.null(newrfs$Br))
          newobject$params$Br = newrfs$Br
        if(!is.null(newrfs$u))
          newobject$lvs <- newrfs$u
        if(!is.null(newrfs$b_lv))
          newobject$params$LvXcoef <- newrfs$b_lv
        if(!is.null(newrfs$r0r))
          newobject$params$row.params.random  <- newrfs$r0r
      }
      if(!any(object$family == "ordinal") || type == "link")
        predSims[r,,] <- predict(newobject, newX = newX, newTR = newTR, newLV = newLV, type = type, level = level, offset = offset, se.fit = FALSE)
      else
        predSims[r,,,] <- predict(newobject, newX = newX, newTR = newTR, newLV = newLV, type = type, level = level, offset = offset, se.fit = FALSE)
    }
    
    if(!any(object$family == "ordinal") || !is.null(alpha) && type == "link" || is.null(alpha)){
    if(is.null(alpha)){
    out <- list(fit = out, ci.sim = predSims) # so that we can access this for predictSR and predictPairwise
    }else if(!is.null(alpha)){
    ci <- apply(predSims, 2:3, quantile, prob = c((1-alpha)/2, 1-(1-alpha)/2))
    out <- list(fit = out, lower = ci[1,,], upper = ci[2,,])
    }
    }else if(any(object$family == "ordinal") && type == "response" && !is.null(alpha)){
      ci <- array(dim=c(2,nrow(out), n,p))
      for(j in 1:p){
        if(object$family[j] != "ordinal"){
          ci[,1,,j] <- apply(predSims[,,j], 2, prob = c((1-alpha)/2, 1-(1-alpha)/2)) 
        }else{
          # here we need to be a bit more careful because of the simplex constraint
          # we take the same approach as in predictSR.gllvm.R
          # we measure the distance to the center of the simplex (which we will consider the prediction on the point estimates)
          # and threshold over the distances
          for(i in 1:n){
            center = pmin(pmax(out[,i,j], .Machine$double.eps), 1-.Machine$double.eps)
            dists = suppressWarnings(hilbert_to_provided_center(pmin(pmax(predSims[,,i,j], .Machine$double.eps),1-.Machine$double.eps), center))
            threshold = quantile(dists, alpha, na.rm = TRUE)
            ci[,,i,j] <- suppressWarnings(apply(predSims[,,i,j][dists <= threshold, ],2,range, na.rm = TRUE))
          }
          if(any(!is.finite(ci)))ci[!is.finite(ci)] <- NA # guaranteed to happen on zeta.struc = "species"
        }
      }
      out <- list(fit = out, lower = ci[1,,,], upper = ci[2,,,])
    }  
  }
  
  return(out)
}
