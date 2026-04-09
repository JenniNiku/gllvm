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
#' @param ordinal.cat integer or \code{NULL}, defaults to \code{NULL}. For ordinal models with \code{type = "response"}, selects a single category to return instead of the full K×n×p array. When set, the return value is an n×p matrix of probabilities for that category. Primarily used internally by \code{predictSR} and \code{goodnessOfFit} to reduce memory use.
#' @param spp integer vector or \code{NULL}, defaults to \code{NULL}. If provided, predictions are computed only for the selected species (column indices into the response matrix), returning an \eqn{n \times \text{length(spp)}} matrix. All parameter-level computations (\code{beta0}, \code{Xcoef}, \code{theta}, etc.) are subsetted accordingly, so no unnecessary computation is performed for the remaining species.
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
predict.gllvm <- function(object, newX = NULL, newTR = NULL, newLV = NULL, type ="link", level = 1, offset = TRUE, se.fit = FALSE, alpha = 0.95, seed = 42, ordinal.cat = NULL, spp = NULL, ...){

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

  spp_idx <- if(!is.null(spp)) spp else seq_len(p)
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

  b0_full <- object$params$beta0
  b0 <- b0_full[spp_idx]
  eta <- matrix(b0, n, length(spp_idx), byrow = TRUE)
  if (!is.null(newTR)) 
    if (nrow(newTR) != p) 
      stop("Number of rows in newTR must match to the number of responses in the original data matrix.")
  if (is.null(colnames(object$y))) {
    colnames(object$y) <- paste("y", 1:p, sep = "")
  }
  if (!is.null(object$X) && is.null(object$TR)) {
    B <- object$params$Xcoef[spp_idx, , drop = FALSE]
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
    eta <- matrix(b0, n, length(spp_idx), byrow = TRUE) + matrix(X.d %*% B, n, p)[, spp_idx, drop = FALSE]
    
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
      eta <- eta + xr %*% object$params$Br[, spp_idx, drop = FALSE]
    }
  }
  if (level == 1 || (level == 0 && (object$num.lv.c + object$num.RR)>0)) {
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
        if ((object$num.lv+object$num.lv.c)>0 && ncol(newLV) != (object$num.lv + object$num.lv.c)) 
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
          RElistLV <- mkReTrms1(bar.f,lv.X, nocorr=corstruc(expandDoubleVerts2(object$lv.formula)), drop.unused.levels = FALSE) #still add find double bars
          # double check column names, because we may now have unobserved combinations of random effect levels in the matrix
          lv.X = t(as.matrix(RElistLV$Zt))
          lv.X <- lv.X[,colnames(lv.X)%in%colnames(object$lv.X.design),drop=FALSE]
          
        }else{
          lv.X <- model.matrix(object$lv.formula, as.data.frame(newdata))[,-1, drop = F]
        }
      } else {
        lv.X <- object$lv.X.design
      }
      theta <- (object$params$theta[spp_idx, 1:(object$num.lv +
                                          (object$num.lv.c + object$num.RR)), drop = F])
      eta <- eta + lvs %*% t(theta)
      if ((object$num.lv.c + object$num.RR) > 0) {
        eta <- eta + lv.X %*% object$params$LvXcoef %*% 
          t((theta[, 1:(object$num.lv.c + object$num.RR), drop = F]))
      }
      if (object$quadratic != FALSE) {
        if(object$num.lv>0){
          # set theta2 as matrix to avoid dimensionality error
          theta2 <- as.matrix(object$params$theta[spp_idx, -c(1:(object$num.lv.c + object$num.RR+object$num.lv)), drop = F])
          theta2 <- (theta2[, (object$num.lv.c + object$num.RR+1):ncol(theta2), drop = F])
          eta <- eta + (lvs[,(ncol(lvs)-object$num.lv+1):ncol(lvs), drop = F])^2 %*% t(theta2)
        }
        if ((object$num.lv.c + object$num.RR) > 0) {
          # set theta2 as matrix to avoid dimensionality error
          theta2 <- as.matrix(object$params$theta[spp_idx, -c(1:(object$num.lv+object$num.lv.c+object$num.RR))])
          theta2C <- abs(theta2[, 1:(object$num.lv.c +
                                       object$num.RR), drop = F])
          lvs <- lvs[,1:(object$num.lv.c+object$num.RR)] + lv.X%*%object$params$LvXcoef
          for (iter in seq_len(length(spp_idx))) {
            eta[, iter] <- eta[, iter] - lvs^2%*%theta2C[iter, ]
          }
        }
      }
    }
  }
  
  if(!is.null(object$params$row.params.fixed)||(!is.null(object$params$row.params.random)&&level>0)){
    if (inherits(object$row.eff, "formula") & is.null(newX)) {
      if(!is.null(object$params$row.params.random))dr = object$TMBfn$env$data$dr0
      if(!is.null(object$params$row.params.fixed))xr = object$TMBfn$env$data$xr
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
          RElistRow <- mkReTrms1(bar.f, mf.new, nocorr=cstruc, drop.unused.levels = FALSE)
          dr <- Matrix::t(RElistRow$Zt)
          # double check column names, because we may now have unobserved combinations of random effect levels in the matrix
          dr <- dr[,colnames(dr)%in%colnames(object$dr),drop=FALSE]
        }
        row.eff <- nobars1_(row.eff)
      }
      # second, fixed effects part
      if(inherits(row.eff, "formula") && length(all.vars(terms(row.eff)))>0){
        xr <- model.matrix(row.eff, newX)[,-1,drop=FALSE]
        xr <- xr[, colnames(xr) %in% names(object$params$row.params.fixed), drop = FALSE]
      }
    }
    
    r0 <- NULL
      if(!is.null(object$params$row.params.random) && level>0){
        object$params$row.params.random = dr%*%object$params$row.params.random # !!!
        r0 <- cbind(r0, as.matrix(object$params$row.params.random))
      } 
      if(!is.null(object$params$row.params.fixed)){
        object$params$row.params.fixed = xr%*%as.matrix(object$params$row.params.fixed)
        r0 <- cbind(r0, as.matrix(object$params$row.params.fixed))
      }
      
    eta <- eta + as.matrix(rowSums(r0))%*%rep(1, length(spp_idx))
    
  }
  
  if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE # backward compatibility
  
  if (object$col.eff$col.eff == "random" && is.null(newX) && is.null(object$TR)) {
    eta <- eta + as.matrix(object$col.eff$spdr%*%object$params$Br[, spp_idx, drop = FALSE])
    if(!is.null(object$params[["B"]]) && length(object$params[["B"]]>0))eta <- eta + as.matrix(object$col.eff$spdr[,names(object$params$B),drop=FALSE]%*%matrix(object$params$B, ncol = length(spp_idx), nrow = length(object$params$B)))
  }else if(object$col.eff$col.eff == "random" && !is.null(newX) && level == 1 && is.null(object$TR)){
    bar.f <- findbars1(object$col.eff$col.eff.formula) # list with 3 terms
    mf <- model.frame(subbars1(reformulate(sprintf("(%s)", sapply(findbars1(object$col.eff$col.eff.formula), deparse1)))),data=data.frame(newdata))
    
    if(length(bar.f)==1 & bar.f[[1]]==bquote(1|1)){
      X.col.eff <- mf <- data.frame(Intercept=rep(1,nrow(object$y)))
    }
    
    RElistSP<- mkReTrms1(bar.f, mf, nocorr=corstruc(expandDoubleVerts2(object$col.eff$col.eff.formula)), drop.unused.levels = FALSE)
    spdr <- Matrix::t(RElistSP$Zt)
    # double check column names, because we may now have unobserved combinations of random effect levels in the matrix
    spdr <- spdr[,colnames(spdr)%in%colnames(object$col.eff$spdr),drop=FALSE]
    
    eta <- eta + as.matrix(spdr%*%object$params$Br[, spp_idx, drop = FALSE])
    if(!is.null(object$params[["B"]]) && length(object$params[["B"]]>0))eta <- eta + as.matrix(spdr[,names(object$params$B),drop=FALSE]%*%matrix(object$params$B, ncol = length(spp_idx), nrow = length(object$params$B)))
  }

  if(!is.null(object$offset)){
    if(offset!=FALSE){
      if(is.matrix(offset)){
        if((NROW(offset) == NROW(eta))){
          eta <- eta+object$offset[, spp_idx, drop = FALSE]
        } else {stop(paste("Incorrect dimension for the 'offset', number of rows should now be ", NROW(eta)))}
      } else if((NROW(object$offset) == NROW(eta))){
        eta <- eta+object$offset[, spp_idx, drop = FALSE]
        } else {warning(paste("Could not include offset values as 'object$offset' has incorrect dimension, set 'offset = FALSE' or include new offset values"))}
    }
  }

  fam_pred <- object$family[spp_idx]
  famgroup1 <- fam_pred %in% c("poisson", "negative.binomial","negative.binomial1", "tweedie", "gamma", "exponential")
  famgroup2 <- fam_pred %in% c("binomial", "beta", "betaH", "orderedBeta", "ordinal", "ZIB", "ZNIB", "beta.binomial")
  famgroup3 <- fam_pred %in% c("ZIP","ZINB")
  famgroup4 <- fam_pred == "gaussian"
  famgroup5 <- fam_pred == "betaH"
  ilinkfun <- vector("list", any(famgroup1)+any(famgroup2)+any(famgroup3)+any(famgroup4)+any(famgroup5))
  pointer  <- NULL
  if(any(famgroup1)){
    ilinkfun[[1]] <- exp
    pointer[famgroup1] <- 1
  }
  if (any(famgroup2)){
    ilinkfun[[any(famgroup1)+1]] <- binomial(link = object$link)$linkinv
    pointer[famgroup2] <- any(famgroup1)+1
  }
  if (any(famgroup3)) {
    ilinkfun[[any(famgroup1)+any(famgroup2)+1]] <- exp
    pointer[famgroup3] <- any(famgroup1)+any(famgroup2)+1
  }
  if (any(famgroup4)) {
    ilinkfun[[any(famgroup1)+any(famgroup2)+any(famgroup3)+1]] <- identity
    pointer[famgroup4] <- any(famgroup1)+any(famgroup2)+any(famgroup3)+1
  }
  if(any(famgroup5)){
    pointer <- c(pointer, rep(unique(pointer[famgroup5]), sum(famgroup5)))
  }
  out <- NULL
  preds <- NULL
  if ("link" %in% type)
    out <- eta
  if ("response" %in% type) {
    out <- matrix(NA_real_, nrow = nrow(eta), ncol = length(spp_idx))
    for (j in unique(pointer)) {
      cols <- which(pointer == j)
      out[, cols] <- ilinkfun[[j]](eta[, cols, drop = FALSE])
    }
    if(any(famgroup3))
      out[, famgroup3] <- out[, famgroup3] * (1 - matrix(object$params$phi[spp_idx][famgroup3], n, sum(famgroup3), byrow = TRUE))
  }
  if ("class" %in% type & any(fam_pred == "binomial")) {
    if(is.null(out)){ out <- eta; out[] <- NA}
    out[, fam_pred == "binomial"] <- round(sapply(which(fam_pred == "binomial"), function(iter) ilinkfun[[pointer[iter]]](eta[, iter])))
  }
  if (is.null(newdata) && is.null(newTR) && is.null(newLV) &&
      "logL" %in% type)
    out <- object$logL
  if (any(fam_pred == "ordinal") && (type == "response" | type == "class")) {
    if(is.null(out)){ out <- eta; out[] <- NA}

    ordi_ind <- which(fam_pred == "ordinal")   # column indices into eta/out (1-based within spp_idx)
    if (object$zeta.struc == "species") {
      k.max <- apply(object$params$zeta[spp_idx, , drop = FALSE], 1, function(x) length(x[!is.na(x)])) + 1
      preds <- array(NA, dim = c(max(k.max, na.rm = TRUE), nrow(eta),
                                 length(spp_idx)), dimnames = list(paste("level", 1:max(k.max, na.rm = TRUE),
                                                           sep = ""), NULL, NULL))
      if(!all(fam_pred == "orderedBeta")) preds[1,,] <- out
      for (iter in ordi_ind) {
        j_full <- spp_idx[iter]  # index into object$params$zeta and object$family
        probK <- matrix(nrow=k.max[iter],ncol=n)
        probK[1:(k.max[iter]-1),] <- ilinkfun[[pointer[iter]]](outer(object$params$zeta[j_full,1:k.max[iter]-1], eta[,iter], function(zeta, eta)zeta-eta))
        if(k.max[iter]>2){
        probK[2:(k.max[iter]-1), ] <- probK[2:(k.max[iter]-1),,drop=FALSE] - probK[1:(k.max[iter]-2),,drop=FALSE]
        }
        probK[k.max[iter],] <- 1 - ilinkfun[[pointer[iter]]](object$params$zeta[j_full,k.max[iter] - 1] - eta[, iter])
        preds[1:k.max[iter],, iter] <- probK
      }
    } else {
      kz <- any(fam_pred == "orderedBeta")*2
      k.max <- length(object$params$zeta) + 1 - kz
      preds <- array(NA, dim = c(k.max, nrow(eta), length(spp_idx)),
                     dimnames = list(paste("level", 1:max(k.max),
                                           sep = ""), NULL, NULL))
      if(!all(fam_pred == "orderedBeta")) preds[1,,] <- out
        for (iter in ordi_ind) {
          probK <- matrix(nrow=k.max,ncol=n)
          probK[1:(k.max-1),] <- ilinkfun[[pointer[iter]]](outer(tail(object$params$zeta,k.max-1), eta[,iter], function(zeta, eta)zeta-eta))
          if(k.max>2){
          probK[2:(k.max-1), ] <- probK[2:(k.max-1),,drop=FALSE] - probK[1:(k.max-2),,drop=FALSE]
          }
          probK[k.max,] <- 1 - ilinkfun[[pointer[iter]]](object$params$zeta[k.max - 1 + kz] - eta[, iter])
          preds[1:k.max,, iter] <- probK
        }
      dimnames(preds)[[3]] <- colnames(object$y)[spp_idx]
    }
    if(type == "response") {
      out <- preds
      if(!is.null(ordinal.cat)) out <- preds[ordinal.cat,,]
    }
    if(type == "class") {
      pred_class <- matrix(NA, dim(preds)[2], dim(preds)[3])
      for (iter in ordi_ind) {
        for (i in 1:nrow(pred_class)) {
          pred_class[i, iter] <- (order(preds[,i,iter], decreasing = TRUE)[1]-1)
        }
      }
      colnames(pred_class) <- colnames(object$y)[spp_idx]
      if(is.null(out)){ out <- eta; out[] <- NA}
      out[, fam_pred == "ordinal"] <- pred_class[, fam_pred == "ordinal"]
    }
  }
  try(rownames(out) <- 1:NROW(out), silent = TRUE)
  if(any(class(out) %in% "dgeMatrix")) out <- as.matrix(out)
  
  if((is.numeric(se.fit) || se.fit) && isFALSE(object$sd)){
    stop("Cannot calculate standard errors for the prediction without standard errors of the parameter estimates. Please refit the model with 'sd.erors = TRUE', or use 'se.gllvm'.")
  }else if(is.numeric(se.fit) || se.fit){
    R <- 1e3
    if(is.numeric(se.fit))R <- se.fit

    params <- simulate.params.gllvm(object, R, seed, level, n = n)
    
    predSims <- array(0.0, dim = c(R, dim(out)))

    for(r in 1:R){
      newobject <- perturb.gllvm(object, params, r, type)
      pred_r <- predict(newobject, newX = newX, newTR = newTR, newLV = newLV,
                        type = type, level = level, offset = offset,
                        se.fit = FALSE, ordinal.cat = ordinal.cat, spp = spp)
      if(!any(object$family == "ordinal") || type == "link" || !is.null(ordinal.cat))
        predSims[r,,] <- pred_r
      else
        predSims[r,,,] <- pred_r
    }
    
    if(!any(object$family == "ordinal") || !is.null(alpha) && type == "link" || is.null(alpha) || !is.null(ordinal.cat)){
    if(is.null(alpha)){
    out <- list(fit = out, ci.sim = predSims) # so that we can access this for predictSR and predictPairwise
    }else if(!is.null(alpha)){
    ci <- apply(predSims, 2:3, quantile, prob = c((1-alpha)/2, 1-(1-alpha)/2))
    out <- list(fit = out, lower = ci[1,,], upper = ci[2,,])
    }
    }else if(any(object$family == "ordinal") && type == "response" && !is.null(alpha)){
      ci <- array(dim=c(2,nrow(out), n, length(spp_idx)))
      for(j in seq_len(length(spp_idx))){
        if(fam_pred[j] != "ordinal"){
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
            ci[,,i,j] <- suppressWarnings(apply(predSims[,,i,j][dists <= threshold, , drop=FALSE],2,range, na.rm = TRUE))
          }
          if(any(!is.finite(ci)))ci[!is.finite(ci)] <- NA # guaranteed to happen on zeta.struc = "species"
        }
      }
      out <- list(fit = out, lower = ci[1,,,], upper = ci[2,,,])
    }  
  }
  
  return(out)
}

# function to simulate new parameters from the asymptotic covariance matrix
# previously in predict.gllvm, but lifted out
# to separately call in predictSR.gllvm
# for improved efficiency
simulate.params.gllvm <- function(object, R, seed = 42, level = 1, n = NULL){
  set.seed(seed)
  
  incl  <- object$Hess$incl
  incla <- object$Hess$incla
  
  if(object$method %in% c("VA","EVA")){
    num.lv    <- object$num.lv
    num.RR    <- object$num.RR
    num.lv.c  <- object$num.lv.c
    num.lv.cor <- object$num.lvcor
    incla <- rep(FALSE, length(object$TMBfn$par))
    if((num.lv.c + num.lv + num.lv.cor) > 0)
      incla[names(object$TMBfn$par) == "u"] <- TRUE
    if("Br" %in% names(object$TMBfn$par))
      incla[names(object$TMBfn$par) == "Br"] <- TRUE
    if((num.RR + num.lv.c) > 0 && !isFALSE(object$randomB))
      incla[names(object$TMBfn$par) == "b_lv"] <- TRUE
    if("r0r" %in% names(object$TMBfn$par))
      incla[names(object$TMBfn$par) == "r0r"] <- TRUE
  }
  if(is.null(incla)) incla <- 0
  
  ffs <- NULL
  Vf  <- NULL
  if(sum(incl) > 0){
    Vf <- vcov(object)
    if(min(diag(Vf)) / max(diag(Vf)) < .Machine$double.eps)
      warning("Seeing some odd things in the fixed effects covariance matrix (variances of effects are not on the same scale), uncertainties may be inaccurate.\n")
    ffs <- try(MASS::mvrnorm(R, object$TMBfn$par[incl], Vf), silent = TRUE)
    if(inherits(ffs, "try-error"))
      stop("Covariance matrix of fixed effects is not semi positive-definite.")
    colnames(ffs) <- names(object$TMBfn$par)[incl]
    
    if(any(colnames(ffs) %in% names(object$TMBfn$env$map))){
      map <- object$TMBfn$env$map[names(object$TMBfn$env$map) %in% colnames(ffs)]
      ffs.new <- NULL
      for(nm in unique(colnames(ffs))){
        if(!nm %in% names(map)){
          ffs.new <- cbind(ffs.new, ffs[, colnames(ffs) == nm, drop = FALSE])
        } else {
          ffs.new <- cbind(ffs.new, ffs[, colnames(ffs) == nm, drop = FALSE][, map[[nm]], drop = FALSE])
        }
      }
      ffs <- ffs.new
      rm(ffs.new)
    }
  }
  
  rfs <- NULL
  if(sum(incla) > 0 && level > 0){
    p <- ncol(object$y)
    if(is.null(n)) n <- nrow(object$y)
    if(object$method == "LA"){
      Vr <- sdrandom(object$TMBfn, Vf, object$Hess$incl, return.covb = TRUE)
    } else {
      Vr    <- CMSEPf(object, return.covb = TRUE)
      renms <- c("r0r","Br","u")
      if(!isFALSE(object$randomB)) renms <- c(renms, "b_lv")
      colnames(Vr) <- row.names(Vr) <- names(object$TMBfn$par[names(object$TMBfn$par) %in% renms])
      if(!is.null(object$Ar))
        Vr[row.names(Vr) == "r0r", colnames(Vr) == "r0r"] <-
        as.matrix(Vr[row.names(Vr) == "r0r", colnames(Vr) == "r0r"] + Matrix::bdiag(object$Ar))
      if(object$col.eff$col.eff == "random"){
        Ab <- switch(object$col.eff$Ab.struct,
                     diagonal      = ,
                     blockdiagonal = Matrix::bdiag(object$Ab),
                     diagonalCL2   = Matrix::bdiag(object$Ab)[
                       order(rep(seq_len(p), times = nrow(object$params$Br))),
                       order(rep(seq_len(p), times = nrow(object$params$Br)))],
                     unstructured  = object$Ab[[1]],
                     MNdiagonal    = ,
                     MNunstructured = kronecker(cov2cor(object$Ab[[2]]), object$Ab[[1]]),
                     object$Ab  # diagonalCL1 / CL1 / CL2
        )
        Vr[row.names(Vr) == "Br", colnames(Vr) == "Br"] <-
          Vr[row.names(Vr) == "Br", colnames(Vr) == "Br"] + Ab
      }
      if((object$num.lv + object$num.lv.c) > 0 && n == nrow(object$y)){
        A   <- lapply(seq(dim(object$A)[1]), function(i) object$A[i, , ])
        d   <- object$num.lv + object$num.lv.c
        idx <- rep(seq_len(n), each = d) + rep(c(0, rep(n, d - 1) * seq_len(d - 1)), times = n)
        Vr[row.names(Vr) == "u", colnames(Vr) == "u"] <-
          as.matrix(Vr[row.names(Vr) == "u", colnames(Vr) == "u"] +
                      Matrix::bdiag(A)[order(idx), order(idx)])
      }
      if(!isFALSE(object$randomB)){
        d     <- object$num.RR + object$num.lv.c
        K     <- nrow(object$params$LvXcoef)
        Ab.lv <- Matrix::bdiag(lapply(seq(dim(object$Ab.lv)[1]), function(x) object$Ab.lv[1, , ]))
        idx   <- seq_len(d * K)
        if(object$randomB == "LV")
          idx <- rep(seq_len(K), each = d) + rep(c(0, rep(K, d - 1) * seq_len(d - 1)), times = K)
        Vr[row.names(Vr) == "b_lv", colnames(Vr) == "b_lv"] <-
          as.matrix(Vr[row.names(Vr) == "b_lv", colnames(Vr) == "b_lv"] +
                      Ab.lv[order(idx), order(idx)])
      }
    }
    if(min(diag(Vr)) / max(diag(Vr)) < .Machine$double.eps)
      warning("Seeing some odd things in the random effects (asymptotic) covariance matrix (variances of effects are not on the same scale), uncertainties may be inaccurate.\n")
    rfs <- try(MASS::mvrnorm(R, object$TMBfn$env$last.par.best[incla], Vr), silent = TRUE)
    if(inherits(rfs, "try-error"))
      stop("Covariance matrix of random effects is not semi positive-definite.")
  }
  
  list(ffs = ffs, rfs = rfs, incl = incl, incla = incla)
}

# re-fills gllvm objects with simulations for se.fit
# previously in predict.gllvm, but lifted out
# to separately call in predictSR.gllvm
# for improved efficiency
perturb.gllvm <- function(object, params, r, type = "response", skeleton = NULL, template = NULL){
  ffs   <- params$ffs
  rfs   <- params$rfs
  incl  <- params$incl
  incla <- params$incla

  num.lv    <- object$num.lv
  num.RR    <- object$num.RR
  num.lv.c  <- object$num.lv.c
  quadratic <- object$quadratic
  p         <- ncol(object$y)

  if(is.null(skeleton)) skeleton <- object$TMBfn$env$parList()

  # blank all params so nothing leaks from the fitted values
  if(!is.null(template)){
    newobject <- template
  } else {
    newobject <- object
    newobject$params <- lapply(object$params, function(x){ if(is.numeric(x)) x[] <- NA; x })
    if(!is.null(newobject$lvs)) newobject$lvs[] <- NA
  }

  if(sum(incl) > 0){
    newpars <- relist_gllvm(ffs[r, ], skeleton)
    
    newobject$params$beta0 <- newpars$b[1, ]
    names(newobject$params$beta0) <- names(object$params$beta0)
    if(nrow(newpars$b) > 1){
      newobject$params$Xcoef <- t(newpars$b[-1, , drop = FALSE])
      rownames(newobject$params$Xcoef) <- rownames(object$params$Xcoef)
      colnames(newobject$params$Xcoef) <- colnames(object$params$Xcoef)
    }

    if((num.lv + num.lv.c) > 0){
      newobject$params$sigma.lv <- abs(newpars$sigmaLV)
      names(newobject$params$sigma.lv) <- names(object$params$sigma.lv)
    }
    
    if((num.lv + num.lv.c + num.RR) > 0){
      theta <- matrix(0, p, num.lv + num.lv.c + num.RR)
      if((num.lv.c + num.RR) > 1)  diag(theta[, seq_len(num.lv.c + num.RR)]) <- 1
      else if((num.lv.c + num.RR) == 1) theta[1, 1] <- 1
      if(num.lv > 1)  diag(theta[, (num.lv.c + num.RR + 1):(num.lv.c + num.RR + num.lv)]) <- 1
      else if(num.lv == 1) theta[1, num.lv.c + num.RR + 1] <- 1
      
      if(num.lv > 0 && (num.lv.c + num.RR) == 0){
        if(p > 1){
          theta[lower.tri(theta[, seq_len(num.lv), drop = FALSE], diag = FALSE)] <- newpars$lambda
          if(quadratic != FALSE)
            theta <- cbind(theta, matrix(-abs(newpars$lambda2), ncol = num.lv, nrow = p, byrow = TRUE))
        } else {
          theta <- if(quadratic == FALSE) as.matrix(1) else c(as.matrix(1), -abs(newpars$lambda2))
        }
      } else if(num.lv == 0 && (num.lv.c + num.RR) > 0){
        if(p > 1){
          theta[lower.tri(theta[, seq_len(num.lv.c + num.RR), drop = FALSE], diag = FALSE)] <- newpars$lambda
          if(quadratic != FALSE)
            theta <- cbind(theta, matrix(-abs(newpars$lambda2), ncol = num.lv.c + num.RR, nrow = p, byrow = TRUE))
        } else {
          theta <- if(quadratic == FALSE) as.matrix(1) else c(as.matrix(1), -abs(newpars$lambda2))
        }
      } else {
        if(p > 1){
          nc1 <- num.lv.c + num.RR
          theta[, seq_len(nc1)][lower.tri(theta[, seq_len(nc1), drop = FALSE], diag = FALSE)] <-
            newpars$lambda[seq_len(sum(lower.tri(theta[, seq_len(nc1), drop = FALSE], diag = FALSE)))]
          theta[, (nc1 + 1):ncol(theta)][lower.tri(theta[, (nc1 + 1):ncol(theta), drop = FALSE], diag = FALSE)] <-
            newpars$lambda[(sum(lower.tri(theta[, seq_len(nc1), drop = FALSE], diag = FALSE)) + 1):length(newpars$lambda)]
          if(quadratic != FALSE)
            theta <- cbind(theta, matrix(-abs(newpars$lambda2), ncol = num.lv + nc1, nrow = p, byrow = TRUE))
        } else {
          theta <- if(quadratic == FALSE) as.matrix(1) else c(as.matrix(1), -abs(newpars$lambda2))
        }
      }
      rownames(theta) <- rownames(object$params$theta)
      colnames(theta) <- colnames(object$params$theta)
      newobject$params$theta <- theta
    }

    if(any(object$family %in% c("ZIP","ZINB"))){
      lp0 <- newpars$lg_phi[object$disp.group]
      newobject$params$phi <- object$params$phi
      newobject$params$phi[object$family %in% c("ZIP","ZINB")] <- (exp(lp0) / (1 + exp(lp0)))[object$family %in% c("ZIP","ZINB")]
    }
    
    if(any(object$family %in% c("orderedBeta","ordinal")) && type == "response"){
      K <- if(any(object$family %in% "ordinal")) max(object$TMBfn$env$data$y, na.rm =TRUE) - min(object$TMBfn$env$data$y, na.rm =TRUE) else 2L
      if(object$zeta.struc == "common"){
        zetaO <- c(rep(TRUE,  2L * any(object$family %in% "orderedBeta")),
                   rep(FALSE, (K - 1L) * any(object$family %in% "ordinal")))
      }
      zetas <- newpars$zeta
      if(object$zeta.struc == "species"){
        zetanew <- matrix(NA, nrow = p, ncol = K); zetanew[, 1] <- 0; idx <- 0
        for(j in which(object$family %in% c("ordinal","orderedBeta"))){
          if(object$family[j] == "ordinal"){
            k <- max(object$y[, j], na.rm = TRUE) - 2
            if(k > 0) zetanew[j, 2:(k + 1)] <- zetas[idx + seq_len(k)]
            idx <- idx + k
            zetanew[j, ] <- cumsum(abs(zetanew[j, ]))
          } else {
            zetanew[j, ] <- c(zetas[idx + 1], exp(zetas[idx + 2])); idx <- idx + 2
          }
        }
      } else {
        zetanew <- NULL
        if(any(object$family == "orderedBeta")){
          zetanew <- c(zetas[1], exp(zetas[2])); names(zetanew) <- c("cutoff0","cutoff1")
        }
        if(any(object$family %in% "ordinal"))
          zetanew <- c(zetanew, 0, cumsum(abs(zetas[!zetaO])))
      }
      newobject$params$zeta <- zetanew
    }
    
    if(!is.null(newpars$B)){
      newobject$params$B <- newpars$B[!is.na(newpars$B)]
      names(newobject$params$B) <- names(object$params$B)
    }
    if(!is.null(newpars$r0f)){
      newobject$params$row.params.fixed <- c(newpars$r0f)
      names(newobject$params$row.params.fixed) <- names(object$params$row.params.fixed)
    }
    if(num.RR > 0 && isFALSE(object$randomB)){
      newobject$params$LvXcoef <- newpars$b_lv
      rownames(newobject$params$LvXcoef) <- rownames(object$params$LvXcoef)
      colnames(newobject$params$LvXcoef) <- colnames(object$params$LvXcoef)
    }
  }

  if(sum(incla) > 0){
    newrfs <- relist_gllvm(rfs[r, ], skeleton)
    if(!is.null(newrfs$Br)){
      newobject$params$Br <- newrfs$Br
      rownames(newobject$params$Br) <- rownames(object$params$Br)
      colnames(newobject$params$Br) <- colnames(object$params$Br)
    }
    if(!is.null(newrfs$u))    newobject$lvs        <- newrfs$u
    if(!is.null(newrfs$b_lv)){
      newobject$params$LvXcoef <- newrfs$b_lv
      rownames(newobject$params$LvXcoef) <- rownames(object$params$LvXcoef)
      colnames(newobject$params$LvXcoef) <- colnames(object$params$LvXcoef)
    }
    if(!is.null(newrfs$r0r)){
      newobject$params$row.params.random <- c(newrfs$r0r)
      names(newobject$params$row.params.random) <- names(object$params$row.params.random)
    }
  }
  
  newobject
}
