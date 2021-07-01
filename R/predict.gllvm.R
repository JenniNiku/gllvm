#' @title Predict Method for gllvm Fits
#' @description Obtains predictions from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param type the type of prediction required. The default (\code{"link"}) is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable. that is, the predictions for the binomial model are predicted probabilities. In case of ordinal data, \code{type = "response"} gives predicted probabilities for each level of ordinal variable.
#' @param newX A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.
#' @param newTR A new data frame of traits for each response taxon. If omitted, the original matrix of traits is used.
#' @param newLV A new matrix of latent variables.  If omitted, the original matrix of latent variables is used.
#' @param level specification for how to predict. Level one attempts to use the predicted site scores from variational approximations of laplace approximation. Level 0 sets the latent variable to zero instead. Defaults to 1.
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
#'@export
#'@export predict.gllvm

predict.gllvm <- function(object, newX = NULL, newTR = NULL, newLV = NULL, type ="link", level = 1, ...){
  r0 <- NULL
  # if(object$row.eff!=FALSE&&(object$num.lv>0|object$num.lv.c>0)&&ncol(newLV)==(object$num.lv+object$num.lv.c+1)){
  #   warning("First column of newLV taken to be row-intercepts. \n")
  #   r0 <- newLV[,1];newLV<-newLV[,-1,drop=F]
  # }
  newdata <- newX
  p <- ncol(object$y)
  n <- max(nrow(object$y),nrow(newdata), nrow(newLV))
  if(!is.null(newdata)) n <- nrow(newdata)
  # if(n!=nrow(object$y)&is.null(newLV)&(object$num.lv.c+object$num.lv)>0)stop("With new predictors with a different number of sites, new latent variables need to be provided through the newLV argument.")
  if(is.null(newdata) && !is.null(object$X) && !is.null(newLV) && (nrow(newLV) != nrow(object$y))) stop("Number of rows in newLV must equal to the number of rows in the response matrix, if environmental variables are included in the model and newX is not included.") 
  if(!is.null(object$X)){formula <- formula(terms(object))}else{formula<-NULL}
  
  if(object$row.eff != FALSE) {
    if(length(object$params$row.params) != nrow(object$y))
      object$params$row.params = c(object$TMBfn$env$data$dr0 %*% object$params$row.params)
  }
  
  b0 <- object$params$beta0
  eta <- matrix(b0, n, p, byrow = TRUE)
  if(!is.null(newTR)) if(nrow(newTR) != p) stop("Number of rows in newTR must match to the number of responses in the original data matrix.")

  if(is.null(colnames(object$y))){
    colnames(object$y) <- paste("y",1:p, sep = "")
  }
    
    
  if(!is.null(object$X) && is.null(object$TR)) {
    B <- object$params$Xcoef
    X.d <- Xnew <- object$X
    if(!is.null(newdata)) {
      if(is.null(object$call$formula)){
        n1 <- colnames(newdata)
        formula1 <- paste("~ ", n1[1], sep = "")
        if(length(n1) > 1){
          for(i1 in 2:length(n1)){
            formula1 <- paste(formula1, n1[i1], sep = "+")
          }}
        formula1 <- paste(formula1, "1", sep = "-")
        formula1 <- formula(formula1)
        Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      }
      formula <- formula(object$formula)
      xb <- as.matrix(model.matrix(formula, data = data.frame(Xnew)))
      X.d <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
      colnames(X.d) <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
    } else {
      X.d <- object$X.design
      }
    eta <- eta + X.d %*% t(B)
  }
  
  if(!is.null(object$X) && !is.null(object$TR)) {
    B <- object$params$B
    X.d <- object$X.design
    if(!is.null(newdata) && !is.null(newTR)) {
      Xnew <- newdata
      TRnew <- newTR
      y1 <- object$y[sample(1:nrow(object$y),nrow(Xnew), replace = TRUE), ]
      # y1 <- object$y[1:nrow(Xnew), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", varying = colnames(data.frame(y1)), v.names = "y", timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if(!is.null(newdata) && is.null(newTR)) {
      formula <- formula(object$formula)
      # if(is.null(object$call$formula)){
        n1 <- colnames(newdata)
        formula1 <- paste("~", n1[1], sep = "")
        if(length(n1) > 1){
          for(i1 in 2:length(n1)){
            formula1 <- paste(formula1, n1[i1], sep = "+")
          }}
        formula1 <- paste(formula1, "1", sep = "-")
        formula1 <- formula(formula1)
        Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      # }
      TRnew <- object$TR
      y1 <- object$y[sample(1:nrow(object$y),nrow(Xnew), replace = TRUE), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", varying = colnames(data.frame(y1)), v.names = "y", timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if(is.null(newdata) && !is.null(newTR)) {
      if(nrow(newTR) != p) stop("Number of rows in TRnew must equal to the number of response variables.")
      formula <- formula(object$formula)
      TRnew <- newTR
      # if(is.null(object$call$formula)){
        n1 <- colnames(newTR)
        formula1 <- paste("~", n1[1], sep = "")
        if(length(n1) > 1){
          for(i1 in 2:length(n1)){
            formula1 <- paste(formula1, n1[i1], sep = "+")
          }}
        formula1 <- paste(formula1, "1", sep = "-")
        formula1 <- formula(formula1)
        TRnew <- as.matrix(model.matrix(formula1, data = data.frame(newTR)))
      # }
      Xnew <- object$X
      y1 <- object$y[sample(1:nrow(object$y),nrow(Xnew), replace = TRUE), ]
      # y1 <- object$y[1:nrow(Xnew), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", varying = colnames(data.frame(y1)), v.names = "y", timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    X.d<- as.matrix(X.d[, colnames(X.d)!= "(Intercept)"])
    eta <- eta + matrix(X.d %*% B, n, p)
    if(!is.null(object$randomX)){
      if(!is.null(newdata)) {
        tr<-try(X.xr <- as.matrix(model.matrix(object$randomX, data = yXT)), silent = TRUE)
        if(inherits(tr, "try-error")) { X.xr <- as.matrix(yXT[, colnames(object$Xrandom)]); colnames(X.xr)<-colnames(object$Xrandom) }
        rnam <- colnames(X.xr)[!(colnames(X.xr) %in% c("(Intercept)"))]
        xr <- as.matrix(X.xr[1:n,!(colnames(X.xr) %in% c("(Intercept)"))])
        # xr <- as.matrix(xr[, rnam]); #as.matrix(X.new[, rnam])
        if(NCOL(xr) == 1) colnames(xr) <- rnam
      } else {
        xr <- object$Xrandom
      }
      eta <- eta + matrix(xr %*% object$params$Br, n, p)
    }
  }
  

  if(level==1){
  if(is.null(newLV) && !is.null(newdata)){ stop("Level 1 predictions cannot be calculated for new X values if new latent variable values are not given. Change to 'level = 0' predictions.")}
  
  if(object$num.lv > 0 |(object$num.lv.c+object$num.RR)>0) {
    
    if(!is.null(newLV)) {
      if(ncol(newLV) != (object$num.lv+object$num.lv.c)) stop("Number of latent variables in input doesn't equal to the number of latent variables in the model.")
      if(!is.null(newdata)){  
        if(nrow(newLV) != nrow(newdata)) stop("Number of rows in newLV must equal to the number of rows in newX, if newX is included, otherwise same as number of rows in the response matrix.") 
      }
      lvs <- t(t(newLV)*object$params$sigma.lv)
    }else{
      if(object$num.RR==0){
        lvs <- t(t(object$lvs)*object$params$sigma.lv)
      }else{
        if(object$num.lv.c>0){
          lvs<- cbind(t(t(object$lvs[,1:object$num.lv.c])*object$params$sigma.lv[1:object$num.lv.c]),matrix(0,ncol=object$num.RR,nrow=n),t(t(object$lvs[,-c(1:object$num.lv.c)])*object$params$sigma.lv[1:object$num.lv]))
        }else if(object$num.lv>0&object$num.lv.c==0){
          lvs<- cbind(matrix(0,ncol=object$num.RR,nrow=n),t(t(object$lvs)*object$params$sigma.lv))
        }else{
          lvs <- matrix(0,ncol=object$num.RR,nrow=n)
        }
      }
    }

    if((object$num.lv.c+object$num.RR)>0&!is.null(newdata)){lv.X <-  as.matrix(model.frame(object$lv.formula,as.data.frame(newdata)))}else{lv.X<-object$lv.X}
    theta <- object$params$theta[,1:(object$num.lv+(object$num.lv.c+object$num.RR))]
      eta <- eta + lvs %*% t(theta)
      if((object$num.lv.c+object$num.RR)>0){
        eta <- eta + lv.X%*%object$params$LvXcoef%*%t(theta[,1:(object$num.lv.c+object$num.RR)])
      }
    if(object$quadratic != FALSE){
      theta2 <- object$params$theta[,-c(1:(object$num.lv+(object$num.lv.c+object$num.RR))),drop=F]
      eta <- eta + lvs^2 %*% t(theta2)
      if((object$num.lv.c+object$num.RR)>0){
        theta2C <- abs(theta2[,1:(object$num.lv.c+object$num.RR),drop=F])
        for(i in 1:n){
          for(j in 1:p){
            eta[i,j]<- eta[i,j] - 2*lvs[i,1:(object$num.lv.c+object$num.RR),drop=F]%*%diag(theta2C[j,])%*%t(lv.X[i,,drop=F]%*%object$params$LvXcoef) - lv.X[i,,drop=F]%*%object$params$LvXcoef%*%diag(theta2C[j,])%*%t(lv.X[i,,drop=F]%*%object$params$LvXcoef)
          }
      }
    }
    }

  }
  }

 
  
  if(object$row.eff %in% c("random", "fixed", "TRUE") && nrow(eta)==length(object$params$row.params)&is.null(r0)) {
    r0 <- object$params$row.params
    eta <- eta + r0
  }else if(!is.null(r0)){
    eta <- eta+r0
  }
  
  if(object$family %in% c("poisson", "negative.binomial", "tweedie", "gamma", "exponential"))
    ilinkfun <- exp
  if(object$family == "binomial" || object$family == "beta")
    ilinkfun <- binomial(link = object$link)$linkinv
  if(object$family == "ordinal")
    ilinkfun <- pnorm
  if(object$family == "ZIP")
    ilinkfun <- function(eta) exp(eta)*(1-matrix(object$params$phi, n, p, byrow = TRUE))
  if(object$family == "gaussian")
    ilinkfun <- gaussian()$linkinv
  
    
  out <- NULL
  preds <- NULL
  if("link" %in% type)
    out <- eta
  if("response" %in% type)
    out <- ilinkfun(eta)
  if(is.null(newdata) && is.null(newTR) && is.null(newLV) && "logL" %in% type)
    out <- object$logL
    
  if (object$family == "ordinal" && type == "response") {
  if (object$zeta.struc == "species"){
    k.max <- apply(object$params$zeta,1,function(x)length(x[!is.na(x)])) + 1
    preds <- array(NA, dim = c(max(k.max), nrow(eta), p), dimnames = list(paste("level", 1:max(k.max), sep = ""), NULL, NULL))
    
    for (i in 1:n) {
      for (j in 1:p) {
        probK <- NULL
        probK[1] <- pnorm(object$params$zeta[j, 1] - eta[i, j], log.p = FALSE)
        probK[k.max[j]] <- 1 - pnorm(object$params$zeta[j, k.max[j] - 1] - eta[i, j])
        if(k.max[j] > 2) {
          j.levels <- 2:(k.max[j] - 1)
          for(k in j.levels) { probK[k] <- pnorm(object$params$zeta[j, k] - eta[i, j]) - pnorm(object$params$zeta[j, k - 1] - eta[i, j]) }
        }
        preds[, i, j] <- c(probK, rep(NA, max(k.max) -k.max[j] ))
      }}
    out <- preds
  }else{
    k.max <- length(object$params$zeta) + 1
    preds <- array(NA, dim = c(k.max, nrow(eta), p), dimnames = list(paste("level", 1:max(k.max), sep = ""), NULL, NULL))
    
    for (i in 1:n) {
      for (j in 1:p) {
        probK <- NULL
        probK[1] <- pnorm(object$params$zeta[1] - eta[i, j], log.p = FALSE)
        probK[k.max] <- 1 - pnorm(object$params$zeta[k.max - 1] - eta[i, j])
          levels <- 2:(k.max - 1)
          for(k in levels) { probK[k] <- pnorm(object$params$zeta[k] - eta[i, j]) - pnorm(object$params$zeta[k - 1] - eta[i, j]) }
        preds[, i, j] <- c(probK)
      }}
    out <- preds
  }
    dimnames(preds)[[3]] <- colnames(object$y)
  }
   try(rownames(out)<-1:NROW(out), silent = TRUE)
  
  return(out)
}


#' @export predict
predict <- function(object, ...)
{
  UseMethod(generic = "predict")
}
