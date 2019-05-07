#' @title Predict Method for gllvm Fits
#' @description Obtains predictions from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param type the type of prediction required. The default (\code{"link"}) is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable. that is, the predictions for the binomial model are predicted probabilities. In case of ordinal data, \code{type = "response"} gives predicted probabilities for each level of ordinal variable.
#' @param newX A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.
#' @param newTR A new data frame of traits for each response taxon. If omitted, the original matrix of traits is used.
#' @param newLV A new matrix of latent variables.  If omitted, the original matrix of latent variables is used.
#' @param ... not used.
#'
#' @details
#' If \code{newX}, \code{newTR} and \code{newLV} are omitted the predictions are based on the data used for fitting the model. Notice that \code{newTR} need to match with the number of species in the original data.
#' Instead, new sites can be specified in \code{newX}. If predictors \code{newX} (and \code{newTR}) are given, and \code{newLV} is not, latent variables are not used in the predictions.
#' 
#' @return A matrix containing requested predictor types.
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- scale(antTraits$env[, 1:3])
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, family = poisson())
#'# fitted values
#'predfit <- predict(fit, type = "response")
#'# linear predictors
#'predlin <- predict(fit)
#'# Predict new sites:
#'# Generate matrix of environmental variables for 10 new sites
#'xnew <- cbind(rnorm(10), rnorm(10), rnorm(10))
#'colnames(xnew) <- colnames(X)
#'predfit <- predict(fit, newX = xnew, type = "response")
#'
#' \donttest{
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
#'predfit <- predict(fitt, newX = xnew, newTR = trnew, type = "response")
#'}
#'@aliases predict predict.gllvm
#'@method predict gllvm
#'@export
#'@export predict.gllvm
predict.gllvm <- function(object, newX = NULL, newTR = NULL, newLV = NULL, type ="link", ...){
  newdata <- newX
  p <- ncol(object$y)
  n <- nrow(object$y)
  if(!is.null(newdata)) n <- nrow(newdata)
  formula <- formula(terms(object))
  
  b0 <- object$params$beta0
  eta <- matrix(b0, n, p, byrow = TRUE)
  if(!is.null(newTR)) if(nrow(newTR) != p) stop("Number of rows in newTR must match to the number of responses in the original data matrix.")
  
  if(!is.null(object$X) && is.null(object$TR)) {
    B <- object$params$Xcoef
    X.d <- Xnew <- object$X
    if(!is.null(newdata)) {
      formula <- formula(object$formula)
      n1 <- colnames(newdata)
      formula1 <- paste("~ ", n1[1], sep = "")
      if(length(n1) > 1){
        for(i1 in 2:length(n1)){
          formula1 <- paste(formula1, n1[i1], sep = "+")
        }}
      formula1 <- formula(formula1)
      Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      xb <- as.matrix(model.matrix(formula, data = data.frame(Xnew)))
      X.d <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
      colnames(X.d) <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      }
    eta <- eta + X.d %*% t(B)
  }
  
  if(!is.null(object$X) && !is.null(object$TR)) {
    B <- object$params$B
    X.d <- object$X.design
    if(!is.null(newdata) && !is.null(newTR)) {
      Xnew <- newdata
      TRnew <- newTR
      y1 <- object$y[1:nrow(Xnew), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", varying = colnames(y1), v.names = "y")
      TR2 <- data.frame(time = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "time")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if(!is.null(newdata) && is.null(newTR)) {
      formula <- formula(object$formula)
      n1 <- colnames(newdata)
      formula1 <- paste("~", n1[1], sep = "")
      if(length(n1) > 1){
        for(i1 in 2:length(n1)){
          formula1 <- paste(formula1, n1[i1], sep = "+")
        }}
      formula1 <- formula(formula1)
      Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      
      TRnew <- object$TR
      y1 <- object$y[1:nrow(Xnew), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", varying = colnames(y1), v.names = "y")
      TR2 <- data.frame(time = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "time")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if(is.null(newdata) && !is.null(newTR)) {
      if(nrow(newTR) != p) stop("Number of rows in TRnew must equal to the number of response variables.")
      formula <- formula(object$formula)
      n1 <- colnames(newTR)
      formula1 <- paste("~", n1[1], sep = "")
      if(length(n1) > 1){
        for(i1 in 2:length(n1)){
          formula1 <- paste(formula1, n1[i1], sep = "+")
        }}
      formula1 <- formula(formula1)
      TRnew <- as.matrix(model.matrix(formula1, data = data.frame(newTR)))
      Xnew <- object$X
      y1 <- object$y[1:nrow(Xnew), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", varying = colnames(y1), v.names = "y")
      TR2 <- data.frame(time = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "time")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    X.d<- as.matrix(X.d[, colnames(X.d)!= "(Intercept)"])
    eta <- eta + matrix(X.d %*% B, n, p)
  }
  
  
  if(object$num.lv > 0) {
    theta <- object$params$theta
    if(is.null(newLV) && is.null(newdata) && is.null(newTR))
      eta <- eta + object$lvs %*% t(theta)
    if(!is.null(newLV)) {
      if(ncol(newLV) != object$num.lv) stop("Number of latent variables in input doesn't equal to the number of latent variables in the model.")
      if(is.null(newdata)==FALSE) #DW addition, 6/5/19: so no error here for intercept models
      {  if(nrow(newLV) != nrow(Xnew)) stop("Number of rows in newLV must equal to the number of rows in newX, if newX is included, otherwise same as number of rows in the response matrix.") }
      lvs <- newLV
      eta <- eta + lvs %*% t(theta)
    }
  }
  
  if(object$row.eff %in% c("random", "fixed", "TRUE")) {
    r0 <- object$params$row.params
    eta <- eta + r0
  }
  
  if(object$family %in% c("poisson", "negative.binomial", "tweedie"))
    ilinkfun <- exp
  if(object$family == "binomial")
    ilinkfun <- binomial(link = object$link)$linkinv
  if(object$family == "ordinal")
    ilinkfun <- pnorm
  if(object$family == "ZIP")
    ilinkfun <- pnorm
  
  out <- NULL
  preds <- NULL
  if("link" %in% type)
    out <- eta
  if("response" %in% type)
    out <- ilinkfun(eta)
  if(is.null(newdata) && is.null(newTR) && is.null(newLV) && "logL" %in% type)
    out <- object$logL
    
  if (object$family == "ordinal" && type == "response") {
    k.max <- ncol(object$params$zeta) + 1
    preds <- array(NA, dim = c(k.max, nrow(eta), p), dimnames = list(paste("level", 1:k.max, sep = ""), NULL, NULL))
    
    for (i in 1:n) {
      for (j in 1:p) {
        probK <- NULL
        probK[1] <- pnorm(object$params$zeta[j, 1] - eta[i, j], log.p = FALSE)
        probK[k.max] <- 1 - pnorm(object$params$zeta[j, k.max - 1] - eta[i, j])
        if(k.max > 2) {
          j.levels <- 2:(k.max - 1)
          for(k in j.levels) { probK[k] <- pnorm(object$params$zeta[j, k] - eta[i, j]) - pnorm(object$params$zeta[j, k - 1] - eta[i, j]) }
        }
        preds[, i, j] <- probK
      }}
    out <- preds
  }
  return(out)
}


#' @export predict
predict <- function(object, ...)
{
  UseMethod(generic = "predict")
}
