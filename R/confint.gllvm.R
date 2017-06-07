#' @title Confidence intervals for model parameters
#' @description Computes confidence intervals for  parameters in a fitted GLLVM model.
#'
#' @param object An object of class 'gllvm'.
#' @param level The confidence level. Scalar between 0 and 1.
#' @param ...	Not used.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- as.matrix(antTraits$env[,4:5])
#'# Fit GLLVM model
#'fit <- gllvm(y = y, X = X, family = "negative.binomial")
#'# 95 % confidence intervals
#'confint(fit, level = 0.95)
#'}
#'@export



confint.gllvm <- function(object, level = 0.95, ...) {
  if(is.logical(object$sd)) stop("Standard errors for parameters haven't been calculated, so confidence intervals can not be calculated.");
  n <- NROW(object$y)
  p <- NCOL(object$y)
  nX <- 0; if(!is.null(object$X)) nX <- dim(object$X)[2]
  nTR <- 0; if(!is.null(object$TR)) nTR <- dim(object$TR)[2]
  num.lv <- object$num.lv
  if (object$family=="negative.binomial") {
    object$params$phi <- NULL
    object$sd$phi <- NULL
  }
  alfa <- (1 - level) / 2
  cilow <- unlist(object$params) + qnorm(alfa) * unlist(object$sd)
  ciup <- unlist(object$params) + qnorm(1 - alfa) * unlist(object$sd)
  M <- cbind(cilow, ciup)

  colnames(M) <- c(paste(alfa * 100, "%"), paste((1 - alfa) * 100, "%"))
  rnames <- names(unlist(object$params))
  cal <- 0
  if (num.lv>0) {
    nr <- rep(1:num.lv, each=p)
    nc <- rep(1:p, num.lv)
    rnames[1:(num.lv * p)] <- paste(paste("theta.LV", nr, sep = ""),nc, sep = ".")
    cal <- cal+num.lv * p
  }
  rnames[(cal+1):(cal+p)] <- paste("Intercept",names(object$params$beta0), sep = ".")
  cal <- cal+p
  if (!is.null(object$TR)) {
    cal <- cal + nX
    if(object$get.fourth){
    nr <- rep(rownames(object$params$fourth), nX)
    nc <- rep(colnames(object$params$fourth), each = nTR)
    rnames[(cal + 1):(cal + nX * nTR)] <- paste(nr, nc, sep = "*")
    cal <- cal + nX * nTR}
    if(object$get.trait) cal <- cal + nTR
  }

  if (is.null(object$TR) && !is.null(object$X)) {
    cnx <- rep(colnames(object$X), each = p)
    rnc <- rep(rownames(object$params$Xcoef), nX)
    newnam <- paste(cnx, rnc, sep = ":")
    rnames[(cal + 1):(cal + nX * p)] <- paste("Xcoef", newnam, sep = ".")
    cal <- cal + nX * p
  }
  if (object$row.eff) {
    rnames[(cal+1):(cal+n)] <- paste("Row.Intercept",1:n, sep = ".")
    cal <- cal + n
  }
  if(object$family=="negative.binomial"){
    s <- dim(M)[1]
    rnames[(cal+1):s] <- paste("inv.phi", names(object$params$beta0), sep = ".")
  }
  if(object$family=="tweedie"){
    s <- dim(M)[1]
    rnames[(cal+1):s] <- paste("Dispersion phi", names(object$params$inv.phi), sep = ".")
  }
  if(object$family=="ZIP"){
    s <- dim(M)[1]
    rnames[(cal+1):s] <- paste("p", names(object$params$p), sep = ".")
  }
  rownames(M) <- rnames
  return(M)
}
