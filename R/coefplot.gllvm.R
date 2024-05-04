#' @title Plot covariate coefficients and confidence intervals
#' @description Plots covariate coefficients and their confidence intervals.
#'
#' @param object an object of class 'gllvm'.
#' @param y.label logical, if \code{TRUE} (default) colnames of y with respect to coefficients are added to plot.
#' @param which.Xcoef vector indicating which covariate coefficients will be plotted. Can be vector of covariate names or numbers. Default is \code{NULL} when all covariate coefficients are plotted.
#' @param order logical, whether or not coefficients are ordered, defaults to \code{TRUE}.
#' @param cex.ylab the magnification to be used for axis annotation relative to the current setting of cex.
#' @param cex.xlab the magnification to be used for axis annotation.
#' @param mfrow same as \code{mfrow} in \code{par}. If \code{NULL} (default) it is determined automatically.
#' @param mar vector of length 4, which defines the margin sizes: \code{c(bottom, left, top, right)}. Defaults to \code{c(4,5,2,1)}.
#' @param xlim.list list of vectors with length of two to define the intervals for an x axis in each covariate plot. Defaults to NULL when the interval is defined by the range of point estimates and confidence intervals
#' @param ...	additional graphical arguments.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui, Sara Taskinen, Bert van der Veen
#'
#' @examples
#'# Extract subset of the microbial data to be used as an example
#'data(microbialdata)
#'X <- microbialdata$Xenv
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, X, formula = ~ pH + Phosp, family = poisson())
#'coefplot(fit)
#'\dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- as.matrix(antTraits$env)
#'# Fit model with environmental covariates
#'fit <- gllvm(y, X, formula = ~ Bare.ground + Shrub.cover,
#'             family = poisson())
#'coefplot.gllvm(fit)
#'
#'# Fit model with all environmental covariates
#'fitx <- gllvm(y, X, family = "negative.binomial")
#'coefplot(fitx, mfrow = c(3,2))
#'coefplot(fitx, which.Xcoef = 1:2)
#'
#'# Fit gllvm model with environmental and trait covariates
#'TR <- antTraits$traits
#'fitT <- gllvm(y = y, X = X, TR = TR, family = "negative.binomial")
#'coefplot(fitT)
#'
#'# Fit  gllvm model with environmental covariances and reduced rank
#'fitRR <- gllvm(y = y, X = X, num.RR = 2, family = "negative.binomial")
#'coefplot(fitRR)
#'}
#'@aliases coefplot coefplot.gllvm
#'@export
#'@export coefplot.gllvm
coefplot.gllvm <- function(object, y.label = TRUE, which.Xcoef = NULL, order = TRUE, cex.ylab = 0.5, cex.xlab = 1.3, mfrow = NULL, mar = c(4,6,2,1), xlim.list = NULL, ...)
{

  if (!any(class(object) == "gllvm"))
    stop("Class of the object isn't 'gllvm'.")
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  if (is.null(object$X) & is.null(object$lv.X.design))
    stop("No X covariates in the model.")
  if(is.null(object$X) && !isFALSE(object$randomB))
    stop("No predictor effects to plot in the model.")

  #Calculate standard errors of species-effects for reduced rank terms
  if(!is.null(object$lv.X.design) && isFALSE(object$randomB)){
    beta <- object$params$theta[,1:(object$num.lv.c+object$num.RR), drop=FALSE]%*%t(object$params$LvXcoef)
    betaSE <- RRse(object)
    object$params$Xcoef<-cbind(object$params$Xcoef,beta)
    object$sd$Xcoef<-cbind(object$sd$Xcoef,betaSE)
  }

  if (is.null(object$TR)) {
    if (is.null(which.Xcoef))
      which.Xcoef <- c(1:NCOL(object$params$Xcoef))
      Xcoef <- (object$params$Xcoef[, which.Xcoef, drop=FALSE])
      sdXcoef <- as.matrix(object$sd$Xcoef[, which.Xcoef, drop=FALSE])
      
    # cnames <- colnames(object$params$Xcoef)[which.Xcoef]
    # Xcoef <- as.matrix(object$params$Xcoef[, which.Xcoef,drop=FALSE])
    cnames <- colnames(Xcoef)

    k <- length(cnames)
    labely <- rownames(Xcoef)
    m <- length(labely)
    Xc <- Xcoef
    if (is.null(mfrow) && k > 1)
      mfrow <- c(1, k)

    if (!is.null(mfrow))
      par(mfrow = mfrow, mar = mar)
    if (is.null(mfrow))
      par(mar = mar)
    for (i in 1:k) {
      Xc <- Xcoef[, i]
      if(nrow(Xcoef)<2) names(Xc) <- rownames(Xcoef)
      lower <- Xc - 1.96 * sdXcoef[, i]
      upper <- Xc + 1.96 * sdXcoef[, i]
      if(order) Xc <- sort(Xc)
      lower <- lower[names(Xc)]
      upper <- upper[names(Xc)]

      col.seq <- rep("black", m)
      col.seq[lower < 0 & upper > 0] <- "grey"

      At.y <- seq(1, m)
      if (!is.null(xlim.list[[i]])) {
        plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = xlim.list[[i]], pch = "x", cex.lab = cex.xlab, ... )
      } else {
        plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = c(min(lower), max(upper)), pch = "x", cex.lab = cex.xlab, ... )
      }

      segments( x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq )
      abline(v = 0, lty = 1)
      if (y.label)
        axis( 2, at = At.y, labels = names(Xc), las = 1, cex.axis = cex.ylab)
    }
  } else{
    Xcoef <- object$params$B
    xnames <- names(object$params$B)
    sdXcoef <- object$sd$B
    if(object$family=="betaH"){
      if(is.null(which.Xcoef)){
        Xcoef <- c(Xcoef,object$params$betaH)
        xnames <- c(xnames,paste("beta",names(object$params$betaH),sep = ":"))
        names(Xcoef)<-xnames
        sdXcoef <- c(sdXcoef,object$sd$betaH)
      } else {
        Xcoef <- c(unlist(object$params[which.Xcoef]))
        xnames <- names(Xcoef)
        sdXcoef <- c(unlist(object$sd[which.Xcoef]))
      }
    }
    m <- length(xnames)
    k <- 1
    if (is.null(mfrow))
      mfrow <- c(1, k)

    par(mfrow = mfrow, mar = mar)
    Xc <- Xcoef
    
    lower <- Xc - 1.96 * sdXcoef
    upper <- Xc + 1.96 * sdXcoef
    if(order) Xc <- sort(Xc)
    lower <- lower[names(Xc)]
    upper <- upper[names(Xc)]

    col.seq <- rep("black", m)
    col.seq[lower < 0 & upper > 0] <- "grey"

    At.y <- seq(1, m)
    if (!is.null(xlim.list)) {
      if (is.list(xlim.list))
        xlim.list <- xlim.list[[1]]
      plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = "Coefficients", xlim = xlim.list, pch = "x", cex.lab = cex.xlab, ... )
    } else {
      plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = "Coefficients", xlim = c(min(lower), max(upper)), pch = "x", cex.lab = cex.xlab, ...)
    }
    segments( x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq )
    abline(v = 0, lty = 1)
    if (y.label)
      axis( 2, at = At.y, labels = names(Xc), las = 1, cex.axis = cex.ylab)

  }
}

#'@export
coefplot <- function(object, ...)
{
  UseMethod(generic="coefplot")
}

