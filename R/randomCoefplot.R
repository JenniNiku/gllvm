#' @title Plot random slope coefficients
#' @description Plots random slopes and their prediction intervals.
#'
#' @param object an object of class 'gllvm'.
#' @param y.label logical, if \code{TRUE} (default) colnames of y with respect to coefficients are added to plot.
#' @param cex.ylab the magnification to be used for axis annotation relative to the current setting of cex.
#' @param mfrow same as \code{mfrow} in \code{par}. If \code{NULL} (default) it is determined automatically.
#' @param mar vector of length 4, which defines the margin sizes: \code{c(bottom, left, top, right)}. Defaults to \code{c(4,5,2,1)}.
#' @param xlim.list list of vectors with length of two to define the intervals for x axis in each covariate plot. Defaults to NULL when the interval is defined by the range of point estimates and confidence intervals
#' @param order logical, if \code{TRUE} (default), coefficients are sorted according to the point estimates
#' @param ...	additional graphical arguments.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui, Sara Taskinen
#'
#' @examples
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- as.matrix(antTraits$env)
#'TR <- antTraits$traits
#'# Fit model with random slopes
#'fitF <- gllvm(y = y, X = X, TR = TR,
#'  formula = ~ Bare.ground + Bare.ground : Webers.length,
#'  family = poisson(), randomX = ~ Bare.ground)
#'randomCoefplot(fitF)
#'
#'@aliases randomCoefplot randomCoefplot.gllvm
#'@export
#'@export randomCoefplot.gllvm
randomCoefplot.gllvm <- function(object, y.label = TRUE, cex.ylab = 0.5, mfrow = NULL, mar = c(4,6,2,1), xlim.list = NULL, order = FALSE, ...)
{
  
  if (any(class(object) != "gllvm"))
    stop("Class of the object isn't 'gllvm'.")
  
  if (is.null(object$Xrandom) || is.null(object$randomX))
    stop("No random covariates in the model.")
  
    which.Xcoef <- c(1:NROW(object$params$Br))
    Xcoef <- as.matrix(t(object$params$Br)[,which.Xcoef])
    cnames <- colnames(object$Xr)[which.Xcoef]
    k <- length(cnames)
    if(is.null(colnames(object$y))) 
      colnames(object$y) <- paste("Y",1:NCOL(object$y), sep = "")
    labely <- colnames(object$y)
    m <- length(labely)
    Xc <- Xcoef
    
    sdXcoef <- t(sqrt(apply(object$Ab,1,diag)))
    
    if (is.null(mfrow) && k > 1)
      mfrow <- c(1, k)
    if (!is.null(mfrow))
      par(mfrow = mfrow, mar = mar)
    if (is.null(mfrow))
      par(mar = mar)
    for (i in 1:k) {
      Xc <- Xcoef[, i]
      lower <- Xc - 1.96 * sdXcoef[, i]
      upper <- Xc + 1.96 * sdXcoef[, i]
      if(order){
        Xc <- sort(Xc)
        lower <- lower[names(Xc)]
        upper <- upper[names(Xc)]
      }
      col.seq <- rep("black", m)
      col.seq[lower < 0 & upper > 0] <- "grey"
      
      At.y <- seq(1, m)
      if (!is.null(xlim.list[[i]])) {
        plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = xlim.list[[i]], pch = "x", cex.lab = 1.3, ... )
      } else {
        plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = c(min(lower), max(upper)), pch = "x", cex.lab = 1.3, ... )
      }
      
      segments( x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq )
      abline(v = 0, lty = 1)
      if (y.label)
        axis( 2, at = At.y, labels = labely, las = 1, cex.axis = cex.ylab)
    }
}

#'@export
randomCoefplot <- function(object, ...)
{
  UseMethod(generic="randomCoefplot")
}

