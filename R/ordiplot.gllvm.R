#' @title Plot latent variables from gllvm model
#' @description Plots latent variables and their corresponding coefficients (biplot).
#'
#' @param object   an object of class 'gllvm'.
#' @param biplot   \code{TRUE} if both latent variables and their coefficients are plotted, \code{FALSE} if only LVs.
#' @param ind.spp  the number of response variables (usually, species) to include on the biplot. The default is none, or all if \code{biplot = TRUE}.
#' @param alpha    a numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot.
#' @param main  main title.
#' @param which.lvs indices of two latent variables to be plotted if number of the latent variables is more than 2. A vector with length of two. Defaults to \code{c(1,2)}.
#' @param jitter   if \code{TRUE}, jittering is applied on points.
#' @param ...	additional graphical arguments.
#'
#' @details
#' Function constructs a scatter plot of two latent variables, i.e. an ordination plot. If only one latent
#' variable is in the fitted model, latent variables are plotted against their corresponding row indices.
#' The latent variables are labeled using the row index of the response matrix y.
#'
#' Coefficients related to latent variables are plotted in the same figure with the latent
#' variables if \code{biplot = TRUE}. They are labeled using the column names of y. The number
#' of latent variable coefficients to be plotted can be controlled by ind.spp. An argument alpha
#' is used to control the relative scaling of the latent variables and their coefficients.
#' If \code{alpha = 0.5}, the latent variables and their coefficients are on the same scale.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# Ordination plot:
#'ordiplot.gllvm(fit)
#'# Biplot with 10 species
#'ordiplot.gllvm(fit, biplot = TRUE, ind.spp = 10)
#'}
#'@export


ordiplot.gllvm <- function(object, biplot = FALSE, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = c(1, 2), jitter = FALSE, ...) {
  if (any(class(object) != "gllvm"))
    stop("Class of the object isn't 'gllvm'.")

  n <- NROW(object$y)
  p <- NCOL(object$y)
  if (!is.null(ind.spp)) {
    ind.spp <- min(c(p, ind.spp))
  } else {
    ind.spp <- p
  }
  if (object$num.lv == 0)
    stop("No latent variables to plot.")

  if (is.null(rownames(object$params$theta)))
    rownames(object$params$theta) = paste("V", 1:p)

  par(mfrow = c(1, 1), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1.2)

  if (object$num.lv == 1) {
    plot(1:n, object$lvs, ylab = "LV1", xlab = "Row index")
  }

  if (object$num.lv > 1) {
    testcov <- object$lvs %*% t(object$params$theta)
    do.svd <- svd(testcov, object$num.lv, object$num.lv)
    choose.lvs <- scale(
      do.svd$u * matrix(do.svd$d[1:object$num.lv] ^ alpha, nrow = n, ncol = object$num.lv, byrow = TRUE),
      center = TRUE, scale = F)
    choose.lv.coefs <- scale(
      do.svd$v * matrix( do.svd$d[1:object$num.lv] ^ (1 - alpha),nrow = p,ncol = object$num.lv,byrow = TRUE),
      center = TRUE,scale = F)

    if (!biplot) {
      plot(object$lvs[, which.lvs],xlab = paste("Latent variable ", which.lvs[1]),
           ylab = paste("Latent variable ", which.lvs[2]),main = main,cex = 1.2 ,type = "n", ...)
      if (!jitter)
        text(object$lvs[, which.lvs], label = 1:n, cex = 1.2)
      if (jitter)
        text(jitter(object$lvs[, which.lvs][, 1]),jitter(object$lvs[, which.lvs][, 2]),label = 1:n,cex = 1.2)
    }

    if (biplot) {
      largest.lnorms <- order(rowSums(choose.lv.coefs[, which.lvs] ^ 2), decreasing = TRUE)[1:ind.spp]
      plot(rbind(choose.lvs[, which.lvs], choose.lv.coefs[, which.lvs]),xlab = paste("Latent variable ", which.lvs[1]),ylab = paste("Latent variable ", which.lvs[2]),main = main,type = "n",
           xlim = 1.1 * range(rbind(choose.lvs[, which.lvs], choose.lv.coefs[, which.lvs])[, 1]),
           ylim = 1.1 * range(rbind(choose.lvs[, which.lvs], choose.lv.coefs[, which.lvs])[, 2]))
      if (!jitter)
        text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2)
      if (jitter)
        text(jitter(choose.lvs[, which.lvs][, 1]),jitter(choose.lvs[, which.lvs][, 2]),label = 1:n,cex = 1.2)
      text(choose.lv.coefs[, which.lvs][largest.lnorms,],label = rownames(object$params$theta[, which.lvs][largest.lnorms,]),col = 4,cex = 0.9)
    }

  }
}
