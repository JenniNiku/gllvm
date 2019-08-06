#' @title Plot latent variables from gllvm model
#' @description Plots latent variables and their corresponding coefficients (biplot).
#'
#' @param object   an object of class 'gllvm'.
#' @param biplot   \code{TRUE} if both latent variables and their coefficients are plotted, \code{FALSE} if only latent variables.
#' @param ind.spp  the number of response variables (usually, species) to include on the biplot. The default is none, or all if \code{biplot = TRUE}.
#' @param alpha    a numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot.
#' @param main  main title.
#' @param which.lvs indices of two latent variables to be plotted if number of the latent variables is more than 2. A vector with length of two. Defaults to \code{c(1,2)}.
#' @param jitter   if \code{TRUE}, jittering is applied on points.
#' @param s.colors colors for sites
#' @param symbols logical, if \code{TRUE} sites are plotted using symbols, if \code{FALSE} (default) site numbers are used
#' @param cex.spp size of species labels in biplot
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
#' #'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'fit <- gllvm(y, family = poisson())
#'# Ordination plot:
#'ordiplot(fit)
#'# Biplot with 10 species
#'ordiplot(fit, biplot = TRUE, ind.spp = 10)
#'
#'@aliases ordiplot ordiplot.gllvm
#'@export
#'@export ordiplot.gllvm
ordiplot.gllvm <- function(object, biplot = FALSE, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = c(1, 2),
                           jitter = FALSE, s.colors = 1, symbols = FALSE, cex.spp = 0.7, ...) {
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

    if (object$num.lv == 1) {
      plot(1:n, object$lvs, ylab = "LV1", xlab = "Row index")
    }

    if (object$num.lv > 1) {
      # testcov <- object$lvs %*% t(object$params$theta)
      # do.svd <- svd(testcov, object$num.lv, object$num.lv)
      # choose.lvs <- do.svd$u * matrix( do.svd$d[1:object$num.lv] ^ alpha,
      #     nrow = n, ncol = object$num.lv, byrow = TRUE )
      # choose.lv.coefs <- do.svd$v * matrix(do.svd$d[1:object$num.lv] ^ (1 - alpha),
      #     nrow = p, ncol = object$num.lv, byrow = TRUE )
      sdi<-(diag(cov(object$lvs)))
      sdt<-(diag(cov(object$params$theta)))
      indexscale <- diag(x = (sdi)^alpha, nrow = length(sdi))
      gammascale <- diag(x = (sdt)^(1-alpha), nrow = length(sdi))
      choose.lvs <- scale(object$lvs, center = FALSE)%*%indexscale;
      choose.lv.coefs <- scale(object$params$theta, center = FALSE)%*%gammascale;

      
      if (!biplot) {
        plot(object$lvs[, which.lvs],
          xlab = paste("Latent variable ", which.lvs[1]),
          ylab = paste("Latent variable ", which.lvs[2]),
          main = main , type = "n", ... )
        if (!jitter)
          if (symbols) {
            points(object$lvs[, which.lvs], col = s.colors, ...)
          } else {
            text(object$lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
          }
        if (jitter)
          if (symbols) {
            points(jitter(object$lvs[, which.lvs][, 1]), jitter(object$lvs[, which.lvs][, 2]), col =
                     s.colors, ...)
          } else {
            text(
              jitter(object$lvs[, which.lvs][, 1]),
              jitter(object$lvs[, which.lvs][, 2]),
              label = 1:n, cex = 1.2, col = s.colors )
          }
      }

      if (biplot) {
        largest.lnorms <- order(apply(object$params$theta ^ 2, 1, sum), decreasing = TRUE)[1:ind.spp]

        plot(
          rbind(choose.lvs[, which.lvs], choose.lv.coefs[, which.lvs]),
          xlab = paste("Latent variable ", which.lvs[1]),
          ylab = paste("Latent variable ", which.lvs[2]),
          main = main, type = "n", ... )

        if (!jitter){
          if (symbols) {
            points(choose.lvs[, which.lvs], col = s.colors, ...)
          } else {
            text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
          }
          text(
            matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)),
            label = rownames(object$params$theta)[largest.lnorms],
            col = 4, cex = cex.spp )
          }
        if (jitter){
          if (symbols) {
            points(jitter(choose.lvs[, which.lvs[1]]), jitter(choose.lvs[, which.lvs[2]]), col =
                     s.colors, ...)
          } else {
            text(
              jitter(choose.lvs[, which.lvs[1]]),
              jitter(choose.lvs[, which.lvs[2]]),
              label = 1:n, cex = 1.2, col = s.colors )
          }
        text(
          jitter(matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)), amount = 0.2),
          label = rownames(object$params$theta)[largest.lnorms],
          col = 4, cex = cex.spp )
        }
      }

    }
  }


#'@export ordiplot
ordiplot <- function(object, ...)
{
  UseMethod(generic = "ordiplot")
}

