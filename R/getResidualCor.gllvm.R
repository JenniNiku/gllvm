#' @title Extract residual correlations from gllvm object
#' @description  Calculates the residual correlation matrix for gllvm model.
#'
#' @param object   an object of class 'gllvm'.
#' @param adjust  The type of adjustment used for  negative binomial and binomial distribution when computing residual correlation matrix. Options are 0 (no adjustment), 1 (the default adjustment) and 2 (alternative adjustment for NB distribution). See details.
#' @param site.index A site index used used in the calculation of a GLLVM with quadratic response model, for which the residual correlations are calculated.
#' @param ... not used
#' 
#' @details
#' Residual correlation matrix is calculated based on the residual covariance matrix, see details from \code{\link{getResidualCov.gllvm}}.
#' 
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#' #'# Extract subset of the microbial data to be used as an example
#'data(microbialdata)
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, family = poisson())
#'fit$logL
#'cr <- getResidualCor(fit)
#'cr[1:5,1:5]
#'\dontrun{
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residual correlations:
#'cr <- getResidualCor(fit)
#'# Plot residual correlations:
#'install.packages("corrplot", "gclus")
#'library(corrplot)
#'library(gclus)
#'corrplot(cr[order.single(cr), order.single(cr)], diag = F,
#'   type = "lower", method = "square", tl.cex = 0.8, tl.srt = 45, tl.col = "red")
#'   }
#'
#'@aliases getResidualCor getResidualCor.gllvm
#'@method getResidualCor gllvm
#'@export
#'@export getResidualCor.gllvm
getResidualCor.gllvm = function(object, adjust = 1, site.index = NULL, ...)
{
  ResCov <- getResidualCov.gllvm(object, adjust = adjust, site.index = site.index)$cov
  Res.sd <- 1 / sqrt(diag(ResCov))
  Res.Cor <- diag(Res.sd) %*% ResCov %*% diag(Res.sd)
  colnames(Res.Cor) <- colnames(object$y)
  rownames(Res.Cor) <- colnames(object$y)
  Res.Cor[abs(Res.Cor) > 1] <- 1 * sign(Res.Cor[abs(Res.Cor) > 1])
  out <- Res.Cor
  return(out)
}

#'@export getResidualCor
getResidualCor <- function(object, adjust)
{
  UseMethod(generic = "getResidualCor")
}
