#' @title Extract residual correlations from gllvm object
#' @description  Calculates the residual correlation matrix for gllvm model.
#'
#' @param object   an object of class 'gllvm'.
#' @param adjust  defaults to 1, when residual covariance is adjusted in the case of negative binomial and binomial distribution. Alternatives are 0, when adjustment is not used, and 2 for negative binomial distribution. See function \code{\link{getResidualCov.gllvm}}.
#' 
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residual correlations:
#'cr <- getResidualCor(fit)
#'\dontrun{
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
getResidualCor.gllvm = function(object, adjust = 1)
{
  ResCov <- getResidualCov.gllvm(object, adjust = adjust)$cov
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
