#' @title Extract residual correlations from gllvm object
#' @description  Calculates the residual correlation from 'gllvm' models that include latent variables.
#'
#' @param object   An object of class 'gllvm'
#'
#' @author David Warton
#'
#' @examples
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit GLLVM model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# residual correlations:
#'get.residual.cor.gllvm(fit)
#'}
#'@export
#'


get.residual.cor.gllvm = function(object)
{
  ResCov <- object$params$theta%*%t(object$params$theta)
  if (object$family == "negative.binomial")
    ResCov <- ResCov + diag(log(object$params$phi+1))
  Res.sd <- 1/sqrt(diag(ResCov))
  Res.Cor <- diag(Res.sd) %*% ResCov %*% diag(Res.sd)
  colnames(Res.Cor) <- colnames(object$y)
  rownames(Res.Cor) <- colnames(object$y)
  out <- list(cor=Res.Cor, cov=ResCov, trace=sum(diag(ResCov)))
  return(out)
}
