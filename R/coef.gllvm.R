#'@export


coef.gllvm <- function(object, ...)
{
  names(object$params)[names(object$params)=="beta0"]="Intercept"
  if(object$row.eff %in% c(TRUE,"fixed")) names(object$params)[names(object$params)=="row.params"]="Row.Intercept"
  if(object$row.eff=="random") names(object$params)[names(object$params)=="row.params"]="Random.Row.Intercept"
  return(object$params)
}
