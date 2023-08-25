#'@export


coef.gllvm <- function(object, parm = NULL, ...)
{
  if(is.null(parm)){
    parm = names(object$params)
  }
  pars <- object$params
  if(any(c("beta0","Intercept")%in%parm)){
    names(pars)[names(pars) == "beta0"] = "Intercept"
    parm[parm%in%c("beta0","Intercept")] <- "Intercept"
  }
  
  
  if (object$row.eff %in% c(TRUE, "fixed") && any(c("row.params","Row.Intercept") %in% parm)){
    names(pars)[names(pars) == "row.params"] = "Row.Intercept"
    parm[parm %in% c("row.params","Row.Intercept")] <- "Row.Intercept"
  }

  
  if (object$row.eff == "random" && any(c("row.params","Row.Intercept","Random.Row.Intercept") %in% parm)){
    names(pars)[names(pars) == "row.params"] = "Random.Row.Intercept"
    parm[parm %in% c("row.params","Row.Intercept","Random.Row.Intercept") ] <- "Random.Row.Intercept"
  }
    
    
  if((object$num.RR+object$num.lv.c)>0 && any(c("LvXcoef","Canonical.coefficients","Canonical.coef","Can.coef","Cancoef","Random.LvXcoef","Random.Canonical.coefficients","Random.Canonical.coef","Random.Can.coef","Random.Cancoef") %in% parm)){
    if(isFALSE(object$randomB)){
      names(pars)[names(pars) == "LvXcoef"] <- "Canonical.coefficients"
      parm[parm%in%c("LvXcoef","Canonical.coefficients","Canonical.coef","Can.coef","Cancoef","Random.LvXcoef","Random.Canonical.coefficients","Random.Canonical.coef","Random.Can.coef","Random.Cancoef")] <- "Canonical.coefficients"
    }
    if(!isFALSE(object$randomB)){
      names(pars)[names(pars) == "LvXcoef"] <- "Random.Canonical.coefficients"
      parm[parm%in%c("LvXcoef","Canonical.coefficients","Canonical.coef","Can.coef","Cancoef","Random.LvXcoef","Random.Canonical.coefficients","Random.Canonical.coef","Random.Can.coef","Random.Cancoef")] <- "Random.Canonical.coefficients"
      }
  }
  if((object$num.lv+object$num.lv.c+object$num.RR)>0 && any(c("theta","loadings","Species.scores") %in% parm)){
    names(pars)[names(pars) == "theta"] <- "Species.scores"
    parm[parm%in%c("theta","loadings","Species.scores")] <- "Species.scores"
    
  }
 if(length(pars[parm])==1){
   return(pars[[parm]])
 }else{
   return(pars[parm])
 }
}
