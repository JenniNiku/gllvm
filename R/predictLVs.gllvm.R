#' @title Predict latent variables for gllvm Fits
#' @description Obtains predictions for latent variables from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param newX A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.
#' @param newY A new response data. Defaults to the dataset used for original model fit.
#' @param ... not used.
#'
#' @details
#' 
#' 
#' @return A matrix containing requested predictor types.
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- scale(antTraits$env[, 1:3])
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, family = poisson())
#'# fitted values
#'predLVs <- predictLVs.gllvm(fit, type = "response")
#'
#'@export
predictLVs.gllvm = function (object, newX = if(is.null(object$X)) NULL else object$X, newY=object$y, ...) 
{
  # create a gllvm object which refits the model using newX and newY:
  assign(as.character(object$call[2]),newY)
  if(is.null(newX)==FALSE)
    assign(as.character(object$call[3]),newX)
  objectTest=eval(object$call)

  # now optimise for LV parameters, keeping all else constant:
  objParam=object$TMBfn$par
  optLVs=optim(objParam,logL4lvs,gr=objectTest$TMBfn$gr,objectTest,object,method="BFGS")
  
  # now find the LV estimates to report:
  ui <- names(optLVs$par)=="u"
  n <- dim(objectTest$y)[1]
  lvs <- (matrix(optLVs$par[ui],n,objectTest$num.lv))
  rownames(lvs) <- rownames(objectTest$y);
  if(objectTest$num.lv>1) colnames(lvs) <- paste("LV", 1:objectTest$num.lv, sep="");

  # and their sds: (assuming method="VA")
  if(objectTest$num.lv>0)
  {
    Au <- optLVs$par[names(optLVs$par)=="Au"]
    A <- array(0,dim=c(n,objectTest$num.lv,objectTest$num.lv))
    for (d in 1:objectTest$num.lv)
    {
      for(i in 1:n)
      {
        A[i,d,d] <- exp(Au[(d-1)*n+i]);
      }
    }
    if(length(Au)>objectTest$num.lv*n)
    {
      k <- 0;
      for (c1 in 1:objectTest$num.lv)
      {
        r <- c1+1;
        while (r <=objectTest$num.lv)
        {
          for(i in 1:n)
          {
            A[i,r,c1] <- Au[objectTest$num.lv*n+k*n+i];
            A[i,c1,r] <- A[i,r,c1];
          }
          k <- k+1; r <- r+1;
        }
      }
    }
  }
  return(list(lvs=lvs,A=A))
}

logL4lvs = function(pars,objectTest,object)
{
  whichNotLVs = names(object$TMBfn$par)!="u" & names(object$TMBfn$par)!="Au"
  pars[whichNotLVs] = object$TMBfn$par[whichNotLVs]
  objectTest$TMBfn$fn(pars)
}

# to check:
# example(gllvm) # escape out of it after the second model, fitv0, has been fitted
# predictLVs.gllvm(fit,newX=fit$X,newY=simulate.gllvm(fit))
# predictLVs.gllvm(fitv0,newY=simulate.gllvm(fitv0))
