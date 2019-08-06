#' @title Predict latent variables for gllvm Fits
#' @description Obtains predictions for latent variables from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param newX A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.
#' @param newY A new response data. Defaults to the dataset used for original model fit.
#' @param ... not used.
#'
#' @details
#' Obtains predictions for latent variables from a fitted generalized linear latent variable model object.
#' 
#' @return A matrix containing requested predictor types.
#' @author  David Warton, Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#'  \donttest{
#' # Load a dataset from the mvabund package
#' data(antTraits)
#' y <- as.matrix(antTraits$abund)
#' X <- scale(antTraits$env[, 1:3])
#' # Fit gllvm model
#' fit <- gllvm(y = y, X, family = poisson())
#' # fitted values
#' predLVs <- predictLVs.gllvm(fit, type = "response")
#' }
#'@aliases predictLVs predictLVs.gllvm
#'@method predictLVs gllvm
#'@export
#'@export predictLVs.gllvm

predictLVs.gllvm <- function (object, newX = if(is.null(object$X)) NULL else object$X, newY=object$y, ...) 
{
  # create a gllvm object which refits the model using newX and newY:
  assign(as.character(object$call[2]),newY)
  if(is.null(newX)==FALSE)
    assign(as.character(object$call[3]),newX)
  objectTest=eval(object$call)

  # now optimise for LV parameters, keeping all else constant:
  whichLVs = names(objectTest$TMBfn$par)=="u" | names(objectTest$TMBfn$par)=="Au"
  # objParam = objectTest$TMBfn$env$last.par.best[whichLVs]
 objParam = objectTest$TMBfn$par[whichLVs]
#  objParam = objectTest$TMBfn$par
#  optLVs=optim(objParam,logL4lvs,gr=objectTest$TMBfn$gr,objectTest,object,method="BFGS")
  optLVs = optim(objParam,logL4lvs,gr=objParamGrad,objectTest,object,method="BFGS")
  
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
  return(list(lvs=lvs,A=A,logL=-optLVs$value))
}

# compute a negative logL using objectTest$TMBfn
logL4lvs = function(pars,objectTest,object)
{
  # build full parameter vector for objectTest$TMBfn, taking params from training fit:
  tpar = getPars(pars,objectTest,object)
  # compute new logL
  return(objectTest$TMBfn$fn(tpar))
}

# compute gradient of negative logL using objectTest$TMBfn
objParamGrad = function(pars,objectTest,object)
{ 
  # build full parameter vector for objectTest$TMBfn, taking params from training fit:
  tpar = getPars(pars,objectTest,object)
  # compute new gradient function
  gr=objectTest$TMBfn$gr(tpar)
  whichLVs = names(objectTest$TMBfn$par)=="u" | names(objectTest$TMBfn$par)=="Au"
  return(gr[whichLVs])
}

# function which builds up a parameter vector for use in objectTest$TMBfn, taking params from training fit (object)
# but setting row effects to mean value and leaving LVs free to estimate
getPars = function(pars,objectTest,object)
{
  # get parameters to use for prediction
  opar = object$TMBfn$par
  tpar = objectTest$TMBfn$par
  
  # set row effects and their variances to their mean value
  r0s=tpar[names(tpar)=="r0"]
  tpar[names(tpar)=="r0"] = mean(r0s)
  lg_Ars=tpar[names(tpar)=="lg_Ar"]
  tpar[names(tpar)=="lg_Ar"] = mean(lg_Ars)
  
  # replace LVs and their estimated se's with input values
  tpar[names(tpar)=="u"] = pars[names(pars)=="u"]
  tpar[names(tpar)=="Au"] = pars[names(pars)=="Au"]
  
  # replace other params with training values
  tpar[names(tpar)=="b"] = opar[names(opar)=="b"]
  tpar[names(tpar)=="B"] = opar[names(opar)=="B"]
  tpar[names(tpar)=="lambda"] = opar[names(opar)=="lambda"]
  tpar[names(tpar)=="lg_phi"] = opar[names(opar)=="lg_phi"]
  tpar[names(tpar)=="log_sigma"] = opar[names(opar)=="log_sigma"]
  return(tpar)
}


predictLogL.gllvm = function (object, newX = if(is.null(object$X)) NULL else object$X, newY=object$y, nRuns=1, ...) 
# function to give predictive log-likelihood for new observations, using parameters from a fitted gllvm object.
# Because TMB returns a sum of logL's, this function loops through the new dataset getting one prediction at a time
# on order to return logLs separately for each new observation
# these values are numerically unstable, try ramping up nRuns to get a better estimate
{
  # ensure newX and newY are dataframes (potherwise error on subsetting to rows)
  newX = data.frame(newX)
  newY = data.frame(newY)
  
  # initialise logLs
  nRows = dim(newY)[1]
  logLs = rep(NA,nRows)
  # get logL for each row using a for loop, attaching each new obs to original dataset:
  logLTry = rep(NA,nRuns)
  for(iRow in 1:nRows)
  {
    for(iRun in 1:nRuns)
    {
      logLTry[iRun] = predictLVs.gllvm(object, rbind(object$X,newX[iRow,]), rbind(object$y,newY[iRow,]), ...)$logL
    }
    logLs[iRow] = max(logLTry)
  }
  logLs = logLs - object$logL #subtract off logL from original dataset
  # and we are done
  return(logLs)
}

#' @export predictLVs
predictLVs <- function(object, ...)
{
  UseMethod(generic = "predictLVs")
}


# to check:
# example(gllvm) # escape out of it after the second model, fitv0, has been fitted
# predictLVs.gllvm(fit,newX=fit$X,newY=simulate(fit))
# predictLVs.gllvm(fitv0,newY=simulate.gllvm(fitv0))
# newy=simulate(fit)
# predictLVs.gllvm(fit,newX=fit$X[1:25,],newY=newy[1:25,])
# predictLogL.gllvm(fit,newX=fit$X[1:25,],newY=newy[1:25,])
