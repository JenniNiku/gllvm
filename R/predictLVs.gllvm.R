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
