# simulate function for gllvm objects.  Note this simulates marginally over LVs.
# an option is to add a conditional argument, which would fix the LVs.
# David Warton
# last modified 6th May 2019

simulate.gllvm = function (object, nsim = 1, seed = NULL, ...) 
{
  # code chunk from simulate.lm to sort out the seed thing:
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  nRows = dim(object$lvs)[1]
  nCols = dim(object$params$theta)[1]
  # generate new latent variables
  lvsNew = matrix(rnorm(nsim*nRows*object$num.lv),ncol=object$num.lv)
  if(is.null(object$X)) # for intercept models, there is a bug in predict fn so loop nsim times as a workaround
  {
     prs = matrix(NA,nsim*nRows,nCols)
     for(iSim in 1:nsim)
     {
       whichRows = (iSim-1)*nRows+1:nRows
       prTemp = predict.gllvm(object,newLV = lvsNew[whichRows,],type="response")
       prs[whichRows,] = prTemp
     }
     dimnames(prs) = list(1:(nsim*nRows),dimnames(prTemp)[[2]]) #match up column names
  }
  else
    prs = predict.gllvm(object,newX=object$X[rep(1:nRows,nsim),], newLV = lvsNew,type="response")
    
  # generate new data
  nTot = nsim*nRows*nCols # total number of values to generate
  if(object$family=="negative.binomial")
    invPhis = matrix(rep(object$params$inv.phi,each=nsim*nRows), ncol=nCols)
  newDat = switch(object$family,"binomial"=rbinom(nTot,size=1,prob=prs),
                  "poisson"=rpois(nTot,prs),
                  "negative.binomial"=rnbinom(nTot,size=invPhis,mu=prs),
                  stop(gettextf("family '%s' not implemented ", object$family), domain = NA))
  # reformat as data frame with the appropriate labels
  newDat = as.data.frame(matrix(newDat,ncol=nCols))
  dimnames(newDat)=dimnames(prs)
  return(newDat)
}
