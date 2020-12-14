#' @title Plot simple 1D curves of a latent variable in a GLLVM
#' @description Plots latent variable and its corresponding coefficients.
#'
#' @param object   an object of class 'gllvm'.
#' @param which.spp  indices of species to plot.
#' @param main  main title.
#' @param which.lv index of latent variable to be plotted.
#' @param s.labels either "rug" or "numeric".
#' @param s.colors colors for sites
#' @param cex.spp size of species labels in biplot
#' @param spp.colors colors for sites, defaults to \code{"blue"}
#' @param xlim
#' @param ylim
#' @param ...	additional graphical arguments.
#'
#' 
#' @author Bert van der Veen
#'
#'@aliases uniplot uniplot.gllvm
#'@export
#'@export uniplot.gllvm
uniplot.gllvm <- function(object, main = NULL, xlim = NULL, s.labels = "rug", s.colors = "black", ylim = NULL, which.lv = 1, spp.colors = NULL,cex.spp=1,which.spp=NULL,...) {
  if (!any(class(object) %in% "gllvm"))
    stop("Class of the object isn't 'gllvm'.")
  n <- NROW(object$y)
  p <- NCOL(object$y)
  quadratic <- object$quadratic
  if(is.null(which.spp))which.spp<-1:p
  if(is.null(spp.colors)){
    spp.colors <- 1:p
    if(!is.null(which.spp)){
      spp.colors<-which.spp
    }
  }else if(length(spp.colors)==1){
    spp.colors <- rep(spp.colors,length(which.spp))
  }else if(length(spp.colors)!=length(which.spp)){
    stop("spp.colors needs to be of length p or 1.")
  }
  if (object$num.lv == 0)
    stop("No latent variables to plot.")
  
  if (is.null(rownames(object$params$theta)))
    rownames(object$params$theta) = paste("V", 1:p)
  
  if(quadratic==FALSE){
    func <- function(beta,x,u){
      return(beta+x*u)
    }
    }else{
      func<-function(x,c,u,d){
        optimum <- -u/(2*d)
        tolerance <- 1/sqrt(-2*d)
        return(c-optimum^2/(2*tolerance^2)+x*u+x^2*d)
      }
    }
  
  if(is.null(main)){
    main <- paste("1D plot of LV ",which.lv,sep="")
  }
  
  if(quadratic!=F){
    pred <- object$lvs[,which.lv,drop=F]%*%t(object$params$theta[,which.lv,drop=F])+object$lvs[,which.lv,drop=F]^2%*%t(object$params$theta[,which.lv+num.lv,drop=F])  
  }else{
    pred <- object$lvs[,which.lv,drop=F]%*%t(object$params$theta[,which.lv,drop=F])
  }
  maxima<-(object$params$beta0+rowSums(-0.5*object$params$theta[,which.lv]/(2*object$params$theta[,-c(1:num.lv)])*object$params$theta[,1:num.lv]))[which.spp]
  if(is.null(ylim)|length(ylim)!=2){
    if(max(maxima)>1000){
      maximum<-max(pred+cex.spp)
    }else{
      maximum<-max(maxima)
    }
    ylim <- c(min(pred[,which.spp]),maximum+cex.spp)
  }
  if(is.null(xlim)|length(xlim)!=2){
    xlim <- range(object$lvs[,which.lv])
    }
  
  plot(c(0,0),main=main,ylim=ylim,xlim=xlim,type="n",xlab=paste("LV",which.lv),ylab="Linear predictor")
  if (s.labels == "numeric") {
    text(x = object$lvs[, which.lv], y = -1, labels = 1:nrow(object$y), col = s.colors)
  } else if (s.labels == "rug") {
    rug(object$lvs[, which.lv], col = s.colors)
  }

  for(j in which.spp){
    if(quadratic!=FALSE){
      curve(func(x,c=maxima[which.spp==j], u=object$params$theta[j,which.lv],d=object$params$theta[j,which.lv+object$num.lv]),from=xlim[1],to=xlim[2],add=TRUE,col=spp.colors[which.spp==j])
      opt <- -object$params$theta[j,which.lv]/(2*object$params$theta[j,which.lv+num.lv])
      maximum <- 0.5*opt*object$params$theta[j,which.lv]
      if (opt < max(object$lvs[, which.lv]) & opt > min(object$lvs[, which.lv])) {
        text(x = opt, y = maximum, labels = colnames(object$y)[j], col = spp.colors[which.spp==j], cex = cex.spp, pos = 3) # should adjust "adj" rather than adding 0.1 to the maximum.
        #segments(x0 = opt, x1 = opt, y0 = min(pred[,j]), y1 = maximum, lty = "dashed", col = spp.colors[j])
      } else if (opt < min(object$lvs[, which.lv])) {
        text(x = min(object$lvs[, which.lv]), y = max(pred[,j]), labels = colnames(object$y)[j], col = spp.colors[which.spp==j], cex = cex.spp, pos = 4)
      } else if (opt > max(object$lvs[, which.lv])) {
        text(x = max(object$lvs[, which.lv]), y = max(pred[,j]), labels = colnames(object$y)[j], col = spp.colors[which.spp==j], cex = cex.spp, pos = 2)
      }
      
    }
    if(quadratic==FALSE){
      opt <- max(pred[,j])
      if(opt!=0){
        curve(func(x,beta=object$params$beta[j],u=object$params$theta[j,which.lv]),from=min(object$lvs[,which.lv]),to=max(object$lvs[,which.lv]),add=TRUE,col=spp.colors[j])
        x<- object$lvs[which(pred[,j]==opt), which.lv]
        text(x =x, y = opt, labels = colnames(object$y)[j], col = spp.colors[j], cex = cex.spp, pos = ifelse(x<0,4,2), offset = 1)  
      }
      
    
    }
    
    
    
  }
}




#'@export uniplot
uniplot <- function(object, ...)
{
  UseMethod(generic = "uniplot")
}

