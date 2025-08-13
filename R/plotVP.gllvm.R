#' @title Plot variance partitioning
#' @description Function \code{plotVarPartitioning()} (alias \code{plotVP()} or just \code{plot()}) plots the results of variance partitioning of a fitted gllvm.
#' 
#' @param x 	a variance partitioning object for a gllvm produced by function varPartitioning.
#' @param main main title
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param legend.text a vector of names for the groups, as a default 'groupnames' from varPartitioning. If FALSE, legend not printed.
#' @param args.legend a list of additional arguments to pass to \code{legend()}.
#' @param mar Margins of the plot. Default \code{c(4,4,6,2)}
#' @param r2scaled logical, whether or not to plot the response specific R2 -scaled variance partitioning. Defaults to FALSE.
#' @param ... additional graphical arguments passed to the barplot function
#' 
#' 
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'@aliases plotVarPartitioning plotVP plot.VP.gllvm
#'@export
#'@export plotVarPartitioning
#'@rdname VP.gllvm 

plotVarPartitioning <- function(x, main = "Variance Partitioning", xlab = "Response", ylab = "Variance proportion", legend.text = NULL, args.legend =list(cex=0.7, x="topright", bty="n", inset =c(0,-0.15)), mar = c(4,4,6,2), r2scaled = FALSE, ...){
  arglist <- list(...)
  fill.args.legend <- function(l1){
    if(!("cex" %in% names(l1))){
      l1$cex = 0.7
    }
    if(!("x" %in% names(l1))){
      l1$x = "topright"
    }
    if(!("bty" %in% names(l1))){
      l1$bty = "n"
    }
    if(!("inset" %in% names(l1))){
      l1$inset = c(0,-0.15)
    }
    
  }

  scaler2 <- rep(1, nrow(x$PropExplainedVarSp))
  
  if(r2scaled){
    if(x$family=="ordinal")
      stop("At the moment, we do not have proper response specific R-squared measure for ordinal model inplemented.")
    if(!("r2species" %in% names(x)))
      stop("Scaled variance partitioning is missing from the object. Please re-calculate variance partitioning with an option 'r2scaled = TRUE'")
    scaler2 <- x$r2species
    if(!("main" %in% names(arglist))) main = "r2-scaled Variance Partitioning"
    if(!("ylab" %in% names(arglist))) ylab = "r2-scaled variance proportion"
  }
  if(x$family == "betaH"){
    scaler2H <- rep(1, nrow(x$PropExplainedVarHurdleSp))
    if(r2scaled){
      scaler2H <- x$r2Hspecies
    }
    par(mfrow=c(1,2), mar = mar)
    if(is.null(legend.text)){
      legend.textH <- paste(colnames(x$PropExplainedVarHurdleSp), ", mean ", round(colMeans(x$PropExplainedVarHurdleSp*scaler2H, na.rm = TRUE), digits = 3)*100, "%", sep = "")
    } else {legend.textH <- legend.text}
    if(is.null(legend.text)) legend.text <- paste(colnames(x$PropExplainedVarSp), ", mean ", round(colMeans(x$PropExplainedVarSp*scaler2, na.rm = TRUE), digits = 3)*100, "%", sep = "")
    barplot(t(x$PropExplainedVarSp*scaler2), legend.text = legend.text, main = paste(main, ": Beta"), xlab = xlab, ylab = ylab, args.legend = args.legend, ...)
    barplot(t(x$PropExplainedVarHurdleSp*scaler2H), legend.text = legend.textH, main = paste(main, ": Hurdle"), xlab = xlab, ylab = ylab, args.legend = args.legend, ...)
  } else {
    par(mar = mar)
    if(is.null(legend.text)) legend.text <- paste(colnames(x$PropExplainedVarSp), ", mean ", round(colMeans(x$PropExplainedVarSp*scaler2), digits = 3)*100, "%", sep = "")
    barplot(t(x$PropExplainedVarSp*scaler2), legend.text = legend.text, main = main, xlab = xlab, ylab = ylab, args.legend = args.legend, ...)
  }
    
}

#'@export plotVP
#'@rdname VP.gllvm 
plotVP <- function(x, ...)
{
  plotVarPartitioning(x, ...)
}

#'@export 
#'@rdname VP.gllvm 
plot.VP.gllvm <- function(x, ...)
{
  plotVarPartitioning(x, ...)
}

