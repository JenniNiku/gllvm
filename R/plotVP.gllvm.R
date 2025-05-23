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
#' @param ... additional graphical arguments passed to the barplot function
#' 
#' 
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'@aliases plotVarPartitioning plotVP plot.VP.gllvm
#'@export
#'@export plotVarPartitioning
#'@rdname VP.gllvm 

plotVarPartitioning <- function(x, main = "Variance Partitioning", xlab = "Response", ylab = "Variance proportion", legend.text = NULL, args.legend =list(cex=0.7, x="topright", bty="n", inset =c(0,-0.15)), mar = c(4,4,6,2), ...){
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
  if(x$family == "betaH"){
    par(mfrow=c(1,2), mar = mar)
    if(is.null(legend.text)){
      legend.textH <- paste(colnames(x$PropExplainedVarHurdleSp), ", mean ", round(colMeans(x$PropExplainedVarHurdleSp), digits = 3)*100, "%", sep = "")
    } else {legend.textH <- legend.text}
    if(is.null(legend.text)) legend.text <- paste(colnames(x$PropExplainedVarSp), ", mean ", round(colMeans(x$PropExplainedVarSp), digits = 3)*100, "%", sep = "")
    barplot(t(x$PropExplainedVarSp), legend.text = legend.text, main = paste(main, ": Beta"), xlab = xlab, ylab = ylab, args.legend = args.legend, ...)
    barplot(t(x$PropExplainedVarHurdleSp), legend.text = legend.textH, main = paste(main, ": Hurdle"), xlab = xlab, ylab = ylab, args.legend = args.legend, ...)
  } else {
    par(mar = mar)
    if(is.null(legend.text)) legend.text <- paste(colnames(x$PropExplainedVarSp), ", mean ", round(colMeans(x$PropExplainedVarSp), digits = 3)*100, "%", sep = "")
    barplot(t(x$PropExplainedVarSp), legend.text = legend.text, main = main, xlab = xlab, ylab = ylab, args.legend = args.legend, ...)
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

