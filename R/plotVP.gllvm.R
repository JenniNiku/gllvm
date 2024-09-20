#' @title Plot variance partitioning
#' @description Function \code{plotVarPartitioning()} (alias \code{plotVP()}) plots the results of variance partitioning of a fitted gllvm.
#' 
#' @param VP 	a variance partitioning object for a gllvm produced by function varPartitioning.
#' @param main main title
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param legend.text a vector of names for the groups, as a default 'groupnames' from varPartitioning. If FALSE, legend not printed.
#' @param ... additional graphical arguments passed to the barplot function
#' 
#' 
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'@aliases plotVarPartitioning plotVP
#'@export
#'@export plotVarPartitioning
#'@rdname varPartitioning.gllvm 

plotVarPartitioning <- function(VP, main = "Variance Partitioning", xlab = "Response", ylab = "Variance proportion", legend.text = NULL, ...){
  arglist <- list(...)
  if(VP$family == "betaH"){
    par(mfrow=c(1,2))
    if(is.null(legend.text)){
      legend.textH <- paste(colnames(VP$PropExplainedVarHurdleSp), ", mean ", round(colMeans(VP$PropExplainedVarHurdleSp), digits = 3)*100, "%", sep = "")
    } else {legend.textH <- legend.text}
    if(is.null(legend.text)) legend.text <- paste(colnames(VP$PropExplainedVarSp), ", mean ", round(colMeans(VP$PropExplainedVarSp), digits = 3)*100, "%", sep = "")
    barplot(t(VP$PropExplainedVarSp), legend.text = legend.text, main = paste(main, ": Beta"), xlab = xlab, ylab = ylab, ...)
    barplot(t(VP$PropExplainedVarHurdleSp), legend.text = legend.textH, main = paste(main, ": Hurdle"), xlab = xlab, ylab = ylab, ...)
  } else {
    if(is.null(legend.text)) legend.text <- paste(colnames(VP$PropExplainedVarSp), ", mean ", round(colMeans(VP$PropExplainedVarSp), digits = 3)*100, "%", sep = "")
    barplot(t(VP$PropExplainedVarSp), legend.text = legend.text, main = main, xlab = xlab, ylab = ylab, ...)
  }
    
}

#'@export plotVP
#'@rdname varPartitioning.gllvm 
plotVP <- function(VP, ...)
{
  plotVarPartitioning(VP, ...)
}