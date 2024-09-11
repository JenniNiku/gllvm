#' @title Plot variance partitioning
#' @description Plots the results of variance partitioning of a fitted gllvm.
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
#' @examples
#' \dontrun{
#'# Calculate variance partitioning for a fitted model 'fit'
#'VP <- varPartitioning(fit)
#'# Plot the result of  variance partitioning
#'plot(VP, col = palette(hcl.colors(5, "viridis")))
#'}
#'@aliases plotVarPartitioning plotVP
#'@export
#'@export plotVarPartitioning

plotVarPartitioning <- function(VP, main = "Variance Partitioning", xlab = "Response", ylab = "Variance proportion", legend.text = NULL, ...){
  arglist <- list(...)
  
  if(is.null(legend.text)) legend.text <- paste(colnames(VP$PropExplainedVarSp), ", mean ", round(colMeans(VP$PropExplainedVarSp), digits = 3)*100, "%", sep = "")
  barplot(t(VP$PropExplainedVarSp), legend.text = legend.text, main = main, xlab = xlab, ylab = ylab, ...)
}

#'@export plotVP
plotVP <- function(VP, ...)
{
  plotVarPartitioning(VP, ...)
}