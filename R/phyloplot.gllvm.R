#' @title Plot phylogenetic random effects from gllvm
#' @description Plots phylogenetic random effects with the phylogeny, and community effects
#'
#' @param object an object of class 'gllvm'.
#' @param tree an object of class ''
#' @param comm.eff logical, defaults to \code{TRUE}. If present in the model, should community effects be plotted? 
#' @param row.eff logical, defaults to \code{FALSE}. If present in the model, should row effects (e.g., community responses to covariates) be included? 
#' @param which.Xcoef List of length 2 with names to subset the effects to plot. The first vector is for the species plot, the second for community effects.
#' @param xlim vector of length two. Limits for the x-axis of the caterpillar plot. Defaults to NULL, in which case the limits are chosen based on the confidence intervals.
#' @param level the confidence level. Scalar between 0 and 1.
#' @param col vector of three colors (defaults to \code{c("#E69F00","white","#009E73")}) passed to \code{\link{colorRampPalette}} for species random effects.
#' @param col.sym logical, defaults to \code{TRUE}. Then, the color scale of the species random effects plot is symmetrical (so that zero is nearly in the middle), so that both the lower and upper limit are determined by the largest absolute value. If \code{FALSE}, the lower and upper limits are determined by the smallest and largest values, respectively.
#' @param tick.length numerical, defaults to 2. Length of ticks for the species-specific random effects plot.
#' @param mar.spec vector of length 4, which defines the margins sizes for the species random effects plot. Defaults to \code{c(3, 2, 0, 0)}.
#' @param mar.phy vector of length 4, which defines the margins sizes for plotting the phylogeny. Defaults to \code{c(0, 2, 2, 0)}.
#' @param mar.comm vector of length 4, which defines the margins sizes for the caterpillar plot. Defaults to \code{c(3, 0.5, 2, 1.5)}.
#' @param cex the magnification to be used for text in the plot. Defaults to 0.6.
#' @param lwd line thickness for the branches in the phylogeny and the confidence intervals in the caterpillar plot. Defaults to 1.
#' @param col.edge character. Color of branches in the phylogeny.
#' @param pch symbol used in the catter pillar plot. Defaults to "x".
#' @param heights vector of length two. Relative row heights, defaults to \code{c(0.55, 0.35)}.
#' @param widths vector of length two. Relative column widths, defaults to \code{c(0.64, 0.10)}.
#' @param phy.place not (yet) in use.
#' @param ...	additional not in use.
#'
#' @details
#' Plots phylogenetically structured random effects together with the phylogeny,
#' and with community-level effects (i.e., effects that are the same across species).
#' If standard errors have been calculated for the model, the prediction intervals for
#' species random effects are checked, and crossed out (i.e., displayed as white) if
#' they cross zero.
#' 
#' @references 
#' van der Veen, B., O'Hara, R.B. (2024). Fast fitting of Fast fitting of phylogenetic mixed effects models. arXiv.
#' 
#' @author Bert van der Veen
#'
#' @examples
#' \dontrun{
#' # Load dataset
#' data(fungi)
#' Y <- fungi$Y
#' # Scale the predictor
#' X <- fungi$X
#' X[,"DBH.CM"] <- scale(X[, "DBH.CM"])
#' tree <- fungi$tree # the tree
#' colMat <- fungi$C # e.g., from ape::vcv(tree)
#' dist <- fungi$dist # e.g., from ape::cophenetic.phylo(tree)
#' order <- gllvm:::findOrder(covMat = colMat, distMat = dist, nn = 15,
#'                            order = order(dist[1:length(tree$tip.label), nrow(dist)],
#'                            decreasing = TRUE))$order
#' order <- tree$tip.label[order]
#' model <- gllvm(y = Y[,order], X = X,
#'                 formula = ~(DBH.CM|1), beta0com = TRUE,
#'                 family = "binomial", num.lv = 0, nn.colMat = 15,
#'                 colMat = list(colMat[order,order], dist = dist[order,order]), 
#'                 colMat.rho.struct = "term")
#' phyloplot(model, tree)
#'}
#'@aliases phyloplot phyloplot.gllvm
#'@export
#'@export phyloplot.gllvm

phyloplot.gllvm <- function(object, tree, comm.eff = TRUE, row.eff = FALSE, which.Xcoef = NULL, xlim = NULL, level = 0.95, col = c("#E69F00","white","#009E73"), col.sym = TRUE, tick.length = 2, mar.spec = c(3, 2, 0, 0), mar.phy = c(0, 2, 2, 0), mar.comm = c(3, 0.5, 2, 1.5), cex = 0.6, lwd = 1, col.edge = "black", pch = "x", heights = c(0.55, 0.35), widths = c(0.64, 0.1),  phy.place = "top", ...){
# add option to change the order of the plot
# graphical pars for every plot
  
  phy.place <- match.arg(phy.place, c("top", "right", "left"))
  if(phy.place %in% c("right", "left"))stop("Sorry, not implemented yet!")
  
  org.par = par(no.readonly = TRUE)
  if(inherits(object, "phylo"))stop("Did you swap the 'object' and 'tree' arguments?")
  if(!inherits(object, "gllvm"))stop("Model must be of class 'gllvm'.")
  if(!inherits(tree, "phylo"))stop("Tree must be of class 'phylo'.")
  if(is.null(object$params$Br))stop("No species random effects present in the model.")
  if(is.null(object$col.eff$colMat))stop("No phylogentic random effects in the model.")
  if((!is.null(which.Xcoef) && !is.list(which.Xcoef) || is.list(which.Xcoef) && length(which.Xcoef) != 2) && comm.eff)stop("'which.Xcoef' should be a list of length 2.")
  sd.err <- isFALSE(object$sd)||is.null(object$sd)
  if(sd.err && comm.eff)stop("No standard errors in the model.")
  
  # strike uncertain REs if possible
  if(sd.err){
  alpha = (1- level)/2
  PEs <- gllvm::getPredictErr(object)
  PIs <- data.frame(cbind(LI = c(object$params$Br)+c(PEs$Br*qnorm(alfa)), UI = c(object$params$Br)+c(PEs$Br*qnorm(1-alfa))))
  object$params$Br[-which((PIs$LI>0 & PIs$UI>0) | (PIs$LI<0 & PIs$UI<0))]<-NA
  }
  
  if(is.null(object$params$B) && (!row.eff || is.null(object$params$row.params.fixed)))comm.eff <- FALSE
  if(comm.eff){
    # Arrange 3 plots
  if(phy.place == "top"){
  layout(cbind(rbind(2,1),c(3,3)), respect = TRUE, heights = heights*dev.size()[2], widths = widths*dev.size()[1])
  
  # start with species REs  
  par(mar = mar.spec)
  if(is.character(which.Xcoef))which.Xcoef <- list(which.Xcoef)
  if(is.list(which.Xcoef) && !is.null(which.Xcoef[[1]])) {
    if(any(!which.Xcoef[[1]]%in%row.names(object$params$Br))){
      stop(paste0("Some species random effects names are not correctly written. This vector should include any of: ", paste0(row.names(object$params$Br), collapse=", "), "."))
    }
    object$params$Br <- object$params$Br[row.names(object$params$Br)%in%which.Xcoef[[1]],,drop=FALSE]
  }
  if(col.sym){
    breaks = seq(-max(abs(object$params$Br)), max(object$params$Br), by = 0.01)
  }else{
    breaks = seq(min(object$params$Br), max(object$params$Br), by = 0.01)
  }
  
  image(1L:ncol(object$y), 1:nrow(object$params$Br), t(object$params$Br[,tree$tip.label, drop=FALSE]), col = colorRampPalette(col)(length(breaks)-1), axes = FALSE, ylim = 0.5+c(0, nrow(object$params$Br)), xlim = 0.5+c(0,ncol(object$y)), xlab = NA, breaks = breaks, ylab = NA)
  mtext(side = 1, text = "Species-specific random effect", padj = 4, cex = cex)
  axis(2, at = 1:nrow(object$params$Br), labels = colnames(t(object$params$Br)), las = 1, cex.axis = cex, lwd = 0)
  par(xpd=TRUE)
  for(i in 1:nrow(object$params$Br)){
    lines(x = c(-tick.length,0), y = c(i,i))
  }
  
  par(mar = mar.phy)
  phylogram.gllvm(tree, direction = "down", cex = cex, scale = 0.95, col = col.edge, lwd = lwd, labels = FALSE)

  # caterpillar plot of community effects
  sds <- coefs <- NULL
  
  if(object$beta0com){
    coefs <- setNames(object$params$beta0[1], "Intercept")
    sds <- object$sd$beta0[1]
  }
  if(!is.null(object$params$B) && !is.null(object$TR)){
    coefs <- c(coefs, object$params$B[!grepl(":", names(object$params$B))])
    sds <- c(sds, object$sd$B[!grepl(":", names(object$params$B))])
  }else if(!is.null(object$params$B)){
    coefs <- c(coefs, object$params$B)
    sds <- c(sds, object$sd$B)
  }
  sds <- sds[sds>0]
  
  if(!is.null(object$row.params.fixed) && row.eff){
    coefs <- c(coefs, object$params$row.params.fixed)
    sds <- c(sds, object$sd$params$row.params.fixed)
  }
  
  CIs <- data.frame(cbind(LI=coefs+sds*qnorm(1-level), UI=coefs+sds*qnorm(level)))
  
  if(is.list(which.Xcoef) && !is.null(which.Xcoef[[2]])) {
    if(any(!which.Xcoef[[1]]%in%names(coefs))){
      stop(paste0("Some community effect names are not correctly written. This vector should include any of: ", paste0(names(coefs), collapse=", "), "."))
    }
    CIs <- CIs[names(coefs)%in%which.Xcoef[[2]],,drop=FALSE]
    coefs <- coefs[which.Xcoef[[2]]]
  }
  
  cols <- rep("black", length(coefs))
  cols[-which((CIs$LI>0 & CIs$UI>0) | (CIs$LI<0 & CIs$UI<0))] <- "grey"
  
  par(mar = mar.comm, xpd = FALSE)
  if(is.null(xlim)){
    xlim <- c(min(CIs$LI)-.1, max(CIs$UI)+.1)
  }
  plot(coefs, y = 1:length(coefs), yaxt = "n", ylab = "", pch = pch, col = cols, xlab = NA, cex.axis = cex, xlim = xlim)
  mtext(side = 1, text = "Community mean response", padj = 4, cex = cex)
  segments(x0 = CIs$LI, y0 = 1:length(coefs), x1 = CIs$UI, y1 = 1:length(coefs), col = cols, lwd = lwd)
  abline(v = 0, lty = "dashed")
  axis(4, at = 1:length(coefs), labels = names(coefs), las = 1, cex.axis = cex, lwd=0)
  }
  }else{
    # Arrange 2 plots
    if(phy.place == "top"){
        # need to reset the right margin on the defaults
        if(missing(mar.phy))mar.phy[4] <- 0.5
        if(missing(mar.spec))mar.spec[4] <- 0.5
      
      layout(rbind(2,1), respect = TRUE, heights = heights*dev.size()[2], widths = sum(widths)*dev.size()[1])
      
      # start with species REs  
      par(mar = mar.spec)
   
      if(is.character(which.Xcoef))which.Xcoef <- list(which.Xcoef)
      if(is.list(which.Xcoef) && !is.null(which.Xcoef[[1]])) {
        if(any(!which.Xcoef[[1]]%in%row.names(object$params$Br))){
          stop(paste0("Some species random effects names are not correctly written. This vector should include any of: ", paste0(row.names(object$params$Br), collapse=", "), "."))
        }
        object$params$Br <- object$params$Br[row.names(object$params$Br)%in%which.Xcoef[[1]],,drop=FALSE]
      }   
      breaks <- seq(min(object$params$Br), max(object$params$Br), by = diff(range(object$params$Br))/15) # arbitrary cut-off for colors
      image(1L:ncol(object$y), 1:nrow(object$params$Br), t(object$params$Br[,tree$tip.label,drop=FALSE]), col = colorRampPalette(col)(length(breaks)-1), axes = FALSE, ylim = 0.5+c(0, nrow(object$params$Br)), xlim = 0.5+c(0,ncol(object$y)), xlab = NA, breaks = breaks, ylab = NA)
      mtext(side = 1, text = "Species-specific random effect", padj = 4, cex = cex)
      axis(2, at = 1:nrow(object$params$Br), labels = colnames(t(object$params$Br)), las = 1, cex.axis = cex, lwd = 0)
      par(xpd=TRUE)
      for(i in 1:nrow(object$params$Br)){
        lines(x = c(-2,0), y = c(i,i))
      }
    
      par(mar = mar.phy)
      
      phylogram.gllvm(tree, direction = "down", cex = cex, scale = 0.95, col = col.edge, lwd = lwd, labels = FALSE)
    }   
    }
  
  on.exit(par(org.par), add = TRUE)
}

# This function plots the phylogenetic tree
# By recursively traversing the branches
# Starting at the tips, we draw segments
# The coordinate of the parent node is calculated as the mean of the coordinates of the child branches
phylogram.gllvm <- function(tree, direction = "right", scale = 0.95, col = "black", lwd = 1, labels = TRUE, ...){

  plot_branch <- function(node, xpos) {
    children <- which(edge[, 1] == node)
    
    if (length(children) == 0) {
      # Terminal node (tip)
      return(node_y[node])
    } else {
      # Internal node
      child_y <- numeric(length(children))
      for (i in seq_along(children)) {
        child <- edge[children[i], 2]
        x_next <- xpos + scaled_edge_length[children[i]]
        
        # Get the y-coordinate for the child node
        child_y[i] <- plot_branch(child, x_next)
        
        # Draw horizontal line to child
        if(direction %in% c("right", "left")){
          segments(xpos, child_y[i], x_next, child_y[i], col = col)  
        }else if(direction == "down"){
          segments(child_y[i], xpos, child_y[i], x_next, col = col, lwd = lwd)  
        }
        
      }
      # Assign and draw vertical line for internal node
      node_y[node] <- mean(child_y)
      
      if(direction %in% c("right", "left")){
        segments(xpos, min(child_y), xpos, max(child_y), col = col, lwd = lwd)  
      }else if(direction == "down"){
        segments(min(child_y), xpos, max(child_y), xpos, col = col, lwd = lwd)  
      }
      return(node_y[node])
    }
  }
  
  calculate_depth <- function(node, current_depth) {
    children <- which(edge[, 1] == node)
    node_depths[node] <<- current_depth
    
    for (child in children) {
      child_node <- edge[child, 2]
      branch_length <- edge.length[child]
      calculate_depth(child_node, current_depth + branch_length)
    }
  }
  
# Extract necessary information
edge <- tree$edge
edge.length <- tree$edge.length
tip.labels <- tree$tip.label

# Set y-coordinates for the tips
y_coords <- setNames(seq_along(tip.labels), tip.labels)

# Initialize y-coordinates for nodes
node_y <- numeric(max(tree$edge))
node_y[1:length(y_coords)] <- y_coords

root_node <- length(tree$tip.label) +1
node_depths <- numeric(nrow(tree$edge)+1)
calculate_depth(root_node, 0) # gets maximum depth

# Calculate the maximum cumulative branch length (root to farthest tip)
max_depth <- max(node_depths)
# Scale branch lengths to fill the plot width
scaled_edge_length <- edge.length / max_depth * scale  # Scaling factor for full width

# Initialize plot area with full width and tight y-limits
xlim = c(0,1)
ylim = c(0.5, length(tip.labels) + 0.5)
if(direction == "left"){
  xlim = c(1,0)
}
# swap x and y
if(direction == "down") {
  xlim = ylim
  ylim = c(1,0)
}
  
plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i", ...)
invisible(plot_branch(root_node, 0))

if(labels){
for (i in seq_along(tip.labels)) {
  # Place the label at the end of its respective branch
  text(scale + 0.01, y_coords[tip.labels[i]], labels = tip.labels[i], pos = 4, cex = 0.3)
}
}
}

#'@export phyloplot
phyloplot <- function(object, ...)
{
  UseMethod(generic = "phyloplot")
}