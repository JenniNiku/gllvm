#makeForm
# Extract correlation structure from formula 
corstruc<-function (term) 
{
  if (length(term) == 1) 
    if(c(term) %in% c("corAR1")){
      return("corAR1")
    } else if(c(term) == "corExp"){
      return("corExp")
    } else if(c(term) == "corCS"){
      return("corCS")
    } else {return("diag")}
  if (length(term) == 2) {
    if(c(term[[1]]) %in% c("corAR1")){
      # term <- corstruc(term[[2]])
      return("corAR1")
    } else if(c(term[[1]]) == "corExp"){
      return("corExp")
    } else if(c(term[[1]]) == "corCS"){
      return("corCS")
    } else if(term[[1]] == "("){return("diag")}
    else return(corstruc(term[[2]])) #term[[2]] <- corstruc(term[[2]])
    
  }
  stopifnot(length(term) >= 3)
  
  for (j in 2:length(term)) {
    term[[j]] <- corstruc(term[[j]])
  }
  term[[1]] = NULL
  unlist(term)
}

# Make random formula following lme4 codes
findbars1<-function (term) {
  fb <- function(term) {
    if (is.name(term) || !is.language(term)) 
      return(NULL)
    if (term[[1]] == as.name("(")) 
      return(fb(term[[2]]))
    stopifnot(is.call(term))
    if (term[[1]] == as.name("|")) 
      return(term)
    if (length(term) == 2) 
      return(fb(term[[2]]))
    c(fb(term[[2]]), fb(term[[3]]))
  }
  expandSlash <- function(bb) {
    makeInteraction <- function(x) {
      if (length(x) < 2) 
        return(x)
      trm1 <- makeInteraction(x[[1]])
      trm11 <- if (is.list(trm1)) 
        trm1[[1]]
      else trm1
      list(substitute(foo:bar, list(foo = x[[2]], bar = trm11)), 
           trm1)
    }
    slashTerms <- function(x) {
      if (!("/" %in% all.names(x))) 
        return(x)
      if (x[[1]] != as.name("/")) 
        stop("unparseable formula for grouping factor", 
             call. = FALSE)
      list(slashTerms(x[[2]]), slashTerms(x[[3]]))
    }
    if (!is.list(bb)) 
      expandSlash(list(bb))
    else unlist(lapply(bb, function(x) {
      if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
        lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | bar, list(foo = x[[2]], bar = trm)))
      else x
    }))
  }
  
  modterm <- expandDoubleVerts1(if (inherits(term,"formula")) 
  # modterm <- expandDoubleVerts1(if (is(term, "formula")) 
    term[[length(term)]]
    else term)
  expandSlash(fb(modterm))
}

subbars1<-function (term) 
{
  if (is.name(term) || !is.language(term)) 
    return(term)
  if (length(term) == 2) {
    if(c((term[[1]])) %in% c("corCS", "corAR1", "corExp"))
      term <- subbars1(term[[2]])
    else term[[2]] <- subbars1(term[[2]])
    return(term)
  }
  stopifnot(length(term) >= 3)
  if (is.call(term) && term[[1]] == as.name("|")) 
    term[[1]] <- as.name("+")
  if (is.call(term) && term[[1]] == as.name("||")) 
    term[[1]] <- as.name("+")
  for (j in 2:length(term)) {
    term[[j]] <- subbars1(term[[j]])
  }
  term
}


mkModMlist <- function (x, frloc) {
  frloc <- factorize(x, frloc)
  ff <- eval(substitute(factor(fac), list(fac = x[[3]])), frloc)
  nl <- length(levels(ff))
  mm <- model.matrix(eval(base::substitute(~foo, list(foo = x[[2]]))), 
                     frloc)
  
  sm <- Matrix::fac2sparse(ff, to = "d", drop.unused.levels = TRUE)
  sm <- Matrix::KhatriRao(sm, t(mm))
  dimnames(sm) <- list(rep(levels(ff), each = ncol(mm)), rownames(mm))
  list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}

#
mkReTrms1 <- function (bars, fr) 
{
  # drop.unused.levels = TRUE; 
  reorder.vars = FALSE
  if (!length(bars)) 
    stop("No random effects terms specified in formula", 
         call. = FALSE)
  stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(fr,"data.frame"))
  
  safeDeparse <- function(x) paste(deparse(x, 500L), collapse = " ")
  barnames <- function(bars) vapply(bars, function(x) safeDeparse(x[[3]]), "")
  
  names(bars) <- vapply(bars, function(x) paste(deparse(x[[3]], 500L), collapse = " "), "")
  term.names <- vapply(bars, safeDeparse, "")
  
  #
  blist <- lapply(bars, mkModMlist, fr) #drop.unused.levels, reorder.vars = reorder.vars)
  nl <- vapply(blist, `[[`, 0L, "nl")
  
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rbind, Ztlist)
  names(Ztlist) <- term.names
  
  ll <- list(Zt = Zt)
  ll
}


factorize <- function (x, frloc, char.only = FALSE) {
  for (i in all.vars(x[[length(x)]])) {
    if (!is.null(curf <- frloc[[i]])) 
      frloc[[i]] <- factor(curf)
  }
  return(frloc)
}


#
expandDoubleVerts1 <- function (term) {
  expandDoubleVert <- function(term) {
    frml <- formula(substitute(~x, list(x = term[[2]])))
    newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
    if (attr(terms(frml), "intercept") != 0) 
      newtrms <- c("1", newtrms)
    as.formula(paste("~(", paste(vapply(newtrms, function(trm) paste0(trm, 
                                                                      "|", deparse(term[[3]])), ""), collapse = ")+("), 
                     ")"))[[2]]
  }
  if (!is.name(term) && is.language(term)) {
    if (term[[1]] == as.name("(")) {
      term[[2]] <- expandDoubleVerts1(term[[2]])
    }
    stopifnot(is.call(term))
    if (term[[1]] == as.name("||")) 
      return(expandDoubleVert(term))
    term[[2]] <- expandDoubleVerts1(term[[2]])
    if (length(term) != 2) {
      if (length(term) == 3) 
        term[[3]] <- expandDoubleVerts1(term[[3]])
    }
  }
  term
}
