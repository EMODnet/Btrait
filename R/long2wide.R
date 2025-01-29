## ====================================================================
## ====================================================================
## Expand descriptor x taxon x density data in long format to wide format, 
## averaging over replicates
## ====================================================================
## ====================================================================

long2wide <- function(data,
                      row,             # descriptor (vector or data.frame)
                      column,          # vector with column names(e.g. taxon) 
                      value,           # density value (vector)
                      averageOver,     # replicate name (vector)
                      taxonomy=NULL,   # taxonomic relationships to 'taxon'
                      subset){         # logical expression which to include

  if (! missing(data)){
    if (! missing(row))
      row <- eval(substitute(data.frame(row)), envir = data, enclos = parent.frame())
    if (! missing(column))
      column <- eval(substitute(data.frame(column)), envir = data, enclos = parent.frame())
    if (! missing(value))
      value <- eval(substitute(data.frame(value)), envir = data, enclos = parent.frame())
    if (! missing(averageOver))
      averageOver <- eval(substitute(data.frame(averageOver)), envir = data, enclos = parent.frame())
    data_name <- substitute(data)
  } else data_name <- NA
  
  if (missing(averageOver))
    averageOver <- NULL
  
  taxon    <- column
  cnTaxon  <- getname(column)
  
  descriptor <- row
  cnDesc  <- getname(row)
  
# value can be missing
  if (is.null(value) | missing(value))
    value <- rep(1, times=length(column))
  
  else if (is.matrix(value) | is.data.frame(value)){
    if (ncol(value) > 1) stop ("'value' should be a vector for 'long2wide'")
    value <- as.vector(value)
  }
  
# value cannot have multiple columns
  if (is.matrix(value) | is.data.frame(value)){
    if (ncol(value) > 1) stop ("'value' should be a vector for 'long2wide'")
    value <- as.vector(value)
  }
  
  descriptor <- checkDescriptor(descriptor, length(value))
  
  # Check if all inputs are the same length
  if (is.list(taxon)) 
    taxon <- unlist(taxon)
  
  nr <- checklength(descriptor, taxon, value, averageOver)
    
  cnDesc   <- getname(descriptor)
  colnames(descriptor) <- cnDesc
  d.column <- 1:ncol(descriptor)

  cnValue  <- getname(value)
  cnames   <- c(cnDesc, cnTaxon, cnValue)

  # Combine multiple columns in descriptor and replicates to average
  Desc     <- OneFactor(descriptor) # factor each row of descriptor data.frame
  
  # number of replicates per descriptor  
  if (! is.null(averageOver)) {
     cnAvg <- colnames(averageOver)
     averageOver  <- OneFactor(as.data.frame(averageOver))
     nAvg <- aggregate(x   = averageOver, 
                       by  = list(desc=Desc), 
                       FUN = function(x) length(unique(x)))
     
  } else {
     averageOver <- rep(1, times=nr)
     nAvg <- data.frame(desc=Desc, x=rep(1, times=nr))
     cnAvg <- NA
     
  }

  # take subset 
  if (! missing(subset)) {
     x <- data.frame(descriptor, taxon=taxon, 
                     averageOver=averageOver, Desc=Desc, value)
     
     if (! missing(data))
       x <- data.frame(x, data)
     
     
     if (! is.null(taxonomy)) {
      nr <- nrow(x)
      x <- cbind(x, 1:nrow(x))
      nc <- ncol(x)
      x <- merge(x, taxonomy, by.x="taxon", by.y=1)
      if (nrow(x) != nr) 
        stop("cannot merge taxonomy with data: not all taxa present")
      x <- x[order(x[,nc]),]
      x <- x[,-nc]
     }
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r)) {
      r <- eval(subset, x, parent.frame(n=2))
    }
    
    if (!is.logical(r)) stop("'subset' must be logical")
    r <- r & !is.na(r)
    if (length(r) != nr) 
      stop ("'subset' evaluation did not provide a selection?")

    x <- x[r,]    # take subset
    
    # value can be a matrix/data.frame
    nV <- which(colnames(x) %in% cnValue)
    if (is.null(nV)) nV <- ncol(x)
    value <- x[,nV]

    taxon       <- x$taxon
    averageOver <- x$averageOver
    Desc        <- droplevels(x$Desc)
    
    # descriptor can have multiple columns
    nD <- which(colnames(x) %in% cnDesc)
    if (is.null(nD)) nD <- 1
    descriptor <- as.data.frame(x[,nD])  
    nr <- length(taxon)
   }
  
  nAvg <- nAvg[nAvg$desc %in% levels(Desc),2]
  
  iunique <- which(!duplicated(Desc))
  Stat   <- data.frame(Desc=Desc[iunique], descriptor[iunique,])
  names(Stat)[-1] <- cnDesc
  row.names(Stat) <- Stat$Desc                # to recover the ordering later
  Stat <- Stat[levels(Desc),]
  descriptor <- Desc

  taxon     <- factor(taxon)
  taxnames  <- levels(taxon)
#  descriptor <- factor(descriptor)

  ntaxa <- length(levels(taxon))
  nstat <- length(levels(descriptor))
  nout  <- as.integer(0)

# remove NAs in value
  if (any(is.na(value))){
    isok <- which(!is.na(value))
    averageOver<- averageOver[isok]
    value <- value[isok]
    taxon <- taxon[isok]
    if (is.matrix(descriptor) | is.data.frame(descriptor))
      descriptor <- descriptor[isok,]
    else
      descriptor <- descriptor[isok]
  }
  nr    <- length(value)
  value <- unlist(value)
  
  RES <- .Fortran("long2wide", ndat=length(value), ntaxa=ntaxa, nstat=nstat, 
                  as.integer(nAvg), desc=as.integer(descriptor),   
                  tax=as.integer(taxon), val=as.double(value), 
                  wide=matrix(0.0, nrow=nstat, ncol=ntaxa),
                  NAOK=FALSE, PACKAGE = "Btrait")
  value.Mean           <- RES$wide
  value.Mean           <- data.frame(Stat[,-1], value.Mean)
  colnames(value.Mean) <- c(cnDesc, taxnames)
  
  row.names(value.Mean) <- NULL
  attributes(value.Mean)$taxon.names      <- taxnames
  attributes(value.Mean)$dataset           <- data_name
  attributes(value.Mean)$names_row         <- cnDesc
  attributes(value.Mean)$names_column      <- cnTaxon
  attributes(value.Mean)$names_value       <- cnValue
  attributes(value.Mean)$names_averageOver <- cnAvg
  
  attributes(value.Mean)$d.column <- d.column
  
  return(value.Mean)
}

## ====================================================================

l2w_trait <- function(trait,        # data set with traits in long format          
                     descriptor,    # descriptor of sample -> columns 
                     taxon,         # vector with taxon names -> rows 
                     value,         # density value (vector)
                     averageOver,   # replicate name (vector)
                     taxonomy=NULL, # taxonomic relationships to 'taxon'
                     subset){
  if (! missing(trait)){
    if (! missing(descriptor))
      descriptor <- eval(substitute(data.frame(descriptor)), envir = trait, enclos = parent.frame())
    if (! missing(taxon))
      taxon <- eval(substitute(data.frame(taxon)), envir = trait, enclos = parent.frame())
    if (! missing(value))
      value <- eval(substitute(data.frame(value)), envir = trait, enclos = parent.frame())
    if (! missing(averageOver))
      averageOver <- eval(substitute(data.frame(averageOver)), envir = trait, enclos = parent.frame())
  }
  
  if (missing(averageOver))
    averageOver <- NULL
  
  if (missing(subset))
    long2wide(row         = taxon, 
              column      = descriptor, 
              value       = value, 
              averageOver = averageOver)
  else
    long2wide(row         = taxon, 
              column      = descriptor, 
              value       = value, 
              averageOver = averageOver,
              taxonomy    = taxonomy, 
              subset      = substitute(subset))
}

## ====================================================================

l2w_density <- function(data,
                       descriptor,     # descriptor of the sample -> columns
                       taxon,          # vector with taxon names  -> rows 
                       value,          # density value (vector)
                       averageOver, # replicate name (vector)
                       taxonomy=NULL,  # taxonomic relationships to 'taxon'
                       subset){   # logical expression which to include

  if (! missing(data)){
    if (! missing(descriptor))
      descriptor <- eval(substitute(data.frame(descriptor)), envir = data, enclos = parent.frame())
    if (! missing(taxon))
      taxon <- eval(substitute(data.frame(taxon)), envir = data, enclos = parent.frame())
    if (! missing(value))
      value <- eval(substitute(data.frame(value)), envir = data, enclos = parent.frame())
    if (! missing(averageOver))
      averageOver <- eval(substitute(data.frame(averageOver)), envir = data, enclos = parent.frame())
  }
  
  if (missing(averageOver))
    averageOver <- NULL
  
  if (missing(subset))
    long2wide(row         = descriptor, 
              column      = taxon, 
              value       = value, 
              averageOver = averageOver)

  else    
    long2wide(row         = descriptor, 
              column      = taxon, 
              value       = value, 
              averageOver = averageOver,
              taxonomy    = taxonomy, 
              subset      = substitute(subset))
}

## ====================================================================
## ====================================================================
## Go from wide to long format
## ====================================================================
## ====================================================================

wide2long <- function(wide, 
                      d.column = 1, # column name/nr with rownames
                      w.names  = NULL, 
                      absences = FALSE){       

  row.names.col <- d.column 

  Nonames <- FALSE
  if (is.null(row.names.col)) 
    Nonames <- TRUE
  else if (length(row.names.col) == 1)
    if (row.names.col == 0) Nonames <- TRUE
  
  if (Nonames) {  # use rownames as names
    wide <- cbind(descriptor=row.names(wide), wide)
    row.names.col <- 1
  }
  irow <- length(row.names.col)
                              
  # sum the same taxon within one descriptor and replicate 
  # (sometimes necessary due to renaming of some taxa)
  
  if (is.null(w.names)) {
    w.names <- data.frame(name=colnames(wide)[-row.names.col])
  }
  w.names <- as.data.frame(w.names)
  nw <- ncol(w.names)
  ii <- (irow+1):(irow+nw)  
  Long <- reshape(wide, 
                  direction ="long", 
                  varying=list(c(1:ncol(wide))[-row.names.col]), 
                  v.names="value")  [,-(irow+3)]
  LL <- data.frame(Long[,1:(irow)], w.names[Long[,irow+1],], 
                   Long[,(irow+2):ncol(Long)])
  names(LL)[1:irow] <- names(Long)[1:irow] #"taxon"
  names(LL)[ii] <- names(w.names)
  names(LL)[ncol(LL)] <- "value"
  if (! absences) 
    LL <- LL[which(LL[,ncol(LL)]!=0),]

  row.names(LL) <- NULL
  LL
}

## ====================================================================

w2l_density <- function(wide, d.column=1, taxon.names=NULL, absences=FALSE){
  if (is.null(taxon.names))
    taxon.names <- attributes(wide)$taxon.names

  wide2long(wide, d.column, w.names=taxon.names, absences=absences)
}

## ====================================================================

w2l_trait <- function(wide, t.column=1, trait.names=NULL, absences=FALSE){
  wide2long(wide, d.column=t.column, w.names=trait.names, absences=absences)
}

## ====================================================================
## ====================================================================
## Add 0s to density data in long format
## ====================================================================
## ====================================================================

add_absences <- function(data,
                        descriptor,      # descriptor name (vector) -> columns
                        taxon,           # vector with taxon names -> rows 
                        value,           # density value (vector)
                        averageOver){ # replicate name (vector)
 
  if (! missing(data)){
    if (! missing(descriptor))
      descriptor  <- eval(substitute(data.frame(descriptor)), 
                          envir = data, enclos = parent.frame())
    if (! missing(taxon))
      taxon       <- eval(substitute(data.frame(taxon)), 
                          envir = data, enclos = parent.frame())
    if (! missing(value))
      value       <- eval(substitute(data.frame(value)), 
                          envir = data, enclos = parent.frame())
    if (! missing(averageOver))
      averageOver <- eval(substitute(data.frame(averageOver)), 
                          envir = data, enclos = parent.frame())
  }
  
  if (missing(averageOver))
    averageOver <- NULL
  
  cnTaxon  <- getname(taxon)
  cnDesc   <- getname(descriptor)
  cnValue  <- getname(value)
  cnames   <- c(cnDesc, cnTaxon, cnValue)
  Z2 <- long2wide(row         = descriptor, 
                  column      = taxon, 
                  value       = value, 
                  averageOver = averageOver)

  ZZ <- wide2long(Z2, 
                  d.column=attributes(Z2)$d.column, absences=TRUE)
 
  rm("Z2")
  colnames(ZZ) <- cnames
  ZZ
}
