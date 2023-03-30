
## ====================================================================
## ====================================================================
## utilities for the Btrait package - Karline soetaert
## ====================================================================
## ====================================================================

# make datasets globally available

if (getRversion() >= "2.15.1")  
  utils::globalVariables(c("Traits_Db", "Traits_irr", "Traits_nioz", 
       "Traits_cefas", "Taxonomy", "Groups", "NSBS", "MWTL"))


## ====================================================================
## helper functions:
## ====================================================================

checkDescriptor <- function(descriptor, nr){
    if (is.null(descriptor) | all (is.na(descriptor))) 
     descriptor <- rep("NA", times=nr)
  
    if (length(unlist(descriptor))==1) 
     descriptor <- rep(descriptor, times=nr)
  
    if (inherits(descriptor, "list")) 
      descriptor <- data.frame(descriptor)

    if (!inherits(descriptor, "data.frame")) 
      descriptor <- as.data.frame(descriptor)
    return(descriptor)
}

## ====================================================================
## Check if all inputs are the same length
## ====================================================================

checklength <- function(descriptor,      # descriptor of the sample (vector or data.frame) -> columns
                        taxon,           # vector with taxon names -> rows 
                        value,           # density value (vector/data.frame/matrix)
                        averageOver=NULL){ # replicate name (vector)

  if (is.data.frame(taxon)) 
    nr <- nrow(taxon)
  else
    nr <- length(taxon)

  if (is.matrix(value) | is.data.frame(value)){
    if (nrow(value) != nr)   stop ("nrow of 'value' and 'taxon' not the same")
  } else if (length(unlist(value)) != nr) stop ("length of 'taxon' and 'value' not the same")
  
  if (! is.null(averageOver)){  
    if (is.matrix(averageOver) | is.data.frame(averageOver)){
       if (nrow(averageOver) != nr)   stop ("nrow of 'averageOver' and 'taxon' not the same")
    } else if (length(unlist(averageOver)) != nr)  stop ("length of 'taxon' and 'averageOver' not the same")
  }  
  
  if (is.matrix(descriptor) | is.data.frame(descriptor)){
    if (nrow(descriptor) != nr)   
       stop ("nrow of 'descriptor' and 'taxon' not the same")
  } else if (length(unlist(descriptor)) != nr) 
       stop ("length of 'descriptor' and 'taxon' not the same")
  invisible (nr)
}

## ====================================================================
## Creates one  factor for each row of a data.frame  
## ====================================================================

OneFactor <- function(frame){  # a data.frame
    frame[is.na(frame)] <- "NA"
    Desc <- as.factor(frame[,1])
    if (ncol(frame) > 1){
      for (i in 2:ncol(frame))
        Desc <- interaction(Desc, as.factor(frame[,i]), drop=TRUE)
    } 
    Desc
}

## ====================================================================
## Creates a data.frame with unique values from two columns
## ====================================================================
# tapply makes a list for each y with unique values of y
# stack concatenates them to a data.frame

unique.df <- function(x, y){
  stack ( tapply(INDEX=as.data.frame(x), 
                 X    =y, 
                 FUN  =function(x) unique(x)))[,2:1]
}

## ====================================================================
## Returns the name(s) of a variable
## ====================================================================

getname <- function(var) {
  cnValue <- NULL
  if (is.matrix(var))
    cnValue <- colnames(var)
  else if (is.list(var))
    cnValue <- names(var)
  if (is.null(cnValue))  
    cnValue <- deparse1(substitute(var))
  
  cnValue
}  

## ====================================================================
## Cleans up a dataframe by removing the column(s) with the rownames,
## and returns only the data; the names are cast as row.names
## ====================================================================

clearRows <- function(x, d.column, name){
  
  cn <- NULL # the name or nr of the column with the descriptor.names
  cnames <- colnames(x)
  
  # check input
  if (is.logical(d.column)) {
    if (d.column) d.column <- 1 else d.column <- NULL
  }
  if (length(d.column) && ! all(d.column == 0)) {
    if (is.character (d.column) & length(d.column) == 1){
      d.column <- which(colnames(x) == d.column)
      cn <- d.column
      if (! length(d.column))  stop("'d.column' not found in matrix ", name)
      rown    <- x[,d.column]
      x       <- x[,-(d.column)] 
      cnames  <- cnames[-d.column]
    } else if (is.character (d.column)){ 
       rown <- d.column  # 
    } else { # one or more numbers
      rown    <- x[,d.column]
      cn      <- colnames(x)[d.column]
      x       <- x[,-(d.column)] 
      cnames  <- cnames[-d.column]
    }  
  } else rown <- row.names(x)  # no separate column
  if (is.vector(x)) {
    x <- matrix(x, ncol=1)
    colnames(x) <- cnames
  }  
  rownNumeric <- is.numeric(rown)

  # if multiple columns, make one descriptor which is a concatenated string 
  if (!is.null(ncol(rown))) rown <- apply(rown, MARGIN=1, FUN=paste, collapse="_")
  if (length(rown) != nrow(x)) stop("'d.column' and matrix '", name, "' not compatible")
  row.names(x) <- rown
  attributes(x)$cn <- cn
  attributes(x)$colnames <- cnames
  attributes(x)$numericRownames <- rownNumeric
  return(x)
}
