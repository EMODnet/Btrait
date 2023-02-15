## ====================================================================
## ====================================================================
## Summary statistics
## ====================================================================
## ====================================================================

getSummary <- function(descriptor, taxon, value, averageOver = NULL, 
                       taxonomy = NULL, subset, what=c("density", "taxa", "occurrence"), 
                       wide.output = FALSE){
  
# value cannot have multiple columns
  if (is.matrix(value) | is.data.frame(value)){
    if (ncol(value) > 1) stop ("'value' should be a vector for 'getSummary'")
    value <- as.vector(value)
  }
  
  descriptor <- checkDescriptor(descriptor, length(taxon))

  # Check if all inputs are the same length
  if (is.list(taxon)) 
    taxon <- unlist(taxon)
  
  nr <- checklength(descriptor, taxon, value, averageOver)

  cnDesc     <- getname(descriptor)
  colnames(descriptor) <- cnDesc
  nc         <- ncol(descriptor)

   # take subset
  if (! missing(subset)) {
    if (! is.null(averageOver))
      x <- data.frame(descriptor, taxon=taxon, averageOver=averageOver, value)
    else
      x <- data.frame(descriptor, taxon=taxon, value)
     
     if (! is.null(taxonomy)) {
      nr <- nrow(x)
      x <- cbind(x, 1:nrow(x))
      nc <- ncol(x)
      x <- merge(x, taxonomy, by.x="taxon", by.y=1)
      if (nrow(x) != nr) stop("cannot merge taxonomy with data: not all taxa present")
      x <- x[order(x[,nc]),]
      x <- x[,-nc]
     }
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r)) stop("'subset' must be logical")
    r <- r & !is.na(r)
    if (length(r) != nr) stop ("'subset' evaluation did not provide a selection?")

    x <- x[r,]    # take subset
    
    # value can only be one number
    nV    <- ncol(x)
    value <- as.data.frame(x[,nV])

    taxon     <- x$taxon
    averageOver <- x$averageOver

    # descriptor can have multiple columns
    nD <- which(colnames(x) %in% cnDesc)
    if (is.null(nD)) nD <- 1
    descriptor <- as.data.frame(x[,nD])  
    nr <- length(taxon)
   }
    
   if (pmatch("den", what, nomatch = FALSE)) {
    density.total <- getDensity(
                         descriptor  = descriptor, 
                         averageOver = averageOver, 
                         taxon       = rep(1, times=nr), # all the same taxon
                         value       = value)
   
    density.total <- density.total[, -which(colnames(density.total)=="taxon")]
    colnames(density.total) <- c(cnDesc, "density")
  
  # second column removed
    if (nc > 1 & wide.output) 
      density.total <- long2wide(
        row    =density.total[,1:(nc-1)], 
        column =density.total[,nc], 
        value  =density.total$density)
    } else density.total <- NULL
  
   if (pmatch("tax", what, nomatch = FALSE) &  ! missing(taxon)) {
    # number of taxa
      density.numtax <- aggregate(taxon, 
                by  = as.data.frame(descriptor[,nc:1]),  # reverse the order to have same output as for density
                FUN = function(x) length(unique(x)))
      density.numtax[,1:nc] <-     density.numtax[,nc:1]  
      colnames(density.numtax) <-  c(cnDesc, "taxa")

      if (nc > 1 & wide.output) 
        density.numtax <- long2wide(
              row    = density.numtax[,1:(nc-1)], 
              column = density.numtax[,nc], 
              value  = density.numtax$taxa)
    } else 
    density.numtax <- NA

  # for all species: number of descriptors over which it is found. 

   if (pmatch("occ", what, nomatch = FALSE) &  ! missing(taxon)) {
    # Combine multiple columns in descriptor and replicates to average
    Desc   <- OneFactor(descriptor)   # returns factor for each row of descriptor (data.frame)
    species.numdescriptor <- aggregate(
               Desc, 
               by  = as.data.frame(taxon),  
               FUN = function(x) length(unique(x)))
    colnames(species.numdescriptor)[2] <-  "occurrence"
    
  } else 
    species.numdescriptor <- NA
  
# for species: mean density over all descriptors, 
   
  list(density    = density.total, 
       taxa       = density.numtax,
       occurrence = species.numdescriptor)
}

## ====================================================================
## Summary statistics: total (summed) values -- TO DO!!!!!!!!!!!
## ====================================================================

getSum <- function(descriptor, value, averageOver = NULL, 
                   subset, wide.output = FALSE){

# value cannot have multiple columns
  if (is.matrix(value) | is.data.frame(value)){
    if (ncol(value) > 1) stop ("'value' should be a vector for 'getTotal'")
    value <- as.vector(value)
  }
 
  descriptor <- checkDescriptor(descriptor, length(value))

  # Check if all inputs are the same length
  nr <- checklength(descriptor, value, value, averageOver)

  cnDesc     <- getname(descriptor)
  colnames(descriptor) <- cnDesc
  nc         <- ncol(descriptor)

   # take subset
  if (! missing(subset)) {
    if (! is.null(averageOver))
      x <- data.frame(descriptor, averageOver=averageOver, value)
    else
      x <- data.frame(descriptor, value)
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r)) stop("'subset' must be logical")
    r <- r & !is.na(r)
    if (length(r) != nr) stop ("'subset' evaluation did not provide a selection?")

    x <- x[r,]    # take subset
    
    # value can only be one number
    nV    <- ncol(x)
    value <- as.data.frame(x[,nV])
    averageOver <- x$averageOver

    # descriptor can have multiple columns
    nD <- which(colnames(x) %in% cnDesc)
    if (is.null(nD)) nD <- 1
    descriptor <- as.data.frame(x[,nD])  
    nr <- length(value)
   }
    
  density.total <- getDensity(
               descriptor  = descriptor, 
               averageOver = averageOver, 
               taxon       = rep(1, times=nr), # treat taxa as if all the same
               value       = value)
  density.total <- density.total[, -which(colnames(density.total)=="taxon")]
  colnames(density.total) <- c(cnDesc, "density")
  
  # second column removed
  if (nc > 1 & wide.output) 
    density.total <- long2wide(
        row    =density.total[,1:(nc-1)], 
        column =density.total[,nc], 
        value  =density.total$density)

  return(density.total)
}

## ====================================================================
## Summary statistics: weighted mean values -- TO DO!!!!!!!!!!!
## ====================================================================

getMean <- function(descriptor, value, weight=NULL,
                     subset, wide.output = FALSE){

# value cannot have multiple columns
  if (is.matrix(value) | is.data.frame(value)){
    if (ncol(value) > 1) stop ("'value' should be a vector for 'getTotal'")
    value <- as.vector(value)
  }
 
  descriptor <- checkDescriptor(descriptor, length(value))

  # Check if all inputs are the same length
  nr <- checklength(descriptor, value, weight, value)

  cnDesc     <- getname(descriptor)
  colnames(descriptor) <- cnDesc
  nc         <- ncol(descriptor)

   # take subset
  if (! missing(subset)) {
    x <- data.frame(descriptor, value, weight)
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r)) stop("'subset' must be logical")
    r <- r & !is.na(r)
    if (length(r) != nr) stop ("'subset' evaluation did not provide a selection?")

    x <- x[r,]    # take subset
    
    # value can only be one number
    nV    <- ncol(x)
    value <- as.data.frame(x[,nV])
    weight <- x$weight

    # descriptor can have multiple columns
    nD <- which(colnames(x) %in% cnDesc)
    if (is.null(nD)) nD <- 1
    descriptor <- as.data.frame(x[,nD])  
    nr <- length(value)
   }
    
  density.total <- getDensity(
               descriptor = descriptor, 
               taxon      = rep(1, times=nr), # treat taxa as if all the same
               value      = value)
  density.total <- density.total[, -which(colnames(density.total)=="taxon")]
  colnames(density.total) <- c(cnDesc, "density")
  
  # second column removed
  if (nc > 1 & wide.output) 
    density.total <- long2wide(
        row    =density.total[,1:(nc-1)], 
        column =density.total[,nc], 
        value  =density.total$density)

  return(density.total)
}

