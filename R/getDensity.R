## ====================================================================
## ====================================================================
## Takes averages from density data
## ====================================================================
## ====================================================================

getDensity <- function(descriptor,      # descriptor of sample, (vector or data.frame) 
                       taxon,           # vector with taxon names 
                       value,           # density value (vector or matrix/data.frame)
                       averageOver=NULL,# replicate name (vector or data.frame)
                       taxonomy=NULL,   # taxonomic relationships to 'taxon'
                       subset,          # logical expression which to include
                       wide.output=FALSE,  # if true, recasts result in wide format
                       full.output=FALSE,  
                       verbose=FALSE  ){ # if true: adds descriptors with 0 value
  
  # value can be a list...
  if (!is.data.frame(value)) 
    value <- as.data.frame(value)
 
 
  nr <- nrow(value)
  
  if (missing(taxon)) 
    taxon <- rep("NA", times=nr)
  cnTaxon <- getname(taxon)
  if (is.list(taxon)) 
    taxon <- unlist(taxon)

  if (missing(descriptor)) 
    descriptor <- rep("NA", times=nr)

  descriptor <- checkDescriptor(descriptor, nr)

  # Check if all inputs are the same length
  nr <- checklength(descriptor, taxon, value, averageOver)
    
  cnDesc  <- getname(descriptor)
  colnames(descriptor) <- cnDesc
  
  cnValue <- getname(value)
  cnames  <- c(cnDesc, cnTaxon, cnValue)

  if (! is.numeric(value)) {
     value <- as.data.frame(lapply(value, as.numeric))
     if (verbose) warning("'value' should be numeric")
  }

  # Combine multiple columns in descriptor and replicates to average
  Desc <- OneFactor(descriptor)   # factor each row of descriptor (data.frame)
  
  # descriptors for which output needs to be returned - defined BEFORE subsetting
  if (full.output) {
    Desc.ori <- levels(Desc)  
    iunique  <- which(!duplicated(Desc))
    Stat     <- data.frame(Desc=Desc[iunique], descriptor[iunique,])
    names(Stat)[-1] <- cnDesc
    row.names(Stat) <- Stat$Desc        # to recover the ordering later
  }
  
  # number of replicates per descriptor  
  if (! is.null(averageOver)) {
     averageOver  <- OneFactor(as.data.frame(averageOver))
     nAvg <- aggregate(x   = averageOver, 
                       by  = list(desc=Desc), 
                       FUN = function(x) length(unique(x)))
  } else { # one replicate
     averageOver <- rep(1, times=nr)
     nAvg <- data.frame(desc=Desc, x=rep(1, times=nr))
  }

  # take subset
  if (! missing(subset)) {
     x <- data.frame(descriptor, taxon=taxon, 
               averageOver=averageOver, Desc=Desc, value)
     
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
    if (sum(r) == 0) stop ("'subset' evaluation did not provide a selection?")

    x <- x[r,]    # take subset
    
    # value can be a matrix/data.frame
    nV <- which(colnames(x) %in% cnValue)
    if (is.null(nV)) nV <- ncol(x)
    value <- as.data.frame(x[,nV])

    taxon     <- x$taxon
    averageOver <- x$averageOver
    Desc      <- droplevels(x$Desc)
    
    # descriptor can have multiple columns
    nD <- which(colnames(x) %in% cnDesc)
    if (is.null(nD)) nD <- 1
    descriptor <- as.data.frame(x[,nD])  
    nr <- length(taxon)
   }

   # adapt nAvg to the subset
   nAvg <- nAvg[nAvg$desc %in% levels(Desc),2]

   # descriptors for which output needs to be returned - defined AFTER subsetting
   if (!full.output) {
    iunique <- which(!duplicated(Desc))
    Stat   <- data.frame(Desc=Desc[iunique], descriptor[iunique,])
    names(Stat)[-1] <- cnDesc
    row.names(Stat) <- Stat$Desc              # to recover the ordering later
   }
   Stat       <- Stat[sort.int(Stat[,1], index.return=TRUE)$ix,]
   descriptor <- Desc

# fortran routine requires sorted input consistent with Desc
   ii         <- sort.int(descriptor, index.return=TRUE)$ix
   descriptor <- descriptor[ii]
   descriptor <- factor(descriptor)  

   averageOver  <- averageOver[ii]

   taxon[is.na(taxon)] <- "UNKNOWN TAXON"
   taxon <- factor(unlist(taxon)[ii])

   if (is.matrix(value) | is.data.frame(value)){
    value <- data.frame(value[ii,])
   } else value <- data.frame(value[ii])

   AVG <- unique.df(descriptor, as.character(taxon))
   colnames(AVG) <- c("descriptor", "taxon")

   for (i in 1:ncol(value)){
    NEW <- density.select.list(descriptor, taxon, unlist(value[,i]), 
                               averageOver, nAvg)
    colnames(NEW)[3] <- paste("value",i, sep="")  # merge requires different column names
    AVG <- merge(AVG, NEW, all.x=TRUE, by=1:2)
   }
  
   if (full.output) {
     nc  <- ncol(AVG)
     DOR <- expand.grid(levels(droplevels(taxon)),Desc.ori)
     AVG <- merge(DOR[,2:1], AVG, by=1:2, all.x=TRUE)
     AVG[is.na(AVG[, nc]), nc] <- 0
   }
   
   # replace "pasted" descriptor names by original columns
   # and restore ordering in AVG
   
   st      <- AVG[,1]
   stNames <- data.frame(Stat[st,-1])
   AVG     <- data.frame(stNames, AVG[,-1])

   colnames(AVG) <- cnames
   row.names(AVG) <- NULL
   if (wide.output) {  #convert to wide format
     if (length(cnValue) > 1) stop ("cannot convert to wide format if there is more than one value")
     ndesc <- length(cnDesc)  # number of columns of descriptors
     AVG <- long2wide(
        row      =AVG[,1:ndesc], 
        column   =AVG[,(ndesc+1)], 
        value    =AVG[,ncol(AVG)])
     colnames(AVG)[1:ndesc] <- cnDesc
   }
   AVG
}

## ====================================================================
## wrapper around fortran routine - one set of values
## ====================================================================

density.select.list <- function(descriptor, taxon, value, averageOver, nAvg){

# all inputs have been checked and are factors with the same length

# remove NAs in value
  if (any(is.na(value))){
    isok       <- which(!is.na(value))
    averageOver<- averageOver[isok]
    value      <- value[isok]
    taxon      <- taxon[isok]
    descriptor <- descriptor[isok]
  }
  ntaxa <- length(levels(taxon))
  nstat <- length(levels(descriptor))
  nout  <- as.integer(0)
  
# summary taking the SUM within descriptor x averageOver 
#         and then average over nAverage
  RES <- .Fortran("long2long", ndat=length(value), ntaxa=ntaxa, nstat=nstat, as.integer(nAvg),   
                  desc=as.integer(descriptor), tax=as.integer(taxon), 
                  val=as.double(value), wide=rep(0.0, times=ntaxa),
                  nout=nout, NAOK=FALSE, PACKAGE = "Btrait")
#  RES <- .Fortran("long2longMean", ndat=length(value), ntaxa=ntaxa, nstat=nstat, as.integer(nAvg),   
#                  desc=as.integer(descriptor), tax=as.integer(taxon), 
#                  as.integer(averageOver), val=as.double(value), wide=rep(0.0, times=ntaxa),
#                  nout=nout, NAOK=FALSE, PACKAGE = "Btrait")
  nout <- RES$nout
  AVG  <- data.frame(descriptor=levels(descriptor)[RES$desc[1:nout]], 
                     taxon     =levels(taxon)     [RES$tax[1:nout]], 
                     value     =RES$val[1:nout])
  AVG # 3 columns
}

## ====================================================================
## ====================================================================
## Estimates proportions from density data
## ====================================================================
## ====================================================================

 getProportion <- function (descriptor, taxon, value, 
                            averageOver = NULL, 
                            taxonomy    = NULL, 
                            verbose     = FALSE) {

   if (!is.data.frame(value))                                           
       value <- as.data.frame(value)                                    
   
   nr <- nrow(value)
   if (ncol(value) > 1)
     stop("cannot calculate proportions for more than one value column")
   
   descriptor <- checkDescriptor(descriptor, nr)    

#   if (is.list(taxon)) 
#    taxon <- unlist(taxon)
      
#   nr         <- checklength(descriptor, taxon, value, averageOver)             
   cnDesc     <- getname(descriptor) 
   NN         <- length(cnDesc)
      
   Mean <- getDensity(descriptor  = descriptor, 
                      value       = value, 
                      taxon       = taxon,
                      averageOver = averageOver, 
                      taxonomy    = taxonomy, 
                      wide.output = FALSE, 
                      verbose     = verbose)
   cn <- ncol(Mean)
 
 # total density             
   Desc  <- Mean[,1:NN]
   Total <- aggregate(Mean[,cn], by =data.frame(Desc), FUN=sum) 
   
 # relative density (fraction)           
   Mean   <- merge(Mean, Total, by=1:NN) # add total density
   Mean$p <- Mean[,cn]/Mean[,ncol(Mean)]        # estimate fraction
   Mean   <- Mean[,-(cn:(cn+1))]
   Mean
 }
