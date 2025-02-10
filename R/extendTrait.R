## ====================================================================
## ====================================================================
## Expand taxon x trait data to include averages at higher levels
## ====================================================================
## ====================================================================

extend_trait <- function(trait,       # data.frame with traits; 1st column:taxon
                         taxonomy,    # data.frame with taxonomic data; 
                         taxon_column=1){ # nr or name of column in trait with taxon names
  
# convert/clean trait data  
  trait <- clearRows(trait, taxon_column, 'trait')  # taxon names become row.names
  cn    <- attributes(trait)$cn                 # name of column holding the taxa
  trnames <- c(cn, attributes(trait)$colnames)
  rn    <- row.names(trait)        # taxa for which traits are already assigned
  
  if (any (!is.numeric(na.omit(unlist(trait)))))
    stop("trait matrix should be numeric-convert categorical variables to fuzzy format")

# taxa converted to factors    
  LEVS   <- unique(c(rn, unlist(taxonomy)))  # levels from taxa (trait and taxonomy)
  TAX    <- sapply(as.data.frame(taxonomy),  # factorize all columns of taxonomy  
                   FUN=function(x) as.integer(factor(x, levels=LEVS)))
  TAX[is.na(TAX)] <- 0L   

# factors in trait names  
  Tnames <- as.integer(factor(rn, levels=LEVS))
  nmax   <- length(LEVS)-length(Tnames)+1    # levels in taxonomy not in trait
  if (nmax == 0) return(NULL)

# all values in trait converted to double - find NAs   
  trait <- sapply(trait, FUN=as.double)
  if (is.vector(trait)) trait <- matrix(ncol=1, data=trait)
  ZZ    <- apply(trait, MARGIN=1, FUN=function(x) sum(is.na(x)))
  isna  <- which (ZZ > 0)
  if (length(isna)) trait[isna,] <- 0
  
  if (is.null(nrow(TAX)) | is.null(nrow(trait))) stop ("cannot proceed: nothing to select from")
# calculate traits for taxa not in trait database  
  RES <- .Fortran("extendtraitf", 
         as.integer(nrow(TAX)), as.integer(ncol(TAX)), taxonomy=as.integer(TAX),
         as.integer(nrow(trait)), as.integer(ncol(trait)), trait=as.double(trait), 
         Tnames=Tnames, nmax=as.integer(nmax), nNew=0L,
         newTrait=matrix(nrow=nmax, ncol=ncol(trait), data=0.0),
         newName=rep(0L, times=nmax), numLevs=rep(0L,times=nmax), 
         Levs=rep(0L, times=nmax))
  
  nNew    <- RES$nNew  # number of new taxa whose traits have been determined
  if (nNew == 0) return(NULL)
  
  newName  <- LEVS[RES$newName[1:nNew]]

  newTrait <- data.frame(taxon=newName, RES$newTrait[1:nNew,])
  colnames(newTrait) <- trnames
  
  attributes(newTrait)$n     <- RES$numLevs[1:nNew]
  attributes(newTrait)$level <- colnames(taxonomy)[RES$Levs[1:nNew]]
  unique(newTrait)  # some may be double 
}
