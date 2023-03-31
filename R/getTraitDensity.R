
## ====================================================================
## ====================================================================
## Merging species density and trait data 
## ====================================================================
## ====================================================================

getTraitDensity <- function(
                      descriptor, taxon, value, 
                      averageOver=NULL, 
                      wide=NULL,         # density data, WIDE format (descriptor x taxon) 
                      d.column=1,        # nr/name of column with descriptor names in wide
                      trait,             # species x trait data, WIDE format    
                      t.column=1,        # nr/name of column with taxon names in trait
                      trait.class=NULL,  # indices to trait levels - vector
                      trait.score=NULL,  # indices to trait values - vector
                      taxonomy=NULL,     # if !NULL, trait is expanded at higher levels
                      scalewithvalue=TRUE,  # rescale with value - 
                      verbose=FALSE)  

{ 

# cast the density data in wide format, if not provided  
  if (is.null(wide)) {
    wide <- long2wide(row         = descriptor, 
                      column      = taxon, 
                      value       = value, 
                      averageOver = averageOver)
    d.column <- attributes(wide)$d.column
  }
  taxon.names <- attributes(wide)$taxon.names
  if (is.null(taxon.names))taxon.names <- colnames(wide)[-1]
  
  # KARLINE: CHECK IF THISIS WANTED...
#  if (! is.null(taxonomy)) {
#    taxonomy <- taxonomy[taxonomy[,1] %in% taxon.names,]
#    if (! nrow(taxonomy)) stop ("density data and taxonomy have nothing in common")
#  }
  x  <- clearRows(wide, d.column, "wide") # remove descriptor columns
  
  DESCs <- wide[,d.column]
  cn    <- attributes(x)$cn  # name(s) of the descriptor columns  
  if (is.null(taxon.names)) taxon.names <- colnames(x)

# check input of trait data
  trait <- clearRows(trait, t.column, "trait") # remove descriptor columns
  row.names.trait <- row.names(trait)
  
  # trait information for the taxa in the data
  trait <- getTrait (taxon    = taxon.names, 
                     trait    = data.frame(taxon=row.names.trait, trait), 
                     taxonomy = taxonomy)
  
  if (any (iun <- which(is.na(trait[,2]))))
    notrait <- trait$taxon[iun]  # species not in trait database
  else 
    notrait <- NA
  
  if (verbose & length(iun)) 
    warning( length(iun), " species are not in the trait database - check attributes()$notrait to find them")
    
  row.names.trait <- trait[,1]
  tnames <- colnames(trait)[-1]
  trait  <- trait[,-1]
  if (is.vector(trait)) {
    trait <- matrix(trait, ncol=1)
    colnames(trait) <- tnames
  }  
  if (nrow(trait) != ncol(x)) {
    stop("dimensions of 'wide' (", nrow(x), ",", ncol(x), ") and 'trait' (", 
        nrow(trait), ",", ncol(trait), ") not compatible- check input or rownames?")
  }
  if (! is.null (row.names.trait) & ! is.null(taxon.names))
    if (!identical(row.names.trait,  taxon.names)) 
      stop("names of 'wide' and 'trait' not compatible")

  NSP             <- as.matrix(trait)
  NSP[is.na(NSP)] <- 0
  X               <- as.matrix(x)
  X[is.na(X)]     <- 0
  
  if (any(!is.numeric(NSP)))
    stop("trait matrix should be numeric - convert categorical variables to fuzzy format")
  cwm <- as.matrix(X) %*% NSP   # density of all traits
  cwm <- data.frame(DESCs, cwm)
  colnames(cwm)[1:length(cn)] <- cn
  
  if (!is.null(trait.class))
    cwm <- fuzzy2crisp(cwm,       # species-trait matrix - in WIDE format
                       t.column    = d.column, 
                       trait.class = trait.class,
                       trait.score = trait.score, # indices to trait values - vector
                       standardize = FALSE)
  if (scalewithvalue) {  
    if (length(iun)) X[, iun] <- 0
    
    if (is.vector(cwm[, -d.column]))
      cwm[, -d.column] <- cwm[,-d.column]/rowSums(X)  # d.column/taxon?
    else 
      cwm[, -d.column]  <- sweep(cwm[,-d.column], 
                                 MARGIN = 1, 
                                 STATS  = rowSums(X, na.rm=TRUE),
                                 FUN    = "/")  # d.column/taxon?
  }
  
  row.names(cwm) <- NULL
  cwm <- as.data.frame(cwm)
  attributes(cwm)$notrait <- notrait
  cwm
}  
