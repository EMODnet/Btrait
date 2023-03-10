
## ====================================================================
## Extracting trait data 
## ====================================================================

getTrait <- function(taxon,            # taxon name (vector)
                     trait,            # data.frame with traits
                     trait.class=NULL, # to obtain crisp values
                     trait.score=NULL,
                     taxonomy=NULL,    # data.frame with taxonomic data
                     t.column=1,       # columnname/nr with taxon in trait data.frame
                     standardize=FALSE, 
                     verbose=FALSE){ 
  
#  trait <- clearRows(trait, t.column, "trait")
  if (t.column == 0) {
    trait <- cbind(taxon=rownames(trait), taxon)
    t.column <- 1
  }
  cn    <- colnames(trait)
  all.x <- TRUE
  if (is.null(taxon))
    stop("'taxon' should be a vector with species names")
  if (all(is.na(taxon)))
    stop("'taxon' should be a vector with species names")
  if (is.null(trait))
    stop("'trait' should be a data.frame with species traits")
  if (! is.null(taxonomy)){  # expand trait with trait at higher taxonomic level
    trait <- rbind(trait, 
                   extendTrait(trait, taxonomy, t.column=t.column) )
  } 
  
  taxon <- unique(taxon)
  
  # check if all traits are already known
  NotMerged <- NULL
  if (! sum(!taxon %in% trait[,t.column]))
     trait.tax <- trait[trait[,t.column]%in% taxon ,]
  
  else {
    if (! is.null(taxonomy)){  # expand trait with traits at higher taxonomic level
      Btax  <- merge(data.frame(taxon=taxon), taxonomy, 
                   by.x = "taxon", 
                   by.y = 1)
      all.x <- FALSE
    } else Btax <- data.frame(taxon=taxon) 
    
    Btax$nr   <- 1:nrow(Btax)

    trait.tax <- merge(Btax, trait, 
                       by.x  = "taxon", 
                       by.y  = t.column, 
                       all.x = all.x) # merged by first taxon
    tname <- NULL
    
    if (! is.null (taxonomy)){  # find closest trait based on taxonomy
     taxnames <- colnames(taxonomy)
     tname <- taxnames[1]

     for (i in 2:ncol(taxonomy)){  # loop through genus, family, ... etc
      NotMerged <- Btax[ !Btax$nr %in% trait.tax$nr,]
      all.x     <- ifelse(i==ncol(taxonomy), TRUE, FALSE)
      Merged    <- merge(NotMerged, trait, 
                         by.x  = taxnames[i], 
                         by.y  = 1, 
                         all.x = all.x)
      trait.tax  <- rbind(trait.tax, Merged[,colnames(trait.tax)])
     }
    }
    
    trait.tax <- trait.tax[sort.int(trait.tax$nr, index.return=TRUE)$ix,]
    trait.tax <- trait.tax[,c("taxon",cn[-1])]
    if (! is.null(tname)) colnames(trait.tax)[1] <- tname
  }
  if (verbose & ! is.null(NotMerged)){
    if (iun <- nrow(NotMerged) > 0)
        warning( "for ", iun, " species, the traits could not be found - they are set = NA - see attributes()$notrait")
  }

  # restore ordering as in taxon
  row.names(trait.tax) <- trait.tax[,1]
  trait.tax            <- trait.tax[taxon,]
  row.names(trait.tax) <- NULL
  
  trait.tax[,1] <- taxon  # in case a taxon was not found
  
  notrait <- trait.tax[is.na(trait.tax[,2]),1]
  if (!is.null(trait.class))
       trait.tax <- fuzzy2crisp(
                     trait.tax, 
                     trait.class = trait.class, # to obtain crisp values
                     trait.score = trait.score, 
                     standardize = standardize)
  
  attributes(trait.tax)$notrait <- notrait
  trait.tax
}
