## ====================================================================
## ====================================================================
## Eztracts information from datasets
## ====================================================================
## ====================================================================

## ====================================================================
## Get the description of a database
## ====================================================================

metadata <- function(object){
  attributes(object)$description
}

## ====================================================================
## Traits of a taxon
## ====================================================================

get_taxon_info <- function(taxon, trait=Traits_nioz, taxonomy=Taxonomy){
  
  wide <- get_trait(taxon    = taxon, 
                    trait    = trait, 
                    taxonomy = taxonomy)
  meta <- metadata(trait)
  long <- w2l_trait(wide, trait_names = meta)
  SUB <- subset(long, subset = !long$modality%in%c("Null", "NULL", "None", "NONE"))
  colnames(SUB)[ncol(SUB)] <- "proportion"
  SUB <- SUB[order(SUB[,1]),]
  zz <- which(colnames(SUB) %in% c("taxon","trait","modality","units","proportion"))
  SUB[,zz]
}

## ====================================================================
## Taxonomic position of a taxon
## ====================================================================

get_taxonomy <- function(taxon, taxonomy=Taxonomy){
  lapply(taxon, FUN = function(tx)
  Taxonomy[unique(which(taxonomy==tx, arr.ind=TRUE)[,1]),])
}

