# Standardize trait modalities to sum = 1 - check if all taxa have at least one  modality

standardize_trait <- function(trait, taxon_column = 1, 
                              trait_class, trait_modality,
                              meta = data.frame(trait_class    = trait_class, 
                                                trait_modality = trait_modality))
{
  
  # length of trait_class should be = ncol(trait) - length(taxon_column)
  if (!is.numeric(taxon_column))
    taxon_column  <- which(colnames(trait) == taxon_column)
  
  spname <- colnames(trait)[taxon_column]
  
  if (taxon_column == 0)
    cn     <- colnames(trait)
  else
    cn     <- colnames(trait)[ -(1:taxon_column)]
  ii     <- which (!cn %in% meta$trait_modality)
  
  if (length(ii))
    stop (" cannot standardize: not all modalities (column names of trait) found in meta")
  
  # to long format (absences so that we can check for at least one modality per trait)
  NT   <- w2l_trait(trait, taxon_column = taxon_column, absences = TRUE)  
  nr   <- nrow(NT)
  
  NT <- merge(NT, meta, by.x = 2, by.y = "trait_modality")
  colnames(NT) <- c("modality", "taxon", "value", "trait_class")
  
  if (sum(NT$value) != sum(trait[, -taxon_column])) 
    stop (" cannot standardize: error when going to long format")
  
  NTsum <- aggregate(x   = list(sum = NT$value), 
                     by  = list(taxon = NT$taxon, trait_class = NT$trait_class), 
                     FUN = sum)
  
  ii <- which(NTsum$sum == 0) 
  if (length(ii)) stop ("some taxa have no traits, i.e. for species ", 
                        paste(head(NTsum[ii,1]), collapse = ", "), 
                        "  and modality ", paste(head(NTsum[ii,2]), collapse = ", "))

  NT <- merge(NT, NTsum, by = c("taxon", "trait_class"))
  NT$relative <- NT$value/NT$sum
  
  # to wide format
  New_traits <- l2w_trait(taxon      = NT$taxon, 
                          descriptor = NT$modality, 
                          value      = NT$relative)
  
  colnames(New_traits)[1] <- spname
  New_traits
}

standardize_taxonomy <- function(taxonomy, verbose = FALSE){
  
  Taxonomy  <- unique(taxonomy)
  ii        <- which(duplicated(Taxonomy[,1]))
  if (length(ii)){
    if (verbose) {
      warning("First column in taxonomy has ", length(ii),
              "duplicated values,\n the first ones being, ",
              paste(head(Taxonomy[ii,1]), collapse = ","))
      warning("duplicates will be removed")
    }
    Taxonomy  <- Taxonomy[-ii, ]
  }

  for (i in 2:ncol(Taxonomy)){
    
    # Fill in blancs (where determination was not up to this taxon)
    
    i2 <- which(Taxonomy[,i] %in% c(""," ", "  ", "   "))
    
    if (length(i2))
      Taxonomy[i2, i] <- Taxonomy[i2, i-1]
  }
  Taxonomy
}
