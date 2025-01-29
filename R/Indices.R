## ====================================================================
## ====================================================================
## Bioturbation and bioirrigation potentials
## ====================================================================
## ====================================================================

get_Db_index <- function(data,
                       descriptor, 
                       taxon, 
                       density, 
                       biomass,
                       averageOver, 
                       subset,
                       trait=Traits_Db,  # species x trait data, Mi and Ri    
                       taxonomy=NULL,    # if !NULL, trait will be expanded at higher levels
                       full.output=FALSE,
                       verbose=FALSE)  {
  if (! missing(data)){
    if (! missing(descriptor))
      descriptor  <- eval(substitute(data.frame(descriptor)), 
                          envir = data, enclos = parent.frame())
    if (! missing(taxon))
      taxon       <- eval(substitute(data.frame(taxon)), 
                          envir = data, enclos = parent.frame())
    if (! missing(density))
      density       <- eval(substitute(data.frame(density)), 
                          envir = data, enclos = parent.frame())
    if (! missing(biomass))
      biomass       <- eval(substitute(data.frame(biomass)), 
                            envir = data, enclos = parent.frame())
    if (! missing(averageOver))
      averageOver <- eval(substitute(data.frame(averageOver)), 
                          envir = data, enclos = parent.frame())
    data_name <- substitute(data)
  } else data_name <- NA
  
  if (missing(averageOver))
    averageOver <- NULL

  if (missing(descriptor)) descriptor <- NA
  
  isnaDesc <- all(is.na(descriptor))
  
  if (isnaDesc) descriptor <- 1:length(taxon)  # to prevent summation over taxon
  
  if (missing(subset))
    BioDens <- get_density(descriptor  = descriptor,
                           averageOver = averageOver, 
                           taxon       = taxon,
                           value       = data.frame(density = density, 
                                                    biomass = biomass),
                           taxonomy    = taxonomy,
                           verbose     = verbose) 
  else 
    BioDens <- get_density(descriptor  = descriptor,
                           averageOver = averageOver, 
                           taxon       = taxon,
                           value       = data.frame(density = density, 
                                                    biomass = biomass),
                           subset     = substitute(subset),
                           taxonomy    = taxonomy,
                           verbose     = verbose) 
  Db <- getPotential(BioDens, trait = trait, taxonomy = taxonomy , 
                      isnaDesc = isnaDesc, 
                      full.output = full.output, type       = "BPc")  
  
  attributes(Db) <- c(attributes(Db), attributes(BioDens)[-(1:3)])
  attributes(Db)$dataset <- data_name
  Db
}
## ====================================================================

get_irr_index <- function(
                         data,
                         descriptor, 
                         taxon, 
                         density, 
                         biomass,
                         averageOver, 
                         subset,
                         trait=Traits_irr,  # species x trait data, FT, BT, ID    
                         taxonomy=NULL,     # if !NULL, trait will be expanded at higher levels
                         full.output=FALSE,
                         verbose=FALSE)  {
  if (! missing(data)){
    if (! missing(descriptor))
      descriptor  <- eval(substitute(data.frame(descriptor)), 
                          envir = data, enclos = parent.frame())
    if (! missing(taxon))
      taxon       <- eval(substitute(data.frame(taxon)), 
                          envir = data, enclos = parent.frame())
    if (! missing(density))
      density       <- eval(substitute(data.frame(density)), 
                            envir = data, enclos = parent.frame())
    if (! missing(biomass))
      biomass       <- eval(substitute(data.frame(biomass)), 
                            envir = data, enclos = parent.frame())
    if (! missing(averageOver))
      averageOver <- eval(substitute(data.frame(averageOver)), 
                          envir = data, enclos = parent.frame())
    data_name <- substitute(data)
  } else data_name <- NA
  
  if (missing(averageOver))
    averageOver <- NULL
  
  if (missing(descriptor)) descriptor <- NA
  
  isnaDesc <- all(is.na(descriptor))

  if (isnaDesc) descriptor <- 1:length(taxon)  # to prevent summation over taxon
  
  if (missing(subset))
    BioDens <- get_density(descriptor  = descriptor,
                           averageOver = averageOver, 
                           taxon       = taxon,
                           value       = data.frame(density = density, 
                                                  biomass = biomass),
                           taxonomy    = taxonomy,
                           verbose     = verbose) 
   else 
     BioDens <- get_density(descriptor  = descriptor,
                            averageOver = averageOver, 
                            taxon       = taxon,
                            value       = data.frame(density = density, 
                                                     biomass = biomass),
                            subset     = substitute(subset),
                            taxonomy    = taxonomy,
                            verbose     = verbose) 
   Irr <- getPotential(BioDens, trait = trait, taxonomy = taxonomy , 
                       isnaDesc = isnaDesc, 
                       full.output = full.output, type       = "IPc")  

  attributes(Irr) <- c(attributes(Irr), attributes(BioDens)[-(1:3)])
  attributes(Irr)$dataset <- data_name
  Irr
  
}
## ====================================================================
## General function
## ====================================================================
  getPotential <- function(BioDens,  
                           trait=Traits_Db,  # species x trait data, FT, BT, ID   
                           taxonomy = NULL,
                           isnaDesc = NULL,
                           full.output = FALSE,
                           type="BPc")  

{ 
    
  Att <- attributes(BioDens)
  cnDesc <- Att$names_descriptor
  cnTaxon <- Att$names_taxon
  names(BioDens)[(1:3)+length(cnDesc)] <- c("taxon", "density", "biomass") 
  
# estimate mean weight of taxa (=biomass/density)
  BioDens$Weight <- BioDens$biomass/BioDens$density  #biomass/ind

# estimate DB/irr traits, also for species in database if (taxonomy not NULL)
  taxname <- names(trait)[1]

  if (type == "BPc") 
    traitnames <- c(taxname, "Ri", "Mi")
  else if (type == "IPc")
    traitnames <- c(taxname, "BT", "FT", "ID")

  Pot.traits <- get_trait(
                     taxon    = unique(BioDens$taxon), 
                     trait    = trait[,traitnames], 
                     taxonomy = taxonomy)

# taxon for which trait could not be derived
  notrait <- attributes(Pot.traits)$notrait  

# merge density/weight with traits
  BioDens  <- merge(BioDens, Pot.traits, by="taxon")

# species BPc/IPc
  if (type == "BPc")
    BioDens$BPc <- with (BioDens, sqrt(Weight)*density*Mi*Ri)
  
  else if (type == "IPc")
    BioDens$IPc <- with(BioDens, Weight^(0.75)*density*BT*FT*ID)

# results: SUM over descriptor and MEAN over taxa
  PC <- list()
  if (!isnaDesc) {
    by <-  BioDens[,cnDesc]  # will be a list (data.frame) if multiple columns
    if (! is.list(by)) by <-  list(by)
    
    PC$descriptor <- aggregate(
                          x   = BioDens[,type],    # data on which to calculate
                          by  = by,                # for each descriptor
                          FUN = sum, na.rm=TRUE)   # take the sum
    colnames(PC$descriptor) <- c(cnDesc, type)
  }  else PC$descriptor <- NA                    

  PC$taxon <- aggregate(x  = BioDens[,type],       # data on which to calculate
                       by  = list(BioDens$taxon),  # for each station
                       FUN = mean            )     # take the sum

  
  colnames(PC$taxon) <- c(cnTaxon, type)
  PC$all <- NA
  if (full.output) {
    PC$all <- BioDens
    if (isnaDesc) 
       PC$all <- PC$all[order(PC$all$descriptor),]
  }
  attributes(PC)$notrait <- notrait
  PC
}  
