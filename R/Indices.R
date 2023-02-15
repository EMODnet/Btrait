## ====================================================================
## ====================================================================
## Bioturbation and bioirrigation potentials
## ====================================================================
## ====================================================================

getDbIndex <- function(descriptor, 
                       taxon, 
                       density, 
                       biomass,
                       averageOver=NULL, 
                       subset,
                       trait=Traits_Db,  # species x trait data, Mi and Ri    
                       taxonomy=NULL,    # if !NULL, trait will be expanded at higher levels
                       full.output=FALSE,
                       verbose=FALSE)  
  getPotential(descriptor = descriptor, 
               taxon      = taxon, 
               density    = density, 
               biomass    = biomass, 
               averageOver= averageOver, 
               subset     = subset,
               trait      = trait, 
               taxonomy   = taxonomy, 
               full.output= full.output,
               verbose    = verbose, 
               type       = "BPc")  

## ====================================================================

getIrrIndex <- function(descriptor, 
                         taxon, 
                         density, 
                         biomass,
                         averageOver=NULL, 
                         subset,
                         trait=Traits_irr,  # species x trait data, FT, BT, ID    
                         taxonomy=NULL,     # if !NULL, trait will be expanded at higher levels
                         full.output=FALSE,
                         verbose=FALSE)  
  getPotential(descriptor = descriptor, 
               taxon      = taxon, 
               density    = density, 
               biomass    = biomass, 
               averageOver= averageOver, 
               subset     = subset,
               trait      = trait, 
               taxonomy   = taxonomy, 
               full.output= full.output,
               verbose    = verbose,
               type       = "IPc")  

## ====================================================================
## General function
## ====================================================================
  getPotential <- function(descriptor, 
                         taxon, 
                         density, 
                         biomass,
                         averageOver=NULL, 
                         subset,
                         trait,         # species x trait data    
                         taxonomy=NULL, # if !NULL, trait will be expanded at higher levels
                         full.output=FALSE,
                         verbose=FALSE, 
                         type="BPc")  

{ 
  if (missing(descriptor)) descriptor <- NA
  
  isnaDesc <- all(is.na(descriptor))
  cnDesc   <- getname(descriptor)
  cnTaxon  <- getname(taxon)
  
  if (isnaDesc) descriptor <- 1:length(taxon)  # to prevent summation over taxon
  BioDens <- getDensity(descriptor  = descriptor,
                        subset      = subset,
                        averageOver = averageOver, 
                        taxon       = taxon,
                        value       = data.frame(density=density, biomass=biomass),
                        taxonomy    = taxonomy,
                        verbose     = verbose) 
  
  names(BioDens)[(1:3)+length(cnDesc)] <- c("taxon", "density", "biomass") 
  
# estimate mean weight of taxa (=biomass/density)
  BioDens$Weight <- BioDens$biomass/BioDens$density  #biomass/ind

# estimate DB/irr traits, also for species in database if (taxonomy not NULL)
  taxname <- names(trait)[1]

  if (type == "BPc") 
    traitnames <- c(taxname, "Ri", "Mi")
  else if (type == "IPc")
    traitnames <- c(taxname, "BT", "FT", "ID")

  Pot.traits <- getTrait(
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
