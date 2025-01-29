## ====================================================================
## ====================================================================
## Summary statistics
## ====================================================================
## ====================================================================


get_summary <- function(data,
                        descriptor, taxon, value, averageOver, 
                        taxonomy = NULL, 
                        subset, 
                        what=c("density", "taxa", "occurrence"), 
                        wide.output = FALSE){

  if (! missing(data)){
    if (! missing(descriptor))
      descriptor  <- eval(substitute(data.frame(descriptor)), 
                          envir = data, enclos = parent.frame())
    if (! missing(taxon))
      taxon       <- eval(substitute(data.frame(taxon)), 
                          envir = data, enclos = parent.frame())
    if (! missing(value))
      value       <- eval(substitute(data.frame(value)), 
                          envir = data, enclos = parent.frame())
    if (! missing(averageOver))
      averageOver <- eval(substitute(data.frame(averageOver)), 
                          envir = data, enclos = parent.frame())
    
    data_name <- substitute(data)
  } else data_name <- NA

   if (pmatch("den", what, nomatch = FALSE)) {
     if(missing(subset))
       density.total <- get_density(
                          descriptor  = descriptor, 
                          averageOver = averageOver, 
                          value       = value)
     else
       density.total <- get_density(
                          descriptor  = descriptor, 
                          averageOver = averageOver, 
                          subset = substitute(subset),
                          value       = value)
     
    density.total <- density.total[, -which(colnames(density.total)=="taxon")]

    } else density.total <- NULL
  
   DD <- NULL
   
   if (any(pmatch(c("tax", "occ"), what, nomatch = FALSE)) &  ! missing(taxon)) {
     if(missing(subset))
       DD <- get_density(
                         descriptor  = descriptor, 
                         averageOver = averageOver, 
                         taxon       = taxon, 
                         value       = value)
     else
       DD <- get_density(
                         descriptor  = descriptor, 
                         averageOver = averageOver, 
                         taxon       = taxon, 
                         subset = substitute(subset),
                         value       = value)
     # number of taxa
     Att   <- attributes(DD)
     Tname <- Att$names_taxon
     Dname <- Att$names_descriptor
     
     nc <- length(Dname)
     
   }
   
   if (pmatch("tax", what, nomatch = FALSE) &  ! missing(taxon)) {
     
     density.numtax <- aggregate(
                DD[, Tname], 
                by  = as.data.frame(DD[, rev(Dname)]),  # reverse the order to have same output as for density
                   FUN = function(x) length(unique(x)))
      
      density.numtax[,1:nc]    <-  density.numtax[,nc:1]  
      colnames(density.numtax) <-  c(Dname, Tname)

      if (nc > 1 & wide.output) 
        density.numtax <- long2wide(
              row    = density.numtax[,1:(nc-1)], 
              column = density.numtax[,nc], 
              value  = density.numtax[,Tname])
    } else 
       density.numtax <- NA

  # for all species: number of descriptors over which it is found. 

   if (pmatch("occ", what, nomatch = FALSE) &  ! missing(taxon)) {
    # Combine multiple columns in descriptor and replicates to average
    Desc   <- data.frame(OneFactor(data.frame(DD[, Dname])))
    species.numdescriptor <- aggregate(
               x   = Desc, 
               by  = as.data.frame(DD[, Tname]),  
               FUN = function(x) length(unique(x)))
    colnames(species.numdescriptor) <-  c(Tname, "occurrence")
    
  } else 
    species.numdescriptor <- NA
  
# for species: mean density over all descriptors, 
   
  list(density    = density.total, 
       taxa       = density.numtax,
       occurrence = species.numdescriptor)
}


