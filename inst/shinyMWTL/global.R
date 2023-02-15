
# ======================
# Extract species data  
# ======================
  
getSpeciesData <- function(DATA     = MWTL, 
                           taxon    = "Amphiura filiformis", 
                           settaxon = "species", 
                           what     = "density", 
                           average  = "date"){

    tx          <- taxon
    averageOver <- NULL

    if (settaxon == "abiotics") {
      Data  <- DATA$abiotics[,c("station", taxon)]
    
    } else if (settaxon == "traits"){
      Data <- getTraitScore(Trait=Traits_nioz, Data = DATA, what=what,
              average=average, tr=tx)
    
    } else if (settaxon == "groups"){
      Data <- getGroupScore(Data = DATA, what=what,
              average=average, tr=tx)
  
    } else if (settaxon == "summary"){
      if (tx == "bioturbation")
        Data <- getDbScore(Data = DATA, average=average)
      else if (tx == "bioirrigation")
        Data <- getIrrScore(Data = DATA, average=average)
      else
        Data <- getSumm(Data = DATA, average=average, what=what, tr=tx)

    } else {
      if (settaxon != "species")   # higher level taxa -not yet implemented
        Data <- merge(DATA$density, Taxonomy)
      else{
        settaxon <- "taxon"  
        Data     <- DATA$density
      }
      colnames(Data)[which(colnames(Data) == what)] <- "val"
    
      if (! is.null(average)) {  # descriptors / variables to average over
        averageOver <- Data[ ,average]
        desc        <- Data$station
      } else {
        desc        <- data.frame(station=Data$station, date=Data$date)
      }    
        Data <- with (Data, 
            long2wide(row=desc, column=taxon, value=val, averageOver=averageOver,
                      subset=(taxon == tx)))
        colnames(Data)[1] <- "station"
    }

    Data    <- merge(DATA$stations, Data, by=1)
    ns      <- ncol(Data)
    present <- which (Data[,ns] > 0)
    Data    <- Data[present, ]
    
    return(Data)
  }  

# ======================
# Trait data
# ======================

getSumm <- function(Data = MWTL, average="date", what="density", tr="density"){
  DATA <- Data$density
  if (is.null(average)) {
    averageOver <- NULL
    descriptor <- data.frame(station=DATA$station, date=DATA$date)
  } else{
    averageOver <- DATA[,average]
    descriptor <- DATA$station
  }
  
  TR <- getSummary(descriptor  = descriptor, 
                   taxon       = DATA$taxon, 
                   value       = DATA[,what], 
                   averageOver = averageOver, 
                   what        = tr)
  TR[[tr]]
}

getDbScore <- function(Traits=Traits_Db, Data = MWTL, average="date"){
  DATA <- Data$density
  if (is.null(average)) {
    averageOver <- NULL
    descriptor <- data.frame(station=DATA$station, date=DATA$date)
  } else{
    averageOver <- DATA[,average]
    descriptor <- DATA$station
  }
  
  TR <- getDbIndex(descriptor  = descriptor, 
                   taxon       = DATA$taxon, 
                   density     = DATA$density, 
                   biomass     = DATA$biomass, 
                   averageOver = averageOver, 
                   taxonomy    = Taxonomy,
                   trait       = Traits)
  TR$descriptor
}

getIrrScore <- function(Traits=Traits_irr, Data = MWTL, average="date"){
  DATA <- Data$density
  if (is.null(average)) {
    averageOver <- NULL
    descriptor <- data.frame(station=DATA$station, date=DATA$date)
  } else{
    averageOver <- DATA[,average]
    descriptor <- DATA$station
  }
  
  TR <- getIrrIndex(descriptor  = descriptor, 
                    taxon       = DATA$taxon, 
                    density     = DATA$density, 
                    biomass     = DATA$biomass, 
                    averageOver = averageOver, 
                    taxonomy    = Taxonomy,
                    trait       = Traits)
  TR
  TR$descriptor
}

getGroupScore <- function(Data = MWTL, what="density",
                          average="date", tr="Small tube dweller"){

  group <- subset(Groups, description==tr)
  stats <- unique(Data$density$station)
  
  DATA <- Data$density
  if (is.null(average)) {
    averageOver <- NULL
    descriptor <- data.frame(descriptor=DATA$station, date=DATA$date)
  } else{
    averageOver <- DATA[,average]
    descriptor <- DATA$station
  }
  
  if (is.matrix(descriptor)|is.data.frame(descriptor))
    cn <- 1:ncol(descriptor)
  else
    cn <- 1
  
  TR <- getDensity(descriptor  = descriptor, 
                   taxon       = DATA$taxon, 
                   value       = DATA[,what], 
                   subset      = taxon %in% group$taxon,
                   averageOver = averageOver)
  by <- TR[,cn]
  if (!is.list(by)) by <- list(by)
  TR <- aggregate(TR$value, by=by, FUN=sum)
  ln <- length(cn)
  colnames(TR)[ln+1] <- tr
  TR
}

getTraitScore <- function(Trait = Traits_nioz, Data = MWTL, what="density",
                          average="date", tr="Substratum depth distribution"){

  metaTrait   <- subset(metadata(Trait), trait==tr)
  subTrait    <- Trait[,c("taxon", metaTrait$colname)]
  
  DATA <- Data$density

  if (is.null(average)) {
    averageOver <- NULL
    descriptor <- data.frame(station=DATA$station, date=DATA$date)
  } else{
    averageOver <- DATA[,average]
    descriptor <- DATA$station
  }

  if (is.matrix(descriptor)|is.data.frame(descriptor))
    cn <- 1:ncol(descriptor)
  else
    cn <- 1
  
  TR <- getTraitDensity(descriptor  = descriptor, 
                        taxon       = DATA$taxon, 
                        value       = DATA[,what], 
                        averageOver = averageOver,
                        trait       = subTrait,
                        trait.class = metaTrait$trait,
                        trait.score = metaTrait$value,
                        taxonomy    = Taxonomy,
                        scalewithvalue = TRUE)
  ln <- length(cn)
  colnames(TR)[ln+1] <- tr
  TR
}

# ======================
# Species data, 
# year as descriptor
# ======================

getSpeciesYear <- function(DATA=MWTL, taxon, settaxon="species", 
                             what="density"){
   DD <- getSpeciesData(DATA, taxon, settaxon, what=what, average=NULL)
   if (is.null(DD$year)) DD$year <- 1900+as.POSIXlt(DD$date, format='%d-%m-%Y')$year
   return(DD)
}

# ======================
# Empty plot - no data
# ======================

plotNone <- function(xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", ...){
  
  plot(x=0, y=0, xlim=xlim, ylim=ylim, type="n", axes=FALSE, 
    xlab=xlab, ylab=ylab, frame.plot=TRUE, ...)
  legend(x="center", legend="NO DATA    ", cex=2)
}

# ======================
# Trait characteristics
# ======================

getTraitVals <- function(Trait=Traits_nioz, taxon="Amphiura filiformis", db="Nioz"){
  Toselect <- c("Maximum size", "Morphology", "Lifespan", "Egg development location",  
  "Larva development location", "Living habit", "Sediment position",  "Feeding mode",              
  "Mobility", "Bioturbation mode", "Ri",  "Mi", "Fti", "BT", "ID", "FT",
    
  "Body mass", "Body length",  "Life span", "Motility", "Feeding type",                
  "Substratum depth distribution", "Ventilation/Pumping",  "Endo-3D structure type", 
  "Epi-3D structure type")

  BioturbMode <-  c("Biodiffusion",  "Downward conveying",  "Upward conveying", 
    "Regeneration", "Biodeposition", "Bioerosion", "Biostabilisation") 
   
  MD   <- metadata(Trait)
  tx   <- taxon
  AFT  <- subset(Trait, taxon==tx)     # select traits for taxon
  if (nrow(AFT) != 0){
  LAFT <- cbind(MD, t(AFT[-1]))        # add metadata
  TTT  <- LAFT[LAFT[,ncol(LAFT)]!=0,]          # remove those that are 0
  Biotur <- subset(TTT, subset=trait %in% BioturbMode)
  TTT  <- subset(TTT, subset=trait %in% Toselect)
  DF <- tapply(TTT$modality, INDEX=TTT$trait, FUN=function(x)x)
  DF <- sapply(DF, FUN=function(x)paste(x, collapse=", "))
  DD <- data.frame(trait=names(DF), value=DF)
  DD <- merge(DD, unique(MD[,c("trait", "units")]),by=1)
  DD$trait[DD$trait== "Substratum depth distribution"] <- "Depth distribution"
  if (! is.null(  Biotur)){
      Biotur <- subset(Biotur, modality != "Null")[,c("trait", "modality", "units")]
      colnames(Biotur)[2] <- "value"
      DD <- rbind(DD, Biotur)
  }
  DD$units[DD$units %in% c("-", " ")] <- ""
  DD$value <- paste(DD$value, DD$units, sep=" ")
  DD$units <- NULL
  G <- subset(Groups, taxon==tx)
  DD <- rbind(c(trait="typology", value=G$description), DD)
  row.names(DD) <- NULL
  } else DD <- "No known traits" 
  DD
}

getTraitSpecs <- function(Traits=Traits_nioz, tr)
    subset(metadata(Traits), trait==tr)[,-1]

getGroupSpecs <- function(Traits=Groups, tr)
    sort(subset(Traits, description==tr)$taxon)

getTaxon <- function(taxon){
  tx <- taxon
  DD <- as.data.frame(t(subset(Taxonomy, taxon==tx)[-1]))
  names(DD) <- tx
  DD
}

GetGroupTaxons <- function(group){

  Nst <- attributes(Groups)$description
  type <- Nst[Nst$description==group, 1]
  GG <- subset(Groups, typology == type)
  return(paste(GG[,1], collapse=",  ", sep=""))
}

