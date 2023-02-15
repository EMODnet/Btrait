\name{Northsea Benthos}
\docType{data}
\alias{NSBS}
\title{
  The ICES northsea benthic survey (NSBS) from 1986, including abiotic conditions.
}

\description{
  The NSBS Northsea macrobenthos data (1986).
  
  Samples are within (-3 and 9) dgE and (51.00 60.75) dgN
  
  The data contain:
  
  \itemize{
    \item Macrofauna species density and biomass (\code{NSBS$density})
  
    \item Abiotic conditions (\code{NSBS$abiotics}), and station positions (\code{NSBS$stations})
  
    \item \code{NSBS$contours}: contourlines for mapping.
  }
}

\usage{
data(NSBS)
}

\format{
==================

**\code{NSNS$density}**: 
This is the main Northsea NSBS benthos data set, containing species 
information for 234 stations sampled in 1985-1986.
 
The data, in long format, are in a \code{data.frame} with the following columns: 
  \itemize{
   \item \code{station}, the NSBS station name (details in NSBS$stations).    \item \code{date}, the sampling date. 
   Note, this is a string; it can be converted to POSIXct by:
    "as.POSIXct(NSBS$density$date, format='\%d-\%m-\%Y')".
    The year can be extracted as 
    \code{1900+as.POSIXlt(NSBS$density$date, format='\%d-\%m-\%Y')$year}.
   \item \code{year}, the year of sampling
   \item \code{taxon}, the taxon name to be used (usually species), and checked against the worms database (details in dataset \code{Taxonomy}).
   \item \code{density}, the number of individuals per m2.
   \item \code{biomass}, the total biomass per m2, in AFDW/m2 (ash-free dry weight).
  }

==================

**\code{NSBS$abiotics}**: 
the abiotic conditions of sampling stations.

NSBS$abiotics is a data.frame with the following columns:   
\itemize{
  \item{station}, the NSBS station name
  \item{depth}, water depth [m]
  \item{mud}, mud content of sediment, fraction [-]
  \item{porosity}, volumetric water content [-]
  \item{permeability}, sediment permeability [m2]
  \item{max.currentspeed}, maximum current speed [m/s]
 }

==================

**\code{NSBS$stations}**: 

The positions of the different stations, in WGS84 format
 \itemize{
  \item{station}, the NSBS station name
  \item{x}, degrees longitude
  \item{y}, degrees latitude
 }
 
==================

 **\code{NSBS$contour}**: 
 
 The data for mapping the contours. The contourlines (x-, y) were derived from 
 GEBCO high-resolution bathymetry, by creating contourlines.
 
 The data set contains:
   \itemize{
      \item{station}, the NSBS station name
      \item{x}: longitude, in [dgE]
      \item{y}: latitude, in [dgN]
      \item{z}: the corresponding depths, in [m]
    }
 
}
\note{
The dataset **\code{Taxonomy}**: 
contains taxonomic information of the original and adjusted taxon in NSBS$density, 
as derived from the World Register of Marine Species (WORMS), using R-package \code{worms}.
}

\details{
NSBS contains the macrofauna data from the 1986 North Sea Benthos Survey, an activity of the Benthos Ecology Working Group of ICES. 

Benthic samples were taken in a standardised way, on a regular grid covering the whole of the North Sea, and analysed by scientists from 10 laboratories. Extensive work was done to standardise taxonomy and identifications across the different laboratories.

Sediment was sampled with a Reineck Boxcorer (0,078 m2). Macrofauna sieved on a 1 mm mesh. 
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{

##-----------------------------------------------------
## Show contents of the data set
##-----------------------------------------------------

metadata(NSBS$abiotics)
metadata(NSBS$density)

##-----------------------------------------------------
## SPECIES data
##-----------------------------------------------------

head(NSBS$density)

# The number of species per station (over all years)
Nspecies <- tapply(X     = NSBS$density$taxon, 
                   INDEX = NSBS$density$station, 
                   FUN   = function(x)length(unique(x)))
summary(Nspecies)

# The number of times a species has been found
Nocc     <- tapply(X     = NSBS$density$station, 
                   INDEX = NSBS$density$taxon, 
                   FUN   = length)
head(sort(Nocc, decreasing = TRUE))     #most often encountered taxa

# total density per station
densstat <- tapply(X     = NSBS$density$density, 
                   INDEX = list(NSBS$density$station), 
                   FUN   = sum)
hist(densstat, n=30)

##-----------------------------------------------------
## ABIOTICS
##-----------------------------------------------------

summary(NSBS$abiotics)

NSBSab <- merge(NSBS$stations, NSBS$abiotics)
with(NSBSab, mapBtrait(x, y, colvar=mud, contours=NSBS$contours, 
                       pch=16, main="mud fraction"))

metadata(NSBS$abiotics)

##-----------------------------------------------------
## COMBINATIONS
##-----------------------------------------------------

NSsp_abi <- merge(NSBS$density, NSBS$abiotics)
ECH      <- subset(NSsp_abi, subset=taxon=="Echinocardium cordatum")
with(ECH, plot(mud, density, pch=16, main="E. cordatum", 
               xlab="mud fraction", ylab="density, ind/m2"))

# add station coordinates
ECH <- merge(ECH, NSBS$stations)

##-----------------------------------------------------
## From long format to wide format (stations x species)
##-----------------------------------------------------

NSwide <- with (NSBS$density, 
         l2wDensity(descriptor = station,    # long2wide for density
                    taxon      = taxon, 
                    value      = density))

PP <- princomp(t(NSwide[,-1]))
\dontrun{
 biplot(PP)
}

##-----------------------------------------------------
## Community mean weighted score for traits.
##-----------------------------------------------------

# Traits estimated for absences, by including taxonomy 

Trait.lab <- metadata(Traits_nioz)

trait.cwm <- getTraitDensity(wide           = NSwide, 
                             trait          = Traits_nioz, 
                             taxonomy       = NSBS$taxonomy,
                             trait.class    = Trait.lab$trait, 
                             trait.score    = Trait.lab$score, 
                             scalewithvalue = TRUE)

head(trait.cwm, n=c(3,4))  

Stations.traits <- merge(NSBS$stations, trait.cwm, 
                         by.x="station", by.y="descriptor")

par(mfrow=c(2,2))

with(Stations.traits, mapBtrait(x, y, colvar=Biodeposition,
                                main="Biodeposition"))
with(Stations.traits, mapBtrait(x, y, colvar=Biodiffusion,
                                main="Biodiffusion"))
with(Stations.traits, mapBtrait(x, y, colvar=Biostabilisation,
                                main="Biostabilisation"))
with(Stations.traits, mapBtrait(x, y, colvar=Burrow.width,
                                main="Burrow width"))

##-----------------------------------------------------
## Community mean weighted score for typological groups.
##-----------------------------------------------------

# Groups is in crisp format -> convert to fuzzy 

Groups.fuz <- crisp2fuzzy(Groups[,c("taxon", "typology")])

head (Groups, n=2)
head (Groups.fuz, n=c(2,5))

group.cwm <- getTraitDensity(wide           = NSwide, 
                             trait          = Groups.fuz, 
                             scalewithvalue = TRUE)

head(group.cwm, n=c(3,4))  

summary(group.cwm)

##-----------------------------------------------------
## Show the depth contours
##-----------------------------------------------------

mapBtrait(contours=NSBS$contours, draw.levels=TRUE, key.levels=TRUE)

# Use a different color scheme
collev <- function(n) c("black", 
      ramp.col(col=c("darkgreen", "darkblue"), n=n-1))
mapBtrait(contours=NSBS$contours, draw.levels=TRUE, col.levels=collev,
  key.levels=TRUE)
}

\references{

The taxonomic information was created using the \code{worms} package:

Jan Holstein (2018). worms: Retrieving Aphia Information from World Register of Marine Species. R package version 0.2.2. https://CRAN.R-project.org/package=worms

Heip, C.H.R.; Basford, D.; Craeymeersch, J.A.; Dewarumez, J.-M.; Dorjes, J.; de Wilde, P.; Duineveld, G.; Eleftheriou, A.; Herman, P.M.J.; Kingston, K.; Niermann, U.; Kunitzer, A.; Rachor, E.; Rumohr, H.; Soetaert, K.; Soltwedel, T. (1992). Trends in biomass, density and diversity of North Sea macrofauna. ICES J. Mar. Sci./J. Cons. int. Explor. Mer 49: 13-22
}

\seealso{
  \link{mapBtrait} for plotting.

  \link{Traits_nioz} for the trait datasets.
  
  \link{getDensity} for functions operating on these data.

  \link{getDbIndex} for extracting bioturbation and bioirrigation indices.

  \link{long2wide} for functions changing the appearance on these data.
}

\keyword{datasets}

