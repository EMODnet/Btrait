\name{Dutch Northsea data}
\docType{data}
\alias{MWTL}
\title{
  The benthic fauna data from the Dutch part of the Northsea, including abiotic conditions.
}

\description{
  The MWTL Northsea macrobenthos data (1995 - 2018)
  
  The dataset contains:
  \itemize{

  \item Macrofauna species density and biomass (MWTL$density)
  
  \item abiotic conditions (MWTL$abiotics), station types (MWTL$types), 
  sediment composition (MWTL$sediment) and station positions (MWTL$stations)
  
  \item MWTL$contour: depth contourlines for mapping.
 }
}

\usage{
data(MWTL)
}

\format{
==================

**\code{MWTL$density}**: 
This is the main Northsea MWTL benthos data set, containing species 
information for 103 stations sampled on a yearly basis from 1995 till 2010, 
after which sampling was less frequent: in 2012, 2015 and 2018.
 
The data, in long format, are in a \code{data.frame} with the following columns: 
  \itemize{
   \item \code{station}, the MWTL station name (details in MWTL$stations).
   \item \code{date}, the sampling date. 
   Note, this is a string; it can be converted to POSIXct by:
    \code{as.POSIXct(MWTL$density$date, format='\%d-\%m-\%Y')}. 
    The year can be extracted as 
    \code{1900+as.POSIXlt(MWTL$density$date, format='\%d-\%m-\%Y')$year}.
   \item \code{year}, the year of sampling
   \item \code{taxon}, the taxon name to be used (usually species); this has been derived from the original taxon in the MWTL data as follows: 
   The original taxon e.g. species is kept, if a minimum of 90\% of individual 
   organisms of this taxon are at the species level, genus otherwise, if not family etc. The taxon name was checked against the worms database (details in dataset \code{Taxonomy}).
   \item \code{density}, the number of individuals per m2.
   \item \code{biomass}, the total biomass per m2, in AFDW/m2 (ash-free dry weight)
   \item \code{taxon.original}, the original taxon name in the data set (see details).
  }

==================

A \code{data.frame} with the following columns: 
\code{taxon}, \code{genus}, \code{family}, \code{order}, \code{class}, 
\code{phylum}, \code{AphiaID}.

=========================================

**\code{MWTL$abiotics} and \code{MWTL$types}**: 
the abiotic conditions of sampling stations,
numerical values (MWTL$abiotics) or typologies (MWTL$types), averaged over all years.

MWTL$abiotics is a data.frame with the following columns:   CHECK IT!
\itemize{
  \item{station}, the MWTL station name.
  \item{current}, the mean current speed, [m/s]
  \item{wave},   wave energy              [Pa]
  \item{disturb},           ???             [???]
  \item{sal},      salinity                  [-]
  \item{pelag.PP}, pelagic Primary Production (spring) [?? mg C m-2 d-1???]
  \item{mud},                        mud             [-]
  \item{gravel},                     gravel          [-]
  \item{sand},                       sand            [-]
  \item{PP},        pelagic PP (sping) [??  mg C m-2 d-1??]
  \item{POM}, particulate organic matter               [\%]
  \item{POC}, particulate organic carbon               [\%]
  \item{Aspect},                     Aspect          [???]
  \item{Curvature},                  Curvature       [???]
  \item{Rugosity},                  Rugosity2       [???]
  \item{Northing},                   Northing        [???]
  \item{Slope},                      Slope           [???]
 }

MWTL$types categorizes the stations into a number of types:
\itemize{
   \item{station}, the MWTL station name.
   \item{depth} (m)
  \itemize{
    \item{Very shallow}: < 10, 
    \item{Shallow}: [10; 20[, 
    \item{Intermediate}: [20; 30[, 
    \item{Deep}: [30; 40[, 
    \item{Very deep}: >= 40
  }
  \item{current} (m/s): 
  \itemize{
   \item{Very low}: < 0.15, 
   \item{Low}: [0.15; 0.20[, 
   \item{Intermediate}: [0.20; 0.25[, 
   \item{High}: [0.25; 0.30[, 
   \item{Very high}: >= 0.30
  } 
  The numbers are Monthly median values, in meters per second averaged from 1996 to 2008
  \item{wave} energy (Pa):  
   \itemize{
    \item{Very low}: < 0.5, 
    \item{Low}: [0.5; 1.0[,
    \item{Intermediate}: [1.0; 1.5[, 
    \item{High}: [1.5; 2.0[, 
    \item{Very high}: >= 2.0
  }; the data are monthly median values in pascals averaged from 1996 to 2008.
  \item{stratification} 
   \itemize{
    \item{PM} : Permanently mixed 
    \item{FI} : Freshwater influence 
    \item{IS} : Intermittently stratified 
    \item{SS} : Seasonally stratified 
    \item{TR} : Transitional
   }
  \item{sediment} 
   \itemize{
   \item{Muddy} includes Mud, Sandy mud, Sandy and slightly gravely mud, Muddy sand; 
   \item{Sandy} means only Sand; 
   \item{Coarse} includes: Gravel and muddy sand, Slightly gravely sand, Gravely sand, Sandy gravel, Gravel and stone; 
   \item{Mixed} includes Gravely and slightly muddy sand, 
   }  
   where clay <8 um, silt: 8-63,  very fine sand: 62-125, fine sand: 125-250, medium sand: 250-500, coarse sand: 500-1000, very coarse sand: 1000-2000, gravel: >2000 um (micrometer).
  \item{BPI_xx} Benthic Terrain Classification parameters, bathymetric position index
  \item{dynamics} 
   \itemize{
    \item{L}: low
    \item{H}: high
   } 
  \item{area} 
   \itemize{
     \item{DOG}: Doggersbank 
     \item{OYS}: Oystergrounds 
     \item{OFF}: Offshore
     \item{COA}: Coastal zone
   }
}

==================

**\code{MWTL$stations}**: 

The positions of the different stations, in WGS84 format
 \itemize{
  \item{station}, the MWTL station name.  
  \item{x}, degrees longitude
  \item{y}, degrees latitude
 }
 
==================

**\code{MWTL$sediment}**: 
records the sediment characteristics (median grain size and silt content)
for the different stations and years.
The data for "BREEVTN11" "BREEVTN16" "WADDKT05" are missing.

This is a \code{data.frame} with:
 \itemize{
  \item{station}, the MWTL station name. 
  \item{year}, the sampling year.  
  \item{D50}, the median grain size of the sediment, in micrometer
  \item{silt}, the silt and clay content (< 63 micrometer), in percentage
 }

De sediment grainsize was determined by laserdiffraction using a Malvern Mastersizer. 
Values denote weight percentages dryweight of the total sediment sample, where big shells and large animals were removed.

==================

 **\code{MWTL$contours}**: 
 
 The data for plotting the depth contours in the area. 
 The contourlines (x-, y) were derived from GEBCO high-resolution bathymetry, 
 by using the \code{contourLines} R-function.
 
 The data set contains:
   \itemize{
      \item{x}: longitude, in [dgE]
      \item{y}: latitude, in [dgN]
      \item{z}: the corresponding depths, [m]
    }
 
}
\note{
The dataset **\code{Taxonomy}**: 
contains taxonomic information of the original and adjusted taxon in MWTL$density, 
as derived from the World Register of Marine Species (WORMS), using R-package \code{worms}.
}
\details{
 The macrobenthos data of the Northsea (MWTL) are commissioned by the 
 "Ministerie van Infrastructuur en Milieu, Rijkswaterstaat Centrale 
 Informatievoorziening (RWS, CIV)". 

MWTL stands for Monitoring Waterstaatkundige Toestand des Lands (dutch).
 
Sediment was sampled with a Reineck Boxcorer (0,078 m2). Macrofauna sieved on a 1 mm mesh. 
All animals determined, except when too much residue or organisms, 
in which case samples were subsampled, so that for molluscs and crustaceans 
at least 100 individuals, and for polychaetes at least 150 individuals were determined. 

Biomass was determined as ash-free dryweight: individuals were dried for >48 hours 
at 65 dgC, and then cooled in an exsiccator for at least 30 minutes and weighed 
(precision 0,01 mg), determining their dryweight. Then they were ashed in an 
oven at 530 dgC (2,5, 4 or 8 hour, depending on the size of organisms). 
Following ashing, they were weighed, after cooling for at least 45 minutes in 
an exsiccator. Bivalvia and Gastropoda =7 mm were ashed without shell, but when 
smaller than 7 mm the shell was not removed. 

AFDW = (dryweight + weight cup) ? (ash weight + weight cup).
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{

##-----------------------------------------------------
## Show contents of the data set
##-----------------------------------------------------

metadata(MWTL$sediment)
metadata(MWTL$abiotics)
metadata(MWTL$types)
metadata(MWTL$density)

##-----------------------------------------------------
## SPECIES data
##-----------------------------------------------------

head(MWTL$density)

# The number of species per station (over all years)

Nspecies <- tapply(X     = MWTL$density$taxon, 
                   INDEX = MWTL$density$station, 
                   FUN   = function(x)length(unique(x)))

summary(Nspecies)

# Per year
Nspyear <- tapply(X     = MWTL$density$taxon, 
                  INDEX = list(MWTL$density$station, MWTL$density$year), 
                  FUN   = function(x)length(unique(x)))

colMeans(Nspyear, na.rm=TRUE)

# The number of times a species has been found
Nocc <- tapply(X     = MWTL$density$station, 
               INDEX = MWTL$density$taxon, 
               FUN   = length)
head(sort(Nocc, decreasing = TRUE))     #most often encountered taxa

# total density
densyear <- tapply(X     = MWTL$density$density, 
                   INDEX = list(MWTL$density$station, MWTL$density$year), 
                   FUN   = sum)
boxplot(densyear, log="y")

##-----------------------------------------------------
## ABIOTICS data
##-----------------------------------------------------

summary(MWTL$abiotics)

MWTLabiot <- merge(MWTL$stations, MWTL$abiotics)
with(MWTLabiot, mapBtrait(x, y, colvar=sand, 
                          pch=16, main="sand fraction"))

# mud, plotted on large Northsea map
with(MWTLabiot, mapBtrait(x, y, colvar=mud, contours=NSBS$contours, 
                          pch=16, main="mud fraction"))

# show the different abiotic data sets
metadata(MWTL$abiotics)

##-----------------------------------------------------
## COMBINATIONS
##-----------------------------------------------------

NSsp_abi <- merge(MWTL$density, MWTL$sediment)
ECH      <- subset(NSsp_abi, subset=taxon=="Echinocardium")

with(ECH, points2D(D50, density, log="xc", colvar=density, pch=16))

# add station coordinates
ECH <- merge(ECH, MWTL$stations)

##-----------------------------------------------------
## From long format to wide format (stations x species)
##-----------------------------------------------------

NSwide <- with (MWTL$density, 
     l2wDensity(descriptor  = station, 
                taxon       = taxon, 
                value       = density, 
                averageOver = year))

PP <- princomp(t(NSwide[,-1]))
\dontrun{
 biplot(PP)
}

##-----------------------------------------------------
## Community weighted mean score.
##-----------------------------------------------------

# Traits estimated for absences, by including taxonomy 

NStrait.lab <- metadata(Traits_nioz)
trait.cwm <- getTraitDensity (wide           = NSwide, 
                              trait          = Traits_nioz, 
                              taxonomy       = Taxonomy,
                              trait.class    = NStrait.lab$trait, 
                              trait.score    = NStrait.lab$score, 
                              scalewithvalue = TRUE)

head(trait.cwm, n=c(3,4))  

Stations.traits <- merge(MWTL$stations,  trait.cwm, 
                         by.x="station", by.y="descriptor")

##-----------------------------------------------------
## Maps
##-----------------------------------------------------

par(mfrow=c(2,2))

with(Stations.traits, mapBtrait(x, y, colvar=Biodeposition,
                                main="Biodeposition", pch=16))
with(Stations.traits, mapBtrait(x, y, colvar=Biodiffusion,
                                main="Biodiffusion", pch=16))
with(Stations.traits, mapBtrait(x, y, colvar=Biostabilisation,
                                main="Biostabilisation", pch=16))
with(Stations.traits, mapBtrait(x, y, colvar=Burrow.width,
                                main="Burrow width", pch=16))

##-----------------------------------------------------
## Show the depth contours
##-----------------------------------------------------

mapBtrait(contours=MWTL$contours, draw.levels=TRUE, key.levels=TRUE)

# Use a different color scheme
collev <- function(n) c("black", 
  ramp.col(col=c("darkgreen", "darkblue"), n=n-1))
mapBtrait(contours=MWTL$contours, draw.levels=TRUE, col.levels=collev,
  key.levels=TRUE)
}

\references{

The taxonomic information was created using the \code{worms} package:

Jan Holstein (2018). worms: Retrieving Aphia Information from World Register of Marine Species. R package version 0.2.2. https://CRAN.R-project.org/package=worms

L. Leewis, E.C. Verduin, R. Stolk ; Eurofins AquaSense Macrozoobenthosonderzoek in de Rijkswateren met boxcorer, jaarrapportage MWTL 2015 : waterlichaam: Noordzee
Publicatiedatum: 31-03-201775 p. Projectnummer Eurofins AquaSense: J00002105. 
Revisie 2, In opdracht van Ministerie van Infrastructuur en Milieu,  Rijkswaterstaat Centrale Informatievoorziening (RWS, CIV)
}

\seealso{
  \link{mapBtrait} for plotting.

  \link{Traits_nioz} for the trait datasets.
  
  \link{getDensity} for functions operating on these data.

  \link{long2wide} for functions changing the appearance on these data.

  \link{getDbIndex} for extracting bioturbation and bioirrigation indices.
}

\keyword{datasets}

