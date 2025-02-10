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
  
  \item MWTL$contours: depth contourlines for mapping.

  \item MWTL$fishing: species trait values that can be used to estimate fishing parameters.

  \item MWTL$sar: station-specific fishing intensities.
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
   organisms of this taxon are at the species level, genus otherwise, if not family etc. The taxon name was checked against the WoRMS database (details in dataset \code{Taxonomy}).
   \item \code{density}, the number of individuals per m2.
   \item \code{biomass}, the total biomass per m2, in AFDW/m2 (ash-free dry weight)
   \item \code{taxon.original}, the original taxon name in the data set (see details).
  }

=========================================

**\code{MWTL$abiotics}**: 
the abiotic conditions of sampling stations.

MWTL$abiotics is a data.frame with the following columns:   
\itemize{
  \item{station}, the NSBS station name
  \item{depth}, water depth, [m]
  \item{D50}, Median grain size, [micrometer]
  \item{mud}, mud content of sediment (<63 um), fraction, [-]
  \item{sand},         sand fraction (64 -2000 um),          [-]
  \item{gravel},          gravel fraction (>2000 um),          [-]
  \item{salinity},                            salinity           
  \item{porosity}, volumetric water content, [-]
  \item{permeability}, sediment permeability, [m2]
  \item{POC},   particulate organic C in sediment,          [percent]
  \item{TN},                 total N in sediment,          [percent]
  \item{surfaceCarbon},   particulate organic C in upper cm,          [percent]
  \item{surfaceNitrogen},                 total N in upper cm,          [percent]
  \item{orbitalVelMean},               mean orbital velocity,        [m/s]
  \item{orbitalVelMax},            maximal orbital velocity,        [m/s]
  \item{tidalVelMean},                 mean tidal velocity,        [m/s]
  \item{tidalVelMax},              maximal tidal velocity,        [m/s]
  \item{bedstress},                    bed shear stress,         [Pa]
  \item{EUNIScode},                           EUNIScode,          [-]
  \item{DRB},                    swept area ratio for dredge, [m2/m2/year]
  \item{OT},               swept area ratio for otter trawl, [m2/m2/year]
  \item{SN},                    swept area ratio for seines, [m2/m2/year]
  \item{TBB},                swept area ratio for beam trawl, [m2/m2/year]
  \item{sar}, total swept area ratio , DRB +OT +SN +TBB, [m2/m2/year]
  \item{subsar},  swept area ratio (fisheries) > 2cm,  [m2/m2/year]
  \item{gpd},              gear penetration depth,         [cm]  
 }
 
 The fishing data (sar, subsar, gpd) were derived from ICES upon request by OSPAR.
 
 They are averaged over (2009-2018). See also MWTL$sar.

==================

**\code{MWTL$fishing}**: 

species trait values that can be used to estimate fishing parameters.

==================

**\code{MWTL$types}**: 
the abiotic conditions of sampling stations, typologies 

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
 
==================

 **\code{MWTL$sar}**: 
 
 Fishing intensity data for the MWTL stations, origin: ICES upon request by OSPAR.  
 The NSBS stations nearest to the ICES data were selected.
 
 A data.frame that contains:
   \itemize{
      \item{station}, the NSBS station name
      \item{year}: the fishing year
      \item{gear}: metier; TBB, OT: beam, otter trawl; DRB: dredge, SN: seine.
      \item{sar}: annual swept area ratios (m2/m2/yr) for the surface (0-2cm).
      \item{subsar}: annual swept area ratios (m2/m2/yr) for the subsurface (>2cm).
      \item{gpd}: estimated gear penetration depths ([cm]), based on metier
    }
}

\note{
The dataset **\code{Taxonomy}**: 
contains taxonomic information of the original and adjusted taxon in MWTL$density, 
as derived from the World Register of Marine Species (WoRMS), using R-package \code{worrms}.
}
\details{
 The *macrobenthos data* of the Northsea (MWTL) are commissioned by the 
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

=====================================================================

The *fishing data* (in abiotics and sar) were derived from ICES upon request by OSPAR.
 
They are averaged over (2009-2018). 

The metiers are Aggregated into beam trawl (TBB), dredge (DRB), demersal seine (SN), and otter trawl (OT), based on the 
metier layers: OT_CRU, OT_DMF, OT_MIX, OT_MIX_CRU, OT_MIX_DMF_BEN, OT_MIX_DMF_PEL, OT_MIX_CRU_DMF,
OT_SPF, TBB_CRU, TBB_DMF, TBB_MOL, DRB_MOL, SDN_DMF, SSC_DMF

The gear penetration depths used were: 

\itemize{
 \item in sand: 3.5, 1.1, 1.1, 1.9 cm for DRB, OT, SN, TBB respectively
 \item in mud : 5.4, 2.0, 2.0, 3.2 cm for DRB, OT, SN, TBB respectively
}

The sediment was assumed to be "mud" when the type "sediment" was "Muddy"

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
metadata(MWTL$sar)

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

with(MWTLabiot, map_key(x, y, colvar = sand, 
                        main = "sand fraction", 
                        pch = 16))

# mud, plotted on large Northsea map
with(MWTLabiot, map_key(x, y, colvar = mud, 
                        contours = NSBS$contours, 
                        main = "mud fraction", 
                        pch = 16))

# show the different abiotic data sets
metadata(MWTL$abiotics)

##-----------------------------------------------------
## COMBINATIONS
##-----------------------------------------------------

NSsp_abi <- merge(MWTL$density, MWTL$sediment)
ECH      <- subset(NSsp_abi, 
                   subset = taxon=="Echinocardium")

with(ECH, points2D(D50, density, colvar = density, 
                   log = "xc", pch = 16))

# add station coordinates
ECH <- merge(ECH, MWTL$stations)

##-----------------------------------------------------
## From long format to wide format (stations x species)
##-----------------------------------------------------

NSwide <- with (MWTL$density, 
     l2w_density(descriptor  = station, 
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
trait.cwm <- get_trait_density (wide           = NSwide, 
                               trait          = Traits_nioz, 
                               taxonomy       = Taxonomy,
                               trait_class    = NStrait.lab$trait, 
                               trait_score    = NStrait.lab$score, 
                               scalewithvalue = TRUE)

head(trait.cwm, n=c(3,4))  

Stations.traits <- merge(MWTL$stations,  trait.cwm, 
                         by.x = "station", by.y = "descriptor")

##-----------------------------------------------------
## Maps
##-----------------------------------------------------

par(mfrow=c(2,2))

with(Stations.traits, 
  map_key(x, y, colvar = Biodeposition,
          main = "Biodeposition", 
          pch = 16))
          
with(Stations.traits, 
  map_key(x, y, colvar = Biodiffusion,
          main = "Biodiffusion", 
          pch = 16))
          
with(Stations.traits, 
  map_key(x, y, colvar = Biostabilisation,
          main = "Biostabilisation", 
          pch = 16))
with(Stations.traits, 
  map_key(x, y, colvar = Burrow.width,
          main = "Burrow width", 
          pch = 16))

##-----------------------------------------------------
## Show the depth contours
##-----------------------------------------------------

map_key(contours = MWTL$contours, draw.levels = TRUE, key.levels = TRUE)

# Use a different color scheme
collev <- function(n) 
   c("black",  ramp.col(col = c("darkgreen", "darkblue"), 
                        n = n-1))
           
map_key(contours = MWTL$contours, 
        draw.levels = TRUE, col.levels = collev,
        key.levels = TRUE)
        

##-----------------------------------------------------
## Fishing data
##-----------------------------------------------------

metadata(MWTL$sar)

# Sum fishing per year, per station

MWTLfish <- tapply(X    = MWTL$sar$sar, 
                  INDEX = list(MWTL$sar$station, MWTL$sar$year), 
                  FUN   = sum)

matplot(x = as.double(colnames(MWTLfish)), 
        y = t(MWTLfish), 
        xlab = "year", ylab = "m2/m2/yr", 
        main = "fishing intensity MWTL stations",
        type = "l", log= "y")
}

\references{

The taxonomic information was created using the \code{worrms} package:

Chamberlain S, Vanhoorne. B (2023). worrms: World Register of Marine Species (WoRMS) Client_. R package
  version 0.4.3, <https://CRAN.R-project.org/package=worrms>.

The MWTL data are described in:

L. Leewis, E.C. Verduin, R. Stolk ; Eurofins AquaSense Macrozoobenthosonderzoek in de Rijkswateren met boxcorer, jaarrapportage MWTL 2015 : waterlichaam: Noordzee
Publicatiedatum: 31-03-201775 p. Projectnummer Eurofins AquaSense: J00002105. 
Revisie 2, In opdracht van Ministerie van Infrastructuur en Milieu,  Rijkswaterstaat Centrale Informatievoorziening (RWS, CIV)

The fishing data were derived from:

ICES Technical Service, Greater North Sea and Celtic Seas Ecoregions, 29 August 2018 sr.2018.14 Version 2: 22 January 2019 https://doi.org/10.17895/ices.pub.4508 OSPAR request on the production of spatial data layers of fishing intensity/pressure.
}

\seealso{
  \link{map_key} for plotting.

  \link{Traits_nioz} for the trait datasets.
  
  \link{get_density} for functions operating on these data.

  \link{long2wide} for functions changing the appearance on these data.

  \link{get_Db_index} for extracting bioturbation and bioirrigation indices.
}

\keyword{datasets}

