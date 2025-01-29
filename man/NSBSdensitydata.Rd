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
    \item Macrofauna species density and biomass (\code{NSBS$density}, \code{NSBS$biomass})
  
    \item Abiotic conditions (\code{NSBS$abiotics}), and station positions (\code{NSBS$stations})
  
    \item \code{NSBS$contours}: contourlines for mapping.
    
    \item NSBS$fishing: species trait values that can be used to estimate fishing parameters.
   
    \item NSBS$sar: station-specific fishing intensities.
    
  }
}

\usage{
data(NSBS)
}

\format{
==================

**\code{NSBS$density}**: 

This is the main Northsea NSBS benthos data set, containing species 
information for 234 stations sampled in 1985-1986.
 
The data, in long format, are in a \code{data.frame} with the following columns: 
  \itemize{
   \item \code{station}, the NSBS station name (details in NSBS$stations).    \item \code{date}, the sampling date (year). 
   \item \code{taxon}, the taxon name to be used (usually species), and checked against the WoRMS database (details in dataset \code{Taxonomy}).
   \item \code{density}, the number of individuals per m2.
  }

==================

**\code{NSBS$biomass}**: 

This data set, contains biomass data  for 214 stations sampled in 1985-1986. 

Biomass information is on higher taxonomic level.
The distinction is made between Crustacea, Polychaeta, Animalia, Mollusca, Echinodermata.
 
The data, in long format, are in a \code{data.frame} with the following columns: 
  \itemize{
   \item \code{station}, the NSBS station name (details in NSBS$stations).    \item \code{date}, the sampling date (year). 
   \item \code{taxon}, the taxon name to be used (5 higher taxa).
   \item \code{biomass}, the biomass, gAFDW/m2.
  }

==================

**\code{NSBS$abiotics}**: 
the abiotic conditions of sampling stations.

NSBS$abiotics is a data.frame with the following columns:   
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
 
 They are averaged over (2009-2018). See also NSBS$sar.

==================

**\code{NSBS$fishing}**: 

species trait values that can be used to estimate fishing parameters.

==================

**\code{NSBS$stations}**: 

The positions of the different stations, in WGS84 format
 \itemize{
  \item{station}, the NSBS station name
  \item{x}, degrees longitude
  \item{y}, degrees latitude
 }
 
==================

 **\code{NSBS$contours}**: 
 
 The data for mapping the contours. The contourlines (x-, y) were derived from 
 GEBCO high-resolution bathymetry, by creating contourlines.
 
 The data set contains:
   \itemize{
      \item{station}, the NSBS station name
      \item{x}: longitude, in [dgE]
      \item{y}: latitude, in [dgN]
      \item{z}: the corresponding depths, in [m]
    }

==================

 **\code{NSBS$sar}**: 
 
 Fishing data for the NSBS stations, origin: ICES upon request by OSPAR.  
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
contains taxonomic information of the original and adjusted taxon in NSBS$density, 
as derived from the World Register of Marine Species (WoRMS), using R-package \code{worrms}.
}

\details{
NSBS contains the *macrofauna data* from the 1986 North Sea Benthos Survey, an activity of the Benthos Ecology Working Group of ICES. 

Benthic samples were taken in a standardised way, on a regular grid covering the whole of the North Sea, and analysed by scientists from 10 laboratories. Extensive work was done to standardise taxonomy and identifications across the different laboratories.

Sediment was sampled with a Reineck Boxcorer (0,078 m2). Macrofauna sieved on a 1 mm mesh. 

=====================================================================

The *fishing data* (in abiotics and sar) were derived from ICES upon request by OSPAR.

They are averaged over (2009-2018). 

The metiers are Aggregated into beam trawl (TBB), dredge (DRB), demersal seine (SN), and otter trawl (OT), based on the 
metier layers: OT_CRU, OT_DMF, OT_MIX, OT_MIX_CRU, OT_MIX_DMF_BEN, OT_MIX_DMF_PEL, OT_MIX_CRU_DMF,
OT_SPF, TBB_CRU, TBB_DMF, TBB_MOL, DRB_MOL, SDN_DMF, SSC_DMF

The gear penetration depths used were the mean values over sand and muddy sediments: 

\itemize{
 \item in sand: 3.5, 1.1, 1.1, 1.9 cm for DRB, OT, SN, TBB respectively
 \item in mud : 5.4, 2.0, 2.0, 3.2 cm for DRB, OT, SN, TBB respectively
}


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

with(NSBSab, 
  map_key(x, y, colvar = mud, 
          contours = NSBS$contours, 
          main = "mud fraction", 
          pch = 16))

metadata(NSBS$abiotics)

##-----------------------------------------------------
## COMBINATIONS
##-----------------------------------------------------

NSsp_abi <- merge(NSBS$density, NSBS$abiotics)

ECH      <- subset(NSsp_abi, 
                   subset = taxon=="Echinocardium cordatum")

with(ECH, 
  plot(mud, density, 
       main = "E. cordatum", 
       xlab = "mud fraction", ylab = "density, ind/m2", 
       pch = 16))

# add station coordinates
ECH <- merge(ECH, NSBS$stations)

##-----------------------------------------------------
## From long format to wide format (stations x species)
##-----------------------------------------------------

NSwide <- with (NSBS$density, 
         l2w_density(descriptor = station,    # long2wide for density
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

trait.cwm <- get_trait_density(wide           = NSwide, 
                               trait          = Traits_nioz, 
                               taxonomy       = NSBS$taxonomy,
                               trait.class    = Trait.lab$trait, 
                               trait.score    = Trait.lab$score, 
                               scalewithvalue = TRUE)

head(trait.cwm, n = c(3, 4))  

Stations.traits <- merge(NSBS$stations, trait.cwm, 
                         by.x = "station", by.y = "descriptor")

par(mfrow=c(2,2))

with(Stations.traits, 
  map_key(x, y, colvar = Biodeposition,
          main = "Biodeposition"))
         
with(Stations.traits, 
  map_key(x, y, colvar = Biodiffusion,
          main = "Biodiffusion"))
         
with(Stations.traits, 
  map_key(x, y, colvar = Biostabilisation,
          main = "Biostabilisation"))
         
with(Stations.traits, 
  map_key(x, y, colvar = Burrow.width,
          main = "Burrow width"))

##-----------------------------------------------------
## Community mean weighted score for typological groups.
##-----------------------------------------------------

# Groups is in crisp format -> convert to fuzzy 

Groups.fuz <- crisp2fuzzy(Groups[,c("taxon", "typology")])

head (Groups, n = 2)
head (Groups.fuz, n = c(2, 5))

group.cwm <- get_trait_density(wide           = NSwide, 
                               trait          = Groups.fuz, 
                               scalewithvalue = TRUE)

head(group.cwm, n=c(3,4))  

summary(group.cwm)

##-----------------------------------------------------
## Show the depth contours
##-----------------------------------------------------

map_key(contours = NSBS$contours, 
        draw.levels = TRUE, key.levels = TRUE)

# Use a different color scheme
collev <- function(n) 
      c("black", ramp.col(col = c("darkgreen", "darkblue"), 
                          n = n-1))
map_key(contours = NSBS$contours, 
        draw.levels = TRUE, col.levels = collev,
        key.levels = TRUE)

##-----------------------------------------------------
## Fishing data
##-----------------------------------------------------

metadata(NSBS$sar)

# Sum fishing per year, per station

NSBSfish <- tapply(X    = NSBS$sar$sar, 
                  INDEX = list(NSBS$sar$station, NSBS$sar$year), 
                  FUN   = sum)

matplot(x = as.double(colnames(NSBSfish)), 
        y = t(NSBSfish), 
        xlab = "year", ylab = "m2/m2/yr", 
        main = "fishing intensity NSBS stations",
        type = "l", log= "y")
}

\references{

The taxonomic information was created using the \code{worrms} package:

Chamberlain S, Vanhoorne. B (2023). worrms: World Register of Marine Species (WoRMS) Client_. R package
  version 0.4.3, <https://CRAN.R-project.org/package=worrms>.

The NSBS data are described in:

Heip, C.H.R.; Basford, D.; Craeymeersch, J.A.; Dewarumez, J.-M.; Dorjes, J.; de Wilde, P.; Duineveld, G.; Eleftheriou, A.; Herman, P.M.J.; Kingston, K.; Niermann, U.; Kunitzer, A.; Rachor, E.; Rumohr, H.; Soetaert, K.; Soltwedel, T. (1992). Trends in biomass, density and diversity of North Sea macrofauna. ICES J. Mar. Sci./J. Cons. int. Explor. Mer 49: 13-22

The fishing data were derived from:

ICES Technical Service, Greater North Sea and Celtic Seas Ecoregions, 29 August 2018 sr.2018.14 Version 2: 22 January 2019 https://doi.org/10.17895/ices.pub.4508 OSPAR request on the production of spatial data layers of fishing intensity/pressure.
}

\seealso{
  \link{map_key} for plotting.

  \link{Traits_nioz} for the trait datasets.
  
  \link{get_density} for functions operating on these data.

  \link{get_Db_index} for extracting bioturbation and bioirrigation indices.

  \link{long2wide} for functions changing the appearance on these data.
}

\keyword{datasets}

