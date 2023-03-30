\name{Density functions}
\alias{addAbsences}
\alias{getProportion}
\alias{getDensity}
\title{
  Functions for working with station x taxon data.
}

\description{
 \code{getDensity} takes a subset of density data, possibly at higher taxonomic level.

 \code{getProportion} estimates proportional values from (a subset of) density data.

 \code{addAbsences} adds taxon absences in density data in long format.
}

\usage{
getDensity(descriptor, taxon, value, averageOver = NULL, 
        taxonomy = NULL, subset, wide.output = FALSE, 
        full.output = FALSE, verbose=FALSE)

getProportion(descriptor, taxon, value, averageOver = NULL, 
        taxonomy = NULL, verbose=FALSE)

addAbsences(descriptor, taxon, value, averageOver = NULL)
}

\arguments{
  \item{descriptor }{name(s) of the descriptor, i.e. *where* the data were taken, e.g. station names. Either a vector, a list, a data.frame or matrix (and with multiple columns). It can be of type numerical, character, or a factor. When a data.frame or a list the "names" will be used in the output; when a vector, the argument name will be used in the output. In theory, this could also be one number or \code{NA}; however, care needs to be taken in case this combined with \code{subset} and \code{averageOver}
}    
  \item{taxon }{vector describing *what* the data are; it gives the taxonomic name (e.g. species). Should be of the same length as (the number of rows of) \code{descriptor}. Can be a list (or data.frame with one column), or a vector. When a data.frame or a list the "name" will be used in the output; when a vector, the argument name will be used.
}
  \item{value }{vector or list that contains the *values* of the data, usually density. Should be of the same length as (the number of rows of) \code{descriptor} and \code{taxon}. For function \code{getDensity}, \code{value} can also be a multi-column data.frame or matrix. 
  }
  \item{averageOver }{vector with *replicates* over which averages need to be taken. Should be of the same length as (the number of rows of) \code{descriptor}.
  }
  \item{subset }{logical expression indicating elements to keep: missing values are taken as FALSE. If NULL, or absent, then all elements are used. Note that the subset is taken *after* the number of samples to average per descriptor is calculated, so this will also work for selecting certain taxa that may not be present in all replicates over which should be averaged.
  }
  \item{taxonomy }{taxonomic information; first column will be matched with \code{taxon}, regardless of its name.
  }
  \item{full.output }{ when TRUE, will also return descriptors for which the value is 0. This can be relevant in case a selection is made for taxonomic composition. Note that the taxon in this case will be undefined.
  }
  \item{wide.output }{ when TRUE, will recast the output in wide format (the default is long format). This only makes sense when \code{value} has multiple columns. In the wide format, each row holds the taxa for a descriptor, and each column holds the values for one taxon and all descriptors.
  }
  \item{verbose }{ when TRUE, may write warnings to the screen.
  }  
}

\value{
\itemize{
 \item \code{getDensity} returns a \code{data.frame} combining descriptor, taxon, and value, in LONG format.
 
 Unless argument \code{full.output} is \code{TRUE}, this will not contain absences 
 (i.e. where value = 0), unless they were already present in the input.

  The \code{value} consists of summed values for (\code{descriptor}, \code{averageOver}, 
  \code{taxon}) combinations, followed by averaging over \code{averageOver}.
}

}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}
\seealso{
\link{MWTL} for the data sets.

\link{mapKey} for simple plotting functions.

\link{long2wide} for estimating converting data from long to wide format and vice versa.

\link{getSummary} for estimating summaries from density data.

\link{getTraitDensity} for functions combining density and traits.

\link{extendTrait} for functions working with traits.

\link{getDbIndex} for extracting bioturbation and bioirrigation indices.
}

\details{

\itemize{
 \item{\code{getDensity}}. In this function a summed \code{value} over 
    \code{taxa} x \code{descriptor} is first calculated, and this value is then 
    divided by the number of replicates (as in \code{averageOver}) per \code{descriptor}.
    
    This is suitable for density and biomass data, but NOT for mean individual weight for instance.
}

}
\keyword{manip}

\examples{

## ====================================================
## A small dataset with replicates
## ====================================================

# 2 stations, 2 replicates for st.a, one replicate for st.b
Bdata.rep <- data.frame(
  station   = c("st.a","st.a","st.a","st.b","st.b","st.b"),
  replicate = c(     1,     1,    2,     1,     1,     1),
  species   = c("sp.1","sp.2","sp.1","sp.3","sp.4","sp.5"),
  density   = c(     1,     2,    3,     3,     1,     3)
)
Bdata.rep

# Proportion of each species
with (Bdata.rep, 
  getProportion (value       = density, 
                 descriptor  = data.frame(station, replicate), 
                 taxon       = species))


##-----------------------------------------------------
## average of replicates
##-----------------------------------------------------

with (Bdata.rep, 
  getDensity (value       = density, 
              descriptor  = station, 
              taxon       = species, 
              averageOver = replicate))

# using named lists to have good column headings
with (Bdata.rep, 
  getDensity (value       = list(density=density), 
              descriptor  = list(station=station), 
              taxon       = species, 
              averageOver = replicate))
                 
# estimating proportions
with (Bdata.rep, 
  getProportion (value    = list(density=density), 
              descriptor  = list(station=station), 
              taxon       = species, 
              averageOver = replicate))
                 
# averaging multiple value columns at once 
# extending the data with biomass - assume no biomass for st b

Bdata.rep$biomass = c(0.1, 0.2, 0.3, NA, NA, NA)
Bdata.rep

with (Bdata.rep, 
  getDensity (value       = list(density=density, biomass=biomass), 
              descriptor  = list(station=station), 
              taxon       = species, 
              averageOver = replicate))

##-----------------------------------------------------
## Select information for one species
##-----------------------------------------------------

with (Bdata.rep, 
  getDensity (subset      = species=="sp.2",
              value       = density, 
              descriptor  = list(station=station), 
              taxon       = species, 
              averageOver = replicate))

# returns also the 0 value 
with (Bdata.rep, 
  getDensity (subset      = species=="sp.2",
              value       = density, 
              descriptor  = list(station=station), 
              taxon       = species, 
              averageOver = replicate, 
              full.output = TRUE))

##-----------------------------------------------------
## Extend the long format with absences
##-----------------------------------------------------

# take averages over replicates
with(Bdata.rep, 
  addAbsences (value      = density, 
               descriptor = cbind(station, replicate), 
               taxon      = species))

# keep replicates
with(Bdata.rep, 
  addAbsences (value       = density, 
               descriptor  = station, 
               taxon       = species, 
               averageOver = replicate))

## ====================================================
## A small dataset without replicates
## ====================================================

Bdata <- data.frame(
  station = c("st.a","st.a","st.b","st.b","st.b","st.c"),
  species = c("sp.1","sp.2","sp.1","sp.3","sp.4","sp.5"),
  density = c(1, 2, 3, 3, 1, 3)
)

## ====================================================
## Small dataset: taxonomy
## ====================================================

Btaxonomy <- data.frame(
  species = c("sp.1","sp.2","sp.3","sp.4","sp.5","sp.6"),
  genus   = c( "g.1", "g.2", "g.2", "g.2", "g.3", "g.4"),
  family  = c( "f.1", "f.1", "f.1", "f.1", "f.2", "f.3"),
  order   = c( "o.1", "o.1", "o.1", "o.1", "o.2", "o.2"),
  class   = c( "c.1", "c.1", "c.1", "c.1", "c.1", "c.1")
  )

##-----------------------------------------------------
## density on higher taxonomic level
##-----------------------------------------------------

# species density for a particular genus 
sp.g2 <- with (Bdata, 
  getDensity(descriptor = station, 
             taxon      = species,
             value      = density,
             taxonomy   = Btaxonomy,
             subset     = genus == "g.2")
              )

# select data for station st.a
Bselect <- with (Bdata, 
   getDensity (value      = density, 
               descriptor = station, 
               taxon      = species, 
               subset     = station=="st.a")
                )
Bselect

# pass taxonomy to select only species that belong to g.1
with (Bdata, 
   getDensity (value      = density, 
               descriptor = station, 
               taxon      = species, 
               taxonomy   = Btaxonomy, 
               subset     = genus=="g.1"))

## ====================================================
## Northsea dataset
## ====================================================

#-----------------------------------------------------
## Occurrence of Abra alba, averaged per station
##-----------------------------------------------------

Abra_alba <- with(MWTL$density, 
   getDensity(subset      = taxon=="Abra alba",
              descriptor  = station,
              averageOver = year,
              taxon       = taxon,
              value       = density))
head(Abra_alba)

# Mean of all Abra species over all stations
# This should be done in two steps.

Abra <- with(MWTL$density, 
   getDensity(subset      = genus=="Abra",
              descriptor  = station,
              averageOver = year,
              taxon       = taxon,
              taxonomy    = Taxonomy,
              value       = density, 
              full.output = TRUE))
head(Abra)

tapply(Abra$value, INDEX=list(Abra$taxon), FUN=mean)

}
