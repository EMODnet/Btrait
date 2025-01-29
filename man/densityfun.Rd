\name{Density functions}
\alias{add_absences}
\alias{get_proportion}
\alias{get_density}
\title{
  Functions for working with station x taxon data.
}

\description{
 \code{get_density} takes a subset of density data, possibly at higher taxonomic level.

 \code{get_proportion} estimates proportional values from (a subset of) density data.

 \code{add_absences} adds taxon absences in density data in long format.
}

\usage{
get_density(data, descriptor, taxon, value, averageOver, 
            taxonomy = NULL, subset, wide.output = FALSE, 
            full.output = FALSE, verbose=FALSE)

get_proportion(data, descriptor, taxon, value, averageOver, 
            taxonomy = NULL, verbose=FALSE)

add_absences(data, descriptor, taxon, value, averageOver)
}

\arguments{
  \item{data }{data.frame to use for extracting the arguments \code{descriptor}, \code{taxon}, \code{value}, \code{averageOver}. Can be missing. 
  }
  \item{descriptor }{ variable(s) *where* the data were taken, e.g. sampling stations. If \code{data} is not missing: one or more column(s) from \code{data}; use \code{cbind} or \code{data.frame} to select more columns. If \code{data} is missing: a vector, a list, a data.frame or a matrix (with one or multiple columns). It can be of type numerical, character, or a factor. In theory, descriptor can also be one number, \code{NA} or missing; however, care needs to be taken in case this combined with \code{subset} and \code{averageOver}.
}    
  \item{taxon }{variables describing *what* the data are;  it gives the taxonomic name (e.g. species). If \code{data} is not missing:  one column from \code{data}.
If \code{data} is missing:  a list (or data.frame with one column), or a vector. When a data.frame or a list the "name" will be used in the output; when a vector, the argument name will be used.
}
  \item{value }{variable that contains the *values* of the data, usually density. If \code{data} is not missing:  one or more column(s) from \code{data}; use \code{cbind} or \code{data.frame} to select more columns. If \code{data} is missing: a vector, a list, a data.frame or a matrix (with one or multiple columns). it should be of the same length (or have the same number of rows) as (the number of rows of) \code{descriptor} and \code{taxon}. Should contain numerical values. Should always be present.
  }
  \item{averageOver }{*replicates* over which averages need to be taken. If \code{data} is not missing:  one or more column(s) from \code{data}; use \code{cbind} or \code{data.frame} to select more columns. Else a vector, a list, a data.frame or a matrix (with one or multiple columns). It can be of type numerical, character, or a factor. Can be absent.
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
 \item \code{get_density} returns a \code{data.frame} combining descriptor, taxon, and value, in LONG format.
  The \code{value} consists of summed values for (\code{descriptor}, \code{averageOver}, 
  \code{taxon}) combinations, followed by averaging over \code{averageOver}.
 
 Unless argument \code{full.output} is \code{TRUE}, this will not contain absences 
 (i.e. where value = 0), unless they were already present in the input.
 
 \item \code{get_proportion} returns a \code{data.frame} combining descriptor, taxon, and proportional, in LONG format.
}

Depending on whether argument \code{data} is passed or not, 
the output columns may be labelled differently:

\itemize{
\item if \code{data} is passed: the original names in data will be kept
\item if \code{data} is not passed: the names will only be kept if explicitly passed.
}

}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}
\seealso{
\link{MWTL} for the data sets.

\link{map_key} for simple plotting functions.

\link{long2wide} for estimating converting data from long to wide format and vice versa.

\link{get_summary} for estimating summaries from density data.

\link{get_trait_density} for functions combining density and traits.

\link{extend_trait} for functions working with traits.

\link{get_Db_index} for extracting bioturbation and bioirrigation indices.
}

\details{

\itemize{

 \item{\code{get_density}}. In this function a summed \code{value} over 
    \code{taxa} x \code{descriptor} is first calculated, and this value is then 
    divided by the number of replicates (as in \code{averageOver}) per \code{descriptor}.
    
    This is suitable for density and biomass data, but NOT for mean individual weight for instance (which cannot be summed).
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

## ====================================================
## Proportion of each species -
## ====================================================

# three ways to extract proportions

# -----------------------------------------------------
# 1. using input argument data 
# -----------------------------------------------------
# note: taxon is called species as in original data set
  
  get_proportion (data        = Bdata.rep,
                  value       = density, 
                  descriptor  = data.frame(station, replicate), 
                  taxon       = species)

  # to force the taxon to be called taxon
  get_proportion (data        = Bdata.rep,
                  value       = density, 
                  descriptor  = data.frame(station, replicate), 
                  taxon       = data.frame(taxon = species))  # named dataframe

# -----------------------------------------------------
# 2. using with() to create an environment 
# -----------------------------------------------------
# note: taxon is called taxon

with (Bdata.rep, 
  get_proportion (value       = density, 
                  descriptor  = data.frame(station, replicate), 
                  taxon       = species)
      )

# force the taxon to be called "species"
with (Bdata.rep, 
  get_proportion (value       = density, 
                  descriptor  = data.frame(station, replicate), 
                  taxon       = data.frame(species = species))
      )          

# -----------------------------------------------------
# 3. using pipes
# -----------------------------------------------------
  Bdata.rep |> get_proportion (value       = density, 
                               descriptor  = data.frame(station, replicate), 
                               taxon       = species)


##-----------------------------------------------------
## average of replicates
##-----------------------------------------------------

  PP <-  get_density(
              data        = Bdata.rep,
              value       = density, 
              descriptor  = station, 
              taxon       = species, 
              averageOver = replicate)


# input arguments are kept in the attributes              
attributes(PP)[-(1:3)]
                 
# averaging multiple value columns at once 
# extending the data with biomass - assume no biomass for st b

Bdata.rep$biomass = c(0.1, 0.2, 0.3, NA, NA, NA)
Bdata.rep

  DD <-  get_density(
              data        = Bdata.rep,
              value       = data.frame(density, biomass), 
              descriptor  = data.frame(station), 
              taxon       = species, 
              averageOver = replicate)

# input arguments              
  attributes(DD)[-(1:3)]

##-----------------------------------------------------
## Select information for one species
##-----------------------------------------------------

  get_density(data        = Bdata.rep, 
              subset      = species=="sp.2",
              value       = density, 
              descriptor  = station, 
              taxon       = species, 
              averageOver = replicate)

# returns also the 0 value 
  get_density(data        = Bdata.rep,
              subset      = species=="sp.2",
              value       = density, 
              descriptor  = station, 
              taxon       = species, 
              averageOver = replicate, 
              full.output = TRUE)

##-----------------------------------------------------
## Extend the long format with absences
##-----------------------------------------------------

# take averages over replicates
Bdata.rep |>
  add_absences (value      = density, 
                descriptor = cbind(station, replicate), 
                taxon      = species)

# keep replicates
  add_absences (data        = Bdata.rep, 
                value       = density, 
                descriptor  = station, 
                taxon       = species, 
                averageOver = replicate)

## ====================================================
## ====================================================
## A small dataset without replicates
## ====================================================
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

# all input:                                 
# ** simple density

get_density(Bdata, 
            descriptor = station, 
            taxon      = species, 
            value      = density)

# without taxon: 
# ** sum per station

get_density(Bdata, 
            descriptor = station, 
            value      = density)

# without taxon, and averaging over species: 
# ** average per station

get_density(Bdata, 
            descriptor  = station, 
            averageOver = species,
            value       = density)

# without descriptor: 
# ** sum per species

get_density(Bdata, 
            taxon      = species, 
            value      = density)

# without descriptor, averaging over stations: 
# ** average per species

get_density(Bdata, 
            taxon       = species, 
            averageOver = station,
            value       = density)

# without descriptor and taxon: 
# ** sum of all

get_density(Bdata, 
            value       = density)

# without descriptor and taxon, averaging over species and stations: 
# ** average of all

get_density(Bdata, 
            averageOver = cbind(station, species),
            value       = density)

##-----------------------------------------------------
## density on higher taxonomic level
##-----------------------------------------------------

# species density for a particular genus 
sp.g2 <- 
  get_density(data       = Bdata,
              descriptor = station, 
              taxon      = species,
              value      = density,
              taxonomy   = Btaxonomy,
              subset     = genus == "g.2")
sp.g2

# select data for station st.a
Bselect <- 
   get_density(data       = Bdata,
               value      = density, 
               descriptor = station, 
               taxon      = species, 
               subset     = station=="st.a")

Bselect

# pass taxonomy to select only species that belong to g.1
   get_density(data       = Bdata,
               value      = density, 
               descriptor = station, 
               taxon      = species, 
               taxonomy   = Btaxonomy, 
               subset     = genus=="g.1")

## ====================================================
## Northsea dataset
## ====================================================

#-----------------------------------------------------
## Occurrence of Abra alba, averaged per station
##-----------------------------------------------------

Abra_alba <- 
   get_density(data        = MWTL$density,
               subset      = taxon=="Abra alba",
               descriptor  = station,
               averageOver = year,
               taxon       = taxon,
               value       = density)
head(Abra_alba)

# Mean of all Abra species over all stations
# This should be done in two steps.

Abra <- 
   get_density(data        = MWTL$density,
               subset      = genus=="Abra",
               descriptor  = station,
               averageOver = year,
               taxon       = taxon,
               taxonomy    = Taxonomy,
               value       = density, 
               full.output = TRUE)
head(Abra)

tapply(Abra$density, INDEX=list(Abra$taxon), FUN=mean)

}
