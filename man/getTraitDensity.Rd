\name{Combining Density and Trait data}
\alias{getTraitDensity}
\title{
  Functions for combining station x taxon and trait data.
}

\description{
 \code{getTraitDensity} combines (descriptor x taxon) density data with (taxon x trait) data to a (descriptor x trait) data set.
}

\usage{
getTraitDensity(descriptor, taxon, value, averageOver = NULL, 
        wide = NULL, d.column = 1, trait, t.column = 1, 
        trait.class = NULL, trait.score = NULL, taxonomy = NULL, 
        scalewithvalue = TRUE, verbose = FALSE)
}

\arguments{
  \item{descriptor }{name(s) of the descriptor, i.e. *where* the data were taken, e.g. station names. Either a vector, a list, a data.frame or matrix (and with multiple columns). It can be of type numerical, character, or a factor. When a data.frame or a list the "names" will be used in the output; when a vector, the argument name will be used in the output.
}    
  \item{taxon }{vector describing *what* the data are; it gives the taxonomic name (e.g. species). Should be of the same length as (the number of rows of) \code{descriptor}. Can be a list (or data.frame with one column), or a vector. When a data.frame or a list the "name" will be used in the output; when a vector, the argument name will be used.
}
  \item{value }{vector or list that contains the *values* of the data, usually density. Should be of the same length as (the number of rows of) \code{descriptor} and \code{taxon}. For function \code{getDensity}, \code{value} can also be a multi-column data.frame or matrix. 
  }
  \item{averageOver }{vector with *replicates* over which averages need to be taken. Should be of the same length as (the number of rows of) \code{descriptor}.
  }
  \item{wide }{density data, in *WIDE* format. This is a data.frame or matrix with (descriptor x taxon) information. If \code{NULL}, this data.frame will be calculated from the \code{descriptor}, \code{taxon}, (\code{replicate}) and \code{value} data. The number of descriptor columns are specified with \code{d.column}. If a \code{data.frame} then the first column usually contains the descriptor name; in this case the dimensions of \code{wide} are (number of descriptors x number of species+1), and \code{d.column=1}. It is also allowed to have the descriptors as row.names of the data.frame -this requires setting \code{d.column=0}.   
  }
  \item{trait }{(taxon x trait) data or (descriptor x trait) data, in *WIDE* format, and containing numerical values only. Traits can be fuzzy coded. The number of columns with taxonomic information is specified with \code{t.column}. The default is to have the first column contain the name of the taxon, and \code{t.column=1}. It is also allowed to have the taxa as \code{row.names} of the data.frame; in this case, \code{t.column=0}. See last example for how to deal with cgtegorical traits.
  }
  \item{trait.class }{indices to trait levels, a vector. The length of this vector should equal the number of columns of \code{trait} or \code{wide} minus the value of \code{t.column}. If present, this -together with \code{trait.score}- will be used to convert the trait matrix from fuzzy to crisp format.
  }
  \item{trait.score }{trait values or scores, a vector. Should be of same length as \code{trait.class}
  }
  \item{scalewithvalue }{when TRUE, will standardize with respect to total density, so that the average trait value is obtained (not the summed value). Note that total density will be estimated only for those taxa whose traits are known.
  }
  \item{verbose }{when TRUE, will write warnings to the screen.
  }
  \item{taxonomy }{taxonomic information (the relationships between the taxa), a data.frame; first column will be matched with \code{taxon}, regardless of its name. This is used to estimate traits of taxa that are not accounted for, and that will be estimated based on taxa at the nearest taxonomic level. See details.
  }
  \item{d.column, t.column }{position(s) or name(s) of the column(s) that holds the descriptor of the density data set (data.frame \code{wide}), and the taxa in the trait data set (data.frame \code{trait}). The default is to have the first column holding the descriptors or taxa. If \code{NULL}, or \code{0}, then there is no separate column with names, so the \code{row.names} of the dataset are used as descriptor or taxon names.
  }
}

\details{

The taxonomy is used to fill in the gaps of the trait information, assuming that closely related taxa will share similar traits. This is done in two steps: 

The trait database is first extended with information on higher taxonomic levels (using *extendTrait*). The traits for a taxonomic level are estimated as the average of the traits at the lower level. For instance, traits on genus level will be averages of known traits of all species in the database belonging to this genus. 

Then, for each taxon that is not present in the trait database, the traits on the closest taxonomic level are used. For instance, for an  unrecorded species, it is first checked if the trait is known on genus level, if not, family level and so on.
}

\value{
\code{getTraitDensity} returns the descriptor x trait density matrix.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}
\seealso{
\link{MWTL} for the data sets.

\link{mapBtrait} for simple plotting functions.

\link{getDensity} for functions working with density data.

\link{getDbIndex} for extracting bioturbation and bioirrigation indices.

\link{extendTrait} for functions working with traits.

\link{getTrait} 
}

\note{
When traits of certain taxa are not found, they are put equal to 0 and 
calculation will proceed; this will be notified only when \code{verbose} is set 
to \code{TRUE}. 
When that happens, the taxa that were ignored can be found in \code{attributes(..)$notrait}
 
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

##-----------------------------------------------------
## From long to wide format, averaging over replicates
##-----------------------------------------------------

Bwide <- with (Bdata.rep, 
  l2wDensity (value      = density, 
             descriptor  = station, 
             averageOver = replicate,
             taxon       = species))
Bwide

##-----------------------------------------------------
## Small dataset: fuzzy-coded traits
##-----------------------------------------------------

# Note: no data for "sp.4"

Btraits <- data.frame(
  taxon   = c("sp.1","sp.2","sp.3","sp.5","sp.6"),
  T1_M1   = c(0     , 0    ,   0  , 0.2  ,     1),
  T1_M2   = c(1     , 0    , 0.5  , 0.3  ,     0),
  T1_M3   = c(0     , 1    , 0.5  , 0.5  ,     0),
  T2_M1   = c(0     , 0    ,   1  , 0.5  ,     1),
  T2_M2   = c(1     , 1    ,   0  , 0.5  ,     0)
)

# The metadata for this trait
Btraits.lab <- data.frame(
  colname  =c("T1_M1","T1_M2","T1_M3","T2_M1","T2_M2"),
  trait    =c("T1"   ,"T1"   ,"T1"   ,"T2"   ,"T2"),
  modality =c("M1"   ,"M2"   ,"M3"   ,"M1"   ,"M2"), 
  score    =c(0      , 0.5   , 1     , 0.2   , 2)
)

##-----------------------------------------------------
## Small dataset: taxonomy
##-----------------------------------------------------

Btaxonomy <- data.frame(
  species = c("sp.1","sp.2","sp.3","sp.4","sp.5","sp.6"),
  genus   = c( "g.1", "g.2", "g.2", "g.2", "g.3", "g.4"),
  family  = c( "f.1", "f.1", "f.1", "f.1", "f.2", "f.3"),
  order   = c( "o.1", "o.1", "o.1", "o.1", "o.2", "o.2"),
  class   = c( "c.1", "c.1", "c.1", "c.1", "c.1", "c.1")
  )

##-----------------------------------------------------
## Community weighted mean score
##-----------------------------------------------------

# if verbose=TRUE, a warning will be given, as there is no trait information
# for sp4 - this affects the trait density for station "b"

cwm.trait <- getTraitDensity (
                        wide           = Bwide,  
                        trait          = Btraits, 
                        trait.class    = Btraits.lab$trait, 
                        trait.score    = Btraits.lab$score, 
                        scalewithvalue = TRUE, 
                        verbose        = FALSE)

cwm.trait
attributes(cwm.trait)$notrait   # species that was ignored

# same, but total trait values
cwm.trait.2 <- getTraitDensity (
                        wide           = Bwide,  
                        trait          = Btraits, 
                        trait.class    = Btraits.lab$trait, 
                        trait.score    = Btraits.lab$score, 
                        scalewithvalue = FALSE, 
                        verbose        = FALSE)

cwm.trait.2

# Traits from all taxa in the dataset, including high-level information 

cwm.trait.2 <- getTraitDensity (
                        wide           = Bwide,  
                        trait          = Btraits, 
                        trait.class    = Btraits.lab$trait, 
                        trait.score    = Btraits.lab$score, 
                        taxonomy       = Btaxonomy,
                        scalewithvalue = TRUE)

cwm.trait.2
attributes(cwm.trait.2)$notrait

# Same but keeping fuzzy scores

cwm.trait.3 <- getTraitDensity (
                        wide           = Bwide,  
                        trait          = Btraits, 
                        taxonomy       = Btaxonomy,
                        scalewithvalue = TRUE)

cwm.trait.3

##-----------------------------------------------------
## categorical traits
##-----------------------------------------------------

# Note: no data for "sp.4"

Bcategory <- data.frame(
  taxon   = c("sp.1","sp.2","sp.3","sp.5","sp.6"),
  C1      = c(   "A",   "B",   "A",   "C",   "C")
)

\dontrun{
# this will not work, as trait should be numerical
cwm.trait.4 <- getTraitDensity (
                        wide           = Bwide,  
                        trait          = Bcategory, 
                        taxonomy       = Btaxonomy,
                        scalewithvalue = TRUE)
}

crisp2fuzzy(Bcategory)

# this will  work, as categorical traits -> numerical
cwm.trait.4 <- getTraitDensity (
                        wide           = Bwide,  
                        trait          = crisp2fuzzy(Bcategory), 
                        taxonomy       = Btaxonomy,
                        scalewithvalue = TRUE)
cwm.trait.4

}
