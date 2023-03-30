\name{Summary statistics for density}
\alias{getSummary}
\title{
  get Summary for station x taxon data (getSummary)
}

\description{
 \code{getSummary} calculates summaries of density data, i.e. total densities and number of taxa.
}

\usage{
getSummary(descriptor, taxon, value, averageOver = NULL, taxonomy = NULL, 
        subset, what=c("density", "taxa", "occurrence"), 
        wide.output = FALSE)
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
  \item{taxonomy }{taxonomic information; first column will be matched with \code{taxon}, regardless of its name.
  }
  \item{subset }{logical expression indicating elements to keep: missing values are taken as FALSE. If NULL, or absent, then all elements are used. Note that the subset is taken *after* the number of samples to average per descriptor is calculated, so this will also work for selecting certain taxa that may not be present in all replicates over which should be averaged.
  }
  \item{what }{the summary statistic to be returned, one of "density": total density per descriptor, "taxa": the number of taxa per descritpor, "occurrence": the number of occurrences of each species.
  }
  \item{wide.output }{when TRUE, will recast the output in wide format (the default is long format). This only makes sense when \code{descriptor} is a matrix or data.frame with more > 1 column; in this case, the last column will make up the columns.+
  }
}

\value{
 \code{getSummary} returns a list, containing the mean densities (\code{density}) 
 and the number of taxa (\code{taxa}) for each descriptor, and the number of descriptors in which the 
 taxa were found (\code{occurrence}).
 
 In case \code{wide.output} is \code{FALSE}, this will be as a data.frame in long 
 format (descriptor, value); when \code{TRUE}, the descriptor will feature both 
 in the rows (all bot the last descriptor) and in the columns (the last descriptor).

}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}
\seealso{.
\link{MWTL} for the data sets

\link{mapKey} for simple plotting functions.

\link{getDensity} for functions working with density data.

\link{getTraitDensity} for functions combining density and traits.
 
\link{getDbIndex} for extracting bioturbation and bioirrigation indices.

\link{extendTrait} for functions working with traits.

\link{getTrait} 
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
## Summary statistics
##-----------------------------------------------------

with (Bdata.rep,  
  getSummary(descriptor  = cbind(station, replicate), 
             value       = density, 
             taxon       = species, 
             wide.output = FALSE))

with (Bdata.rep,  
  getSummary(descriptor  = station, 
             averageOver = replicate,
             value       = density, 
             taxon       = species, 
             wide.output = FALSE))

## ====================================================
## Northsea dataset
## ====================================================

NSsumm <- with(MWTL$density, 
   getSummary(descriptor  = station,
              averageOver = year,
              taxon       = taxon,
              value       = density))
head(NSsumm$density)
head(NSsumm$taxa)
head(NSsumm$occurrence)

}
