\name{Conversion functions}
\alias{long2wide}
\alias{wide2long}
\alias{l2wDensity}
\alias{l2wTrait}
\alias{w2lDensity}
\alias{w2lTrait}
\title{
  Functions for conversion from long to wide format and vice versa.
}

\description{
 \code{long2wide} casts data from long to wide format.

 \code{w2lDensity} casts density data from wide to long format.

 \code{w2lTrait} casts trait data from wide to long format.

 \code{wide2long} casts data from wide to long format.

 \code{l2wDensity} casts density data from long to wide format.

 \code{l2wTrait} casts trait data from long to wide format.

}

\usage{
wide2long(wide, d.column = 1, w.names = NULL, 
          absences = FALSE)

w2lDensity(wide, d.column = 1, taxon.names = NULL, 
           absences = FALSE) 

w2lTrait(wide, t.column = 1, trait.names = NULL, 
           absences = FALSE) 

long2wide(row, column, value, averageOver = NULL, 
        taxonomy = NULL, subset)
        
l2wDensity(descriptor, taxon, value, averageOver = NULL, 
        taxonomy = NULL, subset)
        
l2wTrait(trait, taxon, value, averageOver = NULL, 
        taxonomy = NULL, subset)
}

\arguments{
  \item{wide }{data, in *WIDE* format. For density data, this is a data.frame or matrix with (descriptor x taxon) information (for density data), and The first column usually contains the descriptor name. For trait data this is a data.frame with (taxon x trait) information, and the first column generally contains the names of the taxa. It is also allowed to have the descriptors as row.names of the data.frame -this requires setting \code{d.column=0}.   
  }
  \item{row }{vector or data.frame that contains the data that will be used to label the rows in wide format. This can consist of multiple colums.
}    
  \item{column }{vector with the data that will be used to label the columns in wide format. 
}    
  \item{descriptor }{name(s) of the descriptor, i.e. *where* the data were taken, e.g. station names. Either a vector, a list, a data.frame or matrix (and with multiple columns). It can be of type numerical, character, or a factor. When a data.frame or a list the "names" will be used in the output; when a vector, the argument name will be used in the output.
}    
  \item{taxon }{gives the taxonomic name(s) (e.g. species). Should be of the same length as (the number of rows of) \code{descriptor}. Can be a list (or data.frame with one column), or a vector. When a data.frame or a list the "name" of the column will be used in the output; when a vector, the argument name will be used.
}
  \item{value }{vector or list that contains the *values* of the data, usually density. Should be of the same length as (the number of rows of) \code{descriptor} and \code{taxon}. For function \code{getDensity}, \code{value} can also be a multi-column data.frame or matrix. 
  }
  \item{averageOver }{vector with *replicates* over which averages need to be taken. Should be of the same length as (the number of rows of) \code{descriptor}.
  }
  \item{subset }{logical expression indicating elements to keep: missing values are taken as FALSE. If NULL, or absent, then all elements are used. Note that the subset is taken *after* the number of samples to average per descriptor is calculated, so this will also work for selecting certain taxa that may not be present in all replicates over which should be averaged.
  }
  \item{taxonomy }{taxonomic information; first column will be matched with \code{taxon}, regardless of its name.
  }
  \item{d.column }{position(s) or name(s) of the column(s) that holds the descriptors of the (density) data set, and that should be removed for any calculations. The default is to have the first column holding the descriptors. If \code{NULL}, or \code{0}, then there is no separate column with names, so the \code{row.names} of the dataset are used as descriptor names.
  }
  \item{t.column }{position(s) or name(s) of the column(s) that holds the taxon names of the (trait) data set, and that should be removed for any calculations. The default is to have the first column holding the taxa. If \code{NULL}, or \code{0}, then there is no separate column with names, so the \code{row.names} of the dataset are used as taxon names.
  }
  
  \item{trait }{(taxon x trait) data or (descriptor x trait) data, in WIDE format. Traits can be fuzzy coded. In the default setting, the first column contains the name of the taxon, and \code{t.column=1}. It is also allowed to have the taxa as \code{row.names} of the data.frame - set \code{t.column=0}. 
  }
  \item{w.names, taxon.names, trait.names }{names of the items constituting the columns in the wide dataset. If not given, the columnames (minus d.column) will be used. Input this as a data.frame if you want to set the names of the columns in the long format. 
  }
  \item{absences }{if \code{TRUE} the long format will contains 0's for absences
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\seealso{
\link{MWTL} for the data sets

\link{mapKey} for simple plotting functions

\link{getDensity} for functions working on density data

\link{getSummary} for estimating summaries from density data

\link{getTraitDensity} for functions combining density and traits

\link{getDbIndex} for extracting bioturbation and bioirrigation indices

\link{extendTrait} for functions working with traits

\link{getTrait} 
}

\details{

About \code{long2wide} and \code{wide2long}:

There are two ways in which density data can be inputted: 

\itemize{
 \item{descriptor, taxon, value, replicates, ...} are vectors with density data in *LONG* format: (where, which, replicates (averageOver), value); all these vectors should be of equal length (or NULL). 

 \item{wide} has the density data in *WIDE* format, i.e. as a matrix with the descriptor (and perhaps replicates) in the first column, the taxon as the column names (excluding the first column), and the content of the data is the density.
}

}
\keyword{manip}

\examples{

## ====================================================
## Long to wide format
## ====================================================

##-----------------------------------------------------
## A small dataset with replicates
##-----------------------------------------------------

# 2 stations, 2 replicates for st.a, one replicate for st.b
Bdata.rep <- data.frame(
  station   = c("st.a","st.a","st.a","st.b","st.b","st.b"),
  replicate = c(     1,     1,    2,     1,     1,     1),
  species   = c("sp.1","sp.2","sp.1","sp.3","sp.4","sp.5"),
  density   = c(     1,     2,    3,     3,     1,     3)
)
Bdata.rep

##-----------------------------------------------------
## Go to wide format, average of replicates
##-----------------------------------------------------

with (Bdata.rep, 
  l2wDensity(value      = density, 
            descriptor  = station, 
            taxon       = species, 
            averageOver = replicate))

##-----------------------------------------------------
## Go to wide format, keep replicates
##-----------------------------------------------------

with (Bdata.rep, 
  l2wDensity(value      = density, 
            descriptor = cbind(station, replicate), 
            taxon      = species))

##-----------------------------------------------------
## Go to wide format, ADD replicates
##-----------------------------------------------------

with (Bdata.rep, 
  l2wDensity(value      = density, 
            descriptor = station,  
            taxon      = species))


## ====================================================
## A small dataset without replicates
## ====================================================

Bdata <- data.frame(
  station = c("st.a","st.a","st.b","st.b","st.b","st.c"),
  species = c("sp.1","sp.2","sp.1","sp.3","sp.4","sp.5"),
  density = c(1, 2, 3, 3, 1, 3)
)

##-----------------------------------------------------
## From long to wide format
##-----------------------------------------------------

Bwide <- with (Bdata, 
  l2wDensity (value      = density, 
             descriptor = station, 
             taxon      = species))
Bwide


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

# add genus, family... to the density data

Bdata.ext <- merge(Bdata, Btaxonomy)
head(Bdata.ext)   

# estimate (summed) density on genus level 
Bwide.genus <- with (Bdata.ext, 
  l2wDensity(descriptor = station, 
             taxon      = genus,
             value      = density)
                    )

Bwide.genus

##-----------------------------------------------------
## select part of the data
##-----------------------------------------------------

# return species density for g.2 only
with (Bdata.ext, 
  l2wDensity(value     = density, 
            descriptor = station, 
            taxon      = species, 
            subset     = genus=="g.2")
      )

# create summed values for g.2 only
with (Bdata.ext, 
  l2wDensity(value     = density, 
            descriptor = station, 
            taxon      = genus, 
            subset     = genus=="g.2")
      )

## ====================================================
## From wide to long format
## ====================================================

Bwide <- data.frame(station=c("Sta", "Stb", "Stc"),
                    sp1    =c(    1,     3,     0),
                    sp2    =c(    2,     0,     0),
                    sp3    =c(    0,     0,     3))

# this long format includes the 0 densities
wide2long (wide = Bwide, absences=TRUE)

# this labels the species column appropriately
wide2long (wide = Bwide, w.name=data.frame(species=colnames(Bwide)[-1]))

}
