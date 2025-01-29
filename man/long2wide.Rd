\name{Conversion functions}
\alias{long2wide}
\alias{wide2long}
\alias{l2w_density}
\alias{l2w_trait}
\alias{w2l_density}
\alias{w2l_trait}
\title{
  Functions for conversion from long to wide format and vice versa.
}

\description{
 \code{long2wide} casts data from long to wide format.

 \code{w2l_density} casts density data from wide to long format.

 \code{w2l_trait} casts trait data from wide to long format.

 \code{wide2long} casts data from wide to long format.

 \code{l2w_density} casts density data from long to wide format.

 \code{l2w_trait} casts trait data from long to wide format.

}

\usage{
wide2long(wide, d.column = 1, w.names = NULL, 
          absences = FALSE)

w2l_density(wide, d.column = 1, taxon.names = NULL, 
           absences = FALSE) 

w2l_trait(wide, t.column = 1, trait.names = NULL, 
           absences = FALSE) 

long2wide(data, row, column, value, averageOver, 
        taxonomy = NULL, subset)
        
l2w_density(data, descriptor, taxon, value, averageOver, 
        taxonomy = NULL, subset)
        
l2w_trait(trait, descriptor, taxon, value, averageOver, 
         taxonomy = NULL, subset)
}

\arguments{

  \item{wide }{data, in *WIDE* format. For density data, this is a data.frame or matrix with (descriptor x taxon) information, and the first column usually contains the descriptor name. For trait data this is a data.frame with (taxon x trait) information, and the first column generally contains the names of the taxa. It is also allowed to have the descriptors as row.names of the data.frame -this requires setting \code{d.column=0}.   
  }
  \item{row }{vector or data.frame that contains the data that will be used to label the rows in wide format. This can consist of multiple colums.
}    
  \item{column }{vector with the data that will be used to label the columns in wide format. 
}    
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

\link{map_key} for simple plotting functions

\link{get_density} for functions working on density data

\link{get_summary} for estimating summaries from density data

\link{get_trait_density} for functions combining density and traits

\link{get_Db_index} for extracting bioturbation and bioirrigation indices

\link{extend_trait} for functions working with traits

\link{get_trait} 
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
## Datasets
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
## A small dataset without replicates
##-----------------------------------------------------

Bdata <- data.frame(
  station = c("st.a","st.a","st.b","st.b","st.b","st.c"),
  species = c("sp.1","sp.2","sp.1","sp.3","sp.4","sp.5"),
  density = c(1, 2, 3, 3, 1, 3)
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

## ====================================================
## Long to wide format
## ====================================================

##-----------------------------------------------------
## Go to wide format, average of replicates
##-----------------------------------------------------

# use with() to create an environment -> first column called "descriptor"
  with (Bdata.rep, l2w_density(
              value       = density, 
              descriptor  = station, 
              taxon       = species, 
              averageOver = replicate))

# use data argument -> first column called "station"
  l2w_density(data        = Bdata.rep, 
              value       = density, 
              descriptor  = station, 
              taxon       = species, 
              averageOver = replicate)

##-----------------------------------------------------
## Go to wide format, keep replicates
##-----------------------------------------------------

  l2w_density(data      = Bdata.rep,
              value      = density, 
              descriptor = cbind(station, replicate), 
              taxon      = species)

##-----------------------------------------------------
## Go to wide format, ADD replicates
##-----------------------------------------------------

  l2w_density(data       = Bdata.rep,
              value      = density, 
              descriptor = station,  
              taxon      = species)


##-----------------------------------------------------
## Go to wide format, AVERAGE over replicates
##-----------------------------------------------------

  l2w_density(data        = Bdata.rep,
              value       = density, 
              descriptor  = station,  
              averageOver = replicate,
              taxon       = species)


##-----------------------------------------------------
## density on higher taxonomic level
##-----------------------------------------------------

# add genus, family... to the density data

Bdata.ext <- merge(Bdata, Btaxonomy,
                   by = "species")
head(Bdata.ext)   

# estimate (summed) density on genus level 
Bwide.genus <- l2w_density(
              data       = Bdata.ext, 
              descriptor = station, 
              taxon      = genus,
              value      = density)

Bwide.genus

##-----------------------------------------------------
## select part of the data
##-----------------------------------------------------

# return species density for g.2 only
  l2w_density(data       = Bdata.ext,  
              value      = density, 
              descriptor = station, 
              taxon      = species, 
              subset     = Bdata.ext$genus=="g.2")
    

# create summed values for g.2 only
  l2w_density(data       = Bdata.ext,  
              value      = density, 
              descriptor = station, 
              taxon      = genus, 
              subset     = Bdata.ext$genus=="g.2")

## ====================================================
## From wide to long format
## ====================================================

  Bwide <- data.frame(station = c("Sta", "Stb", "Stc"),
                      sp1     = c(    1,     3,     0),
                      sp2     = c(    2,     0,     0),
                      sp3     = c(    0,     0,     3))

# this long format includes the 0 densities
  wide2long (wide     = Bwide, 
             absences = TRUE)

# this does not include the absences, and renames the species
  wide2long (wide = Bwide, 
             w.name = paste("Species", 1:3, sep="_"))

## ====================================================
## From wide trait data to long format
## ====================================================

head(Traits_nioz, n = c(3, 5))

T_long   <- w2l_trait(Traits_nioz)
T_long_a <- w2l_trait(Traits_nioz, absences = TRUE)

head(T_long)
head(T_long_a)

## ====================================================
## From long trait data to wide format
## ====================================================

head(T_long)

# go back from long to wide format
T_wide   <- l2w_trait(trait      = T_long,
                      taxon      = taxon, 
                      descriptor = name, 
                      value      = value)
head(T_wide, n=c(5,5))

# other way around
T_wide_2  <- long2wide(data      = T_long,
                      row        = name, 
                      column     = taxon, 
                      value      = value)

head(T_wide_2, n=c(5,5))
}
