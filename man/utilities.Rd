\name{General utilities for the Btrait package}
\alias{getTaxonInfo}
\alias{getTaxonomy}
\alias{metadata}
\title{

  getTaxonInfo gets trait information for a taxon
  
  getTaxonomy gets taxonomic information for a taxon
  
  metadata gets the description of the content of a database.
}

\description{
 \code{getTaxonInfo} extracts trait information for a certain taxon from a trait database.

 \code{metadata} returns a description of an object, i.e the meaning of the columns, units, etc...
}

\usage{
getTaxonInfo(taxon, trait = Traits_nioz, taxonomy = Taxonomy)  

getTaxonomy(taxon, taxonomy=Taxonomy)

metadata(object)        
}

\arguments{
  \item{taxon }{name of the taxon (or vector with names) for which the traits need to be listed.
}
  \item{trait }{a trait database in wide format; it should have taxon as first column, and it should contain a description of the trait database in its attribute ("description"). See \link{Traits_nioz} for a suitable structure. 
  }
  \item{taxonomy }{taxonomic information; first column will be matched with \code{taxon} from the trait database, regardless of its name.
  }
  \item{object }{any object that has a data.frame \code{description} in its \code{attributes}.
  }
}

\value{
 \code{getTaxonInfo} returns a data.frame with the traits of the species whose modalities are not Null, NULL or None. 
 
 \code{getTaxonomy} returns the occurrences of \code{taxon} within \code{taxonomy}. 

 \code{metadata} returns a data.frame with a description of the contents of the database.
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
## Description of the Traits_nioz database
## ====================================================
metadata(Traits_nioz)

## ====================================================
## Traits for certain taxa
## ====================================================

# Taxonomy of amphiura and related species in the databases.
getTaxonomy("Amphiura filiformis")
getTaxonomy("Amphiura")

getTaxonInfo("Amphiura filiformis")
getTaxonInfo("Amphiura filiformis", trait=Traits_cefas)

# there is no direct information on the traits of this species in Traits_cefas
getTaxonInfo("Amphiura filiformis", trait=Traits_cefas, NULL)

## ====================================================
## Northsea dataset
## ====================================================

# show its contents
metadata(MWTL)          
metadata(MWTL$density)

}
