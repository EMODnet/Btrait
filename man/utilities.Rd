\name{General utilities for the Btrait package}
\alias{get_taxon_info}
\alias{get_taxonomy}
\alias{metadata}
\title{

  \code{get_taxon_info} gets trait information for a taxon
  
  \code{get_taxonomy} gets taxonomic information for a taxon
  
  \code{metadata} gets the description of the content of a database.
}

\description{
 \code{get_taxon_info} extracts trait information for a certain taxon from a trait database.

 \code{metadata} returns a description of an object, i.e the meaning of the columns, units, etc...
}

\usage{
get_taxon_info(taxon, trait = Traits_nioz, taxonomy = Taxonomy)  

get_taxonomy(taxon, taxonomy=Taxonomy)

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
 \code{get_taxon_info} returns a data.frame with the traits of the species whose modalities are not Null, NULL or None. 
 
 \code{get_taxonomy} returns the occurrences of \code{taxon} within \code{taxonomy}. 

 \code{metadata} returns a data.frame with a description of the contents of the database.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}
\seealso{.
\link{MWTL} for the data sets.

\link{map_key} for simple plotting functions.

\link{get_density} for functions working with density data.

\link{get_Db_index} for extracting bioturbation and bioirrigation indices.

\link{extend_trait} for functions working with traits.

\link{get_trait} 

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
get_taxonomy("Amphiura filiformis")
get_taxonomy("Amphiura")

get_taxon_info("Amphiura filiformis")
get_taxon_info("Amphiura filiformis", 
               trait = Traits_cefas)

# there is no direct information on the traits of this species in Traits_cefas
get_taxon_info("Amphiura filiformis", 
               trait = Traits_cefas, NULL)

## ====================================================
## Northsea dataset
## ====================================================

# show its contents
metadata(MWTL)          
metadata(MWTL$density)

}
