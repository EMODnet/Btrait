\name{Standardizing data sets}
\alias{standardize_taxonomy}
\alias{standardize_trait}
\title{
 standarization of taxonomic and trait datasets
}

\description{
  \code{standardize_taxonomy} checks taxonomic datasets - fills in blands, checks for duplicates
  
  \code{standardize_trait} checks trait data (in wide format) so that modalities for one trait sum to 1
}

\usage{
standardize_trait(trait, taxon_column = 1, 
                  trait_class, trait_modality,
                  meta = data.frame(trait_class    = trait_class, 
                                    trait_modality = trait_modality))
standardize_taxonomy(taxonomy, verbose = FALSE)

}

\arguments{
  \item{trait }{(taxon x trait) data or (descriptor x trait) data, in *WIDE* format, and containing numerical values only. Traits can be fuzzy coded. The number of columns with taxonomic information is specified with \code{taxon_column}. The default is to have the first column contain the name of the taxon, and \code{taxon_column=1}. It is also allowed to have the taxa as \code{row.names} of the data.frame; in this case, \code{taxon_column=0}.
  }  
  \item{taxon_column }{position(s) or name(s) of the column(s) that holds the taxon names of the (trait) data set, and that should be removed for any calculations. The default is to have the first column holding the taxa. If \code{NULL}, or \code{0}, then there is no separate column with names, so the \code{row.names} of the dataset are used as taxon names.
  }
  \item{trait_class }{indices to trait levels, a vector. The length of this vector should equal the number of columns of \code{trait} or \code{wide} minus the value of \code{taxon_column}. The order should be conform the order in the trait columns (not checked). 
  }
  \item{trait_modality }{names of the modalities of the trait, a vector. Should be of same length as \code{trait_class}; it will be matched to the column names of \code{trait}
  }
  \item{meta }{The trait metadata, containing the trait class and modalities.
  }
  \item{taxonomy }{taxonomic information (the relationships between the taxa), a data.frame; first column will be matched with \code{taxon} (in the data.frames \code{data} and \code{trait}). The subsequent columns should have increasing taxonomic level. (e.g. the column order should be *species, genus, family, ...*. The taxonomic relations from \code{data} and \code{trait} are used to estimate traits of taxa that are not accounted for, and that will be estimated based on taxa at the nearest taxonomic level. See details.
  }
  \item{verbose }{when TRUE, will write warnings to the screen.
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
## standardizing a trait database
## ====================================================

Btraits <- data.frame(
  taxon   = c("sp.1","sp.2","sp.3","sp.5","sp.6"),
  T1_M1   = c(0     , 0    ,   0  , 1  ,     0.5),
  T1_M2   = c(1     , 0    ,   1  , 2  ,     0),
  T1_M3   = c(0     , 1    ,   1  , 3  ,     0),
  T2_M1   = c(0     , 0    ,   1  , 1  ,     1),
  T2_M2   = c(1     , 1    ,   0  , 1  ,     0)
)

# The metadata for this trait
Meta <- data.frame(
  trait    = c("T1"   , "T1"   , "T1"   , "T2"   , "T2"),
  modality = c("T1_M1", "T1_M2", "T1_M3", "T2_M1", "T2_M2")
)

Btraits
Meta

standardize_trait(Btraits, 
                  trait_class    = Meta$trait,
                  trait_modality = Meta$modality)

# This will give an error

\dontrun{

# The metadata lack one modality
Meta2 <- data.frame(
  trait    = c("T2"   , "T2",    "T1"   , "T1"       ),
  modality = c("T2_M1", "T2_M2", "T1_M1", "T1_M2" )
)

# gives an error
standardize_trait(Btraits, 
                  trait_class    = Meta2$trait,
                  trait_modality = Meta2$modality)


# The traits for species sp1 and sp 2 do not have modality for trait 1:

standardize_trait(Btraits[,-ncol(Btraits)], # remove last column
                  trait_class    = Meta$trait,
                  trait_modality = Meta$modality)

}


## ====================================================
## Standardizing taxonomy
## ====================================================


Btaxonomy <- data.frame(
  species = c("sp.1","sp.2","sp.3","sp.4","sp.5","sp.6"),
  genus   = c( "g.1", " ", "g.2", "g.2", "g.3", "g.4"),
  family  = c( "f.1", " ", "f.1", "f.1", "f.2", "f.3"),
  order   = c( "o.1", "o.1", "o.1", "", "o.2", "o.2"),
  class   = c( "c.1", "c.1", "c.1", "c.1", "c.1", "c.1")
  )

Btaxonomy
standardize_taxonomy(Btaxonomy)
}
