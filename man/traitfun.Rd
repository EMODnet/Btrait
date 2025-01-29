\name{Trait and taxonomic functions}
\alias{get_trait}
\alias{extend_trait}
\alias{fuzzy2crisp}
\alias{crisp2fuzzy}
\title{
  Functions for working with trait and taxonomic data.
}

\description{
 \code{get_trait} gets trait information for (unique) taxa.

 \code{fuzzy2crisp} translate fuzzy traits into crisp traits.

 \code{crisp2fuzzy} translate crisp traits into fuzzy traits. Note: this makes sense only for categorical traits.

 \code{extend_trait} expands taxon x trait data based on their taxonomic composition.
}

\usage{
get_trait(taxon, trait, trait.class=NULL, trait.score=NULL, 
          taxonomy = NULL, t.column = 1, 
          standardize = FALSE, verbose = FALSE)

extend_trait(trait, taxonomy, t.column = 1) 

fuzzy2crisp(trait, trait.class, trait.score, t.column = 1, 
        standardize = TRUE)

crisp2fuzzy (trait, t.column = 1)
}

\arguments{
  \item{taxon }{vector with the taxon names for which to get the trait information. 
 }
  \item{trait}{ (taxon x trait) data or (descriptor x trait) data, in WIDE format. 
  Traits can be fuzzy coded. In the default setting, the first column contains the name of the taxon, and \code{t.column=1}. It is also allowed to have the taxa as \code{row.names} of the data.frame - set \code{t.column=0}. 
  }
  \item{trait.class }{indices to trait classes, a vector. The length of this vector should equal the number of columns of \code{trait} minus the descriptor columns  (nrow(\code{trait}) - \code{t.column}). If not \code{NULL}, the fuzzy traits will be converted to crisp values.
  }
  \item{trait.score }{trait values or scores, a vector. Should be of same length as \code{trait.class}.
  }
  \item{standardize }{when TRUE, will standardize the trait modalities, so that for each trait class they sum to 1. Only used when \code{trait.class} is given a value.
  }
  \item{verbose }{when TRUE, will write warnings to the screen.
  }
  \item{taxonomy }{taxonomic information (the relationships between the taxa), a data.frame; first column will be matched with \code{taxon}, regardless of its name. This is used to estimate traits of taxa that are not accounted for, and that will be estimated based on taxa at the nearest taxonomic level. See details.
  }
  \item{t.column }{position(s) or name(s) of the column(s) that holds the taxon of the data set. The default is to have the first column holding the taxon. If \code{NULL}, or \code{0}, then there is no separate column with names, so the \code{row.names} of the dataset are used as the taxon names.
  }
}
\value{

\code{extend_trait} returns the trait matrix, augmented with higher taxonomic levels,

\code{get_trait} returns the taxon x trait matrix, in fuzzy format (if the \code{trait} database was in fuzzy format). However, when \code{trait.class} is not NULL, the matrix will be in crisp format.

\code{fuzzy2crisp} returns the trait matrix, in crisp format (not fuzzy coded),

\code{crisp2fuzzy} returns the trait matrix, in fuzzy format. If the original,fuzzy trait matrix was categorical, the output will include a description of each column, in its attribute \code{description}, use \code{metadata(...)} to extract this.

}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}
\seealso{
\link{MWTL} for the data sets.

\link{map_key} for simple plotting functions.

\link{get_density} for functions working with density data.

\link{get_summary} for estimating summaries from density data.

\link{get_trait_density} for functions combining density and traits.

\link{get_Db_index} for extracting bioturbation and bioirrigation indices.
}
\details{
The \code{taxonomy} is used to fill in the gaps of the trait information, assuming that closely related taxa will share similar traits. This is done in several steps:

In function \code{extend_trait} the traits are extended with information on higher taxonomic levels, provided that information is not yet in the trait database. The traits for a taxonomic level are estimated as the average of the traits at the lower level. For instance, traits on genus level will be averages of known traits of all species in the database belonging to this genus. 

In function *get_trait*, the trait database is first extended with information on higher taxonomic levels (using *extend_trait*). Then, for each taxon that is not present in the trait database, the traits on the closest taxonomic level are used. For instance, for an  unrecorded species, it is first checked if the trait is known on genus level, if not, family level and so on.


for \code{get_trait}, only the traits for unique taxa are returned.
Traits will be \code{NA} if traits were not present in the trait database and that 
could not be derived based on taxonomic closeness. The list of taxa whose 
traits remain unknown can be found in \code{attributes(..)$notrait}.

It is best, for subsequent calculations, to remove the \code{NA}s in the result, 
or to put them to 0 (when the remaining unassigned taxa are deemed to be acceptable).
}
\keyword{manip}

\examples{

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

## ====================================================
## Small dataset: fuzzy-coded traits
## ====================================================

Btraits <- data.frame(
  species = c("sp.1","sp.2","sp.3","sp.5","sp.6"),
  T1_M1   = c(0     , 0    ,   0  , 0.2  ,     1),
  T1_M2   = c(1     , 0    , 0.5  , 0.3  ,     0),
  T1_M3   = c(0     , 1    , 0.5  , 0.5  ,     0),
  T2_M1   = c(0     , 0    ,   1  , 0.5  ,     1),
  T2_M2   = c(1     , 1    ,   0  , 0.5  ,     0)
)

# The 'metadata' of this trait database
Btraits.lab <- data.frame(
  colname  =c("T1_M1","T1_M2","T1_M3","T2_M1","T2_M2"),
  trait    =c("T1"   ,"T1"   ,"T1"   ,"T2"   ,"T2"),
  modality =c("M1"   ,"M2"   ,"M3"   ,"M1"   ,"M2"), 
  score    =c(0      , 0.5   , 1     , 0.2   , 2)
)

##-----------------------------------------------------
## Show traits
##-----------------------------------------------------

# sp.4 is not in Btraits 
get_trait (taxon      = "sp.4", 
          trait      = Btraits)

# sp.4 traits derived from taxonomic tree (g.2)
get_trait (taxon      = "sp.4", 
          trait      = Btraits, 
          taxonomy   = Btaxonomy)

get_trait (taxon      = c("g.2"), 
          trait      = Btraits, 
          taxonomy   = Btaxonomy)

# g.2 is derived as mean of sp.2, sp.3 and sp.4 
get_trait (taxon      = c("sp.2", "sp.3", "sp.4"), 
          trait      = Btraits, 
          taxonomy   = Btaxonomy)

get_trait (taxon      = c("g.2"), 
          trait      = Btraits, 
          taxonomy   = Btaxonomy)

##-----------------------------------------------------
## crisp values of traits
##-----------------------------------------------------

# for categoric values: the most abundant category is taken
 C1 <- fuzzy2crisp(trait       = Btraits, 
                   trait.class = Btraits.lab$trait, 
                   trait.score = Btraits.lab$modality, 
                   standardize = TRUE)
 C1
 
# the reverse returns a binary-coded value
 (C2fuz <- crisp2fuzzy(trait       = C1))
 metadata(C2fuz)
 
# for numeric (or binary) values: the weighted mean is calculated
 C2 <- fuzzy2crisp(trait       = Btraits, 
                   trait.class = Btraits.lab$trait, 
                   trait.score = Btraits.lab$score, 
                   standardize = TRUE)
 C2
 (C2f <- crisp2fuzzy(C2))  # this has no effect
 
##-----------------------------------------------------
## Extend traits with higher level information
##-----------------------------------------------------

Btraits.ext <- extend_trait(trait    = Btraits, 
                            taxonomy = Btaxonomy)

Btraits.all <- rbind(Btraits, Btraits.ext)
Btraits.all

#same, but in crisp format
fuzzy2crisp(trait       = Btraits.all, 
            trait.score = Btraits.lab$score, 
            trait.class = Btraits.lab$trait)

# In one go, and including sp.4
get_trait(taxon    = c("sp.1", "sp.2", "sp.3", "sp.4", "sp.5", "sp.6"),
          trait    = Btraits, 
          taxonomy = Btaxonomy, 
          trait.score = Btraits.lab$score, 
          trait.class = Btraits.lab$trait)
          
          
TraitSmall <- data.frame(
                species = c("Alcyonium acaule",  "Alcyonium coralloides"),
                T1.M1   = c(0.9, 1), 
                T1.M2   = c(0.1, 0), 
                T1.M3   = c(0,   0))
TaxSmall <- data.frame(
              species = c("Alcyonium acaule",  "Alcyonium coralloides"), 
              Genus   = c("Alcyonium", "Alcyonium"), 
              Family  = c("Alcyoniidae", "Alcyoniidae"), 
              order   = c("Alcyonacea", "Alcyonacea"))

T1 <- extend_trait(trait    = TraitSmall, 
                   taxonomy = TaxSmall)         
}
