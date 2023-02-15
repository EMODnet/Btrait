\name{Bioturbation and Bioirrigation potential indices}
\alias{getDbIndex}
\alias{getIrrIndex}
\title{
  Functions for obtaining the bioturbation and bio-irrigation potential index.
}

\description{
 \code{getDbIndex} estimates the BPc index sensu Querios et al., 2013 and Solan et al., 2004.

 \code{getIrrIndex} estimates the IPc index sensu Wrede et al., 2018.
}

\usage{
getDbIndex (descriptor, taxon, density, biomass, averageOver = NULL,    
       subset, trait = Traits_Db, taxonomy = NULL, 
       full.output = FALSE, verbose = FALSE)

getIrrIndex (descriptor, taxon, density, biomass, averageOver = NULL,    
       subset, trait = Traits_irr, taxonomy = NULL, 
       full.output = FALSE, verbose = FALSE)
}

\arguments{
  \item{descriptor }{name(s) of the descriptor, i.e. *where* the data were taken, e.g. station names. Either a vector, a list, a data.frame or matrix (and with multiple columns). It can be of type numerical, character, or a factor. When a data.frame or a list the "names" will be used in the output; when a vector, the argument name will be used in the output.
}    
  \item{taxon }{vector describing *what* the data are; it gives the taxonomic name (e.g. species). Should be of the same length as (the number of rows of) \code{descriptor}. Can be a list (or data.frame with one column), or a vector. When a data.frame or a list the "name" will be used in the output; when a vector, the argument name will be used.
}
  \item{density, biomass }{vector or list that contains the *density* and *biomass* of the data. Should be of the same length as (the number of rows of) \code{descriptor} and \code{taxon}. 
  }
  \item{averageOver }{vector with *replicates* over which averages need to be taken. Should be of the same length as (the number of rows of) \code{descriptor}.
  }
  \item{subset }{logical expression indicating elements to keep from the density data: missing values are taken as FALSE. If NULL, or absent, then all elements are used. Note that the subset is taken *after* the number of samples to average per descriptor is calculated, so this will also work for selecting certain taxa that may not be present in all replicates over which should be averaged.
  }
  \item{trait }{(taxon x trait) data, in *WIDE* format, and containing numerical values only. The first column should contain the name of the taxa. For function \code{getDbIndex}, also the columns named \code{Mi} and \code{Ri}, denoting the mobility and reworking mode (values between 1-4 and 1-5 respectively) should be present. For function \code{getIrrIndex}, the columns named \code{BT} (burrowtype, 1-3), \code{FT} (feeding type, 1-3), and \code{ID} (injection depth, 1-4) should be present. Good choices are \link{Traits_Db} and \link{Traits_irr} (the defaults).
  }
  \item{verbose }{when TRUE, will write warnings to the screen.
  }
  \item{full.output }{when TRUE, will output the full data.frame with the descriptor x taxon indices (called \code{all}. See \code{value}.
  }
  \item{taxonomy }{taxonomic information (the relationships between the taxa), a data.frame; first column will be matched with \code{taxon}, regardless of its name. This is used to estimate traits of taxa that are not accounted for, and that will be estimated based on taxa at the nearest taxonomic level. See details.
  }
}

\details{
 The algorithm first calls function \code{getDensity}, to obtain the (depending on \code{averageOver} averaged) taxon densities and biomass per descriptor. The weight is estimated from biomass and density.

Then, for each taxon in the obtained dataset, the required traits are extracted from the trait database using function \code{getTrait}. 

The two data.frames are then merged (based on taxon), so that for each descriptor x taxon occurrence, the density, weight and required traits.    

The bioturbation or bioirrigation Index is then estimated by using the appropriate formula,

Finally, the sums of the taxon indices are taken per descriptor, and the averages are estimated for the taxa, using the R-function \code{aggregate}.
}

\value{
Both return a \code{list} with the following elements:

 \itemize{
 \item{descriptor} a data.frame with two columns, the descriptor, and the index (BPc or IPc), 
    which consist of the *summed* values over all taxa
 \item{taxon} a data.frame with two columns, the taxon name, and the index (BPc or IPc), 
    which is *averaged* over all descriptors
 \item{all}, will only be present if \code{full.output} is \code{TRUE}: the full dataset on which the indices were estimated.    
 }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\seealso{
\link{MWTL} for the data sets

\link{mapBtrait} for simple plotting functions

\link{getDensity} for functions working with density data

\link{getTraitDensity} for functions operating on density and trait data.

\link{extendTrait} for functions working with traits

\link{getTrait} 
}

\references{
Queiros, Ana M., Silvana N. R. Birchenough, Julie Bremner, Jasmin A. Godbold, Ruth E. Parker,  Alicia Romero-Ramirez, Henning Reiss, Martin Solan, Paul J. Somerfield, Carl Van Colen,  Gert Van Hoey, Stephen Widdicombe, 2013.  A bioturbation classification of European marine infaunal invertebrates.  Ecology and Evolution 3 (11), 3958-3985

Solan M, Cardinale BJ, Downing AL, Engelhardt KAM, Ruesink JL, Srivastava DS. 2004. Extinction and ecosystem function in the marine benthos. Science 306:1177-80.

A. Wrede, J.Beermann, J.Dannheim, L.Gutow, T.Brey, 2018. Organism functional traits and ecosystem supporting services - A novel approach to predict bioirrigation. Ecological indicators, 91, 737-743.

}

\note{
Equations:
\itemize{
\item
The formula for estimating the bioturbation Index for taxon i (as in Querios et al., 2013)is:

BPc_i = sqrt(Wi) * density_i * Ri*Mi,

where \code{Mi} and \code{Ri}, denote the mobility and reworking mode (values between 1-4 and 1-5 respectively). See \link{Traits_Db} for what these numbers mean.
 
\item
The formula for estimating the bioirrigation Index for taxon i (as in Wrede et al., 2018) is:

IPc_i = (Wi)^(0.75) * density_i * BTi*FTi*IDi,

where \code{BT} is burrowtype (1-3), \code{FT}  is feeding type (1-3), and \code{ID} is injection depth (BF1-4).
See \link{Traits_irr} for what these numbers mean.
}

The stations Index is the sum of all species indices.
}
\keyword{manip}

\examples{

##-----------------------------------------------------
## The bioturbation potential for one species
##-----------------------------------------------------

# Amphiura filiformis, for increasing density 
DbAmp <- getDbIndex(
              taxon   = rep("Amphiura filiformis", times=10), 
              density = 1:10, 
              biomass = (1:10)*4.5e-3, 
              full.output=TRUE, 
              trait   = Traits_Db)
head(DbAmp$all)

##-----------------------------------------------------
## The bioirrigation potential for one species
##-----------------------------------------------------

# Amphiura filiformis, in dutch part of the northsea 
IrrAmp <- with (MWTL$density, 
     getIrrIndex(
              descriptor  = station,         
              taxon       = taxon, 
              subset      = taxon == "Amphiura filiformis", 
              averageOver = year,
              density     = density, 
              biomass     = biomass, 
              full.output = TRUE, 
              trait       = Traits_irr))

# irrigation activity per station
head(IrrAmp$descriptor) 

# add coordinates of the stations
IRR.amphiura <- merge(MWTL$stations, IrrAmp$descriptor, by=1)

# create a map
with(IRR.amphiura, mapBtrait(x=x, y=y, colvar=IPc, pch=18,
                     contours=MWTL$contours, draw.levels=TRUE))
                     
IrrAmp$taxon  # average irrigation activity

##-----------------------------------------------------
## The bioturbation potential index for communities
##-----------------------------------------------------
# BPc = sqrt(weight) * density * Mi * Ri

# BPc for the Dutch part of the Northsea, in 1995
BPC <- with (MWTL$density,
      getDbIndex(descriptor = list(station=station),
                     subset     = (year == 1995),
                     taxon      = taxon, 
                     density    = density, 
                     biomass    = biomass, 
                     trait      = Traits_Db,      
                     taxonomy   = Taxonomy))  

# There is one taxon for which trait could not be derived
attributes(BPC)$notrait  

# Total BPC per station
head(BPC$descriptor)

# Average BPC per taxon (only where taxon is present)
head(BPC$taxon)

# Plot the results (after adding the coordinates)
BPC_MWTL <- merge(MWTL$stations, BPC$descriptor, by=1)

with (BPC_MWTL, mapBtrait(x=x, y=y, colvar=BPc, pch=18, 
       contours=MWTL$contours, draw.levels=TRUE, main= "BPc"))

# The 10 main bioturbators in the data set:
head(BPC$taxon[order(BPC$taxon$BPc, decreasing=TRUE), ], n=10)

##-----------------------------------------------------
## The bioirrigation Index
##-----------------------------------------------------
# IPc = (weight_i)^0.75 * density * FTi * BTi * IDi

# IPc for the NSBS station OESTGDN19 in the OysterGrounds
IPC <- with (MWTL$density,
      getIrrIndex(descriptor  = list(station=station),
                  subset      = (station == "OESTGDN19"),
                  taxon       = taxon, 
                  averageOver = year,
                  density     = density, 
                  biomass     = biomass, 
                  trait       = Traits_irr,      
                  taxonomy    = Taxonomy,
                  full.output = TRUE))  

# The 10 main bioirrigators in the data set, and why :
head(IPC$all[order(IPC$all$IPc, decreasing=TRUE), ], n=10)

}
