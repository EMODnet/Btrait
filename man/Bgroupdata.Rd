\name{Benthic typology data}
\docType{data}
\alias{Groups}
\title{
  Benthic typological groups
}

\description{
  Macrofauna main typological groups.
}

\format{
**\code{Groups}** contains 11 typological species groups representing sea floor functions, as derived from a cluster analysis on species traits by Beauchard et al. (submitted).

\itemize{
 \item{Epi3D, }{Epibenthic 3D, Builders of large and complex 3D-structures on the sediment; high biodeposition ability}.
 
 \item{Foul, }{Small fouling species of various substrata, such as barnacles and epibenthic tube builders; can create extended organic mats; high biodeposition ability}
 
 \item{Neutral, }{species deprived of the ability to alter sediments; also characteristic of biodeposition; also slowly mobile species from hard substrata}
 
 \item{ShalShel, }{Shallow shell, mostly sub-sediment surface buried bivalves; dominantly biodiffusors, shell providers; low to high biodepositors}
 
 \item{SurfDiff, }{Surficial biodiffusors, with variable degrees of burying ability; no other specific attribute}
 
 \item{SmalTub, }{Small tube dweller; mostly superficial tubicolous occurring in wide lawns of tubes or tube protrusions ; commonly biostabilisers, advectors (downward and upward conveying, moderate ventilation ability, with some degree of biodeposition}
 
 \item{DeepTub, }{Deep tube dweller; similar to smalTub, but les tubicolous and more burrow dwellers (I- or J-shaped), and deeper burrowers; higher ventilation ability}
 
 \item{MinBiot, }{Minor bioturbator, larger and less sedentary, deeper burrowing ability than DeepTub; mostly advectors}
 
 \item{SesBio, }{Sessile bioturbator, deep burrowing ability, great bioirrigation potential due to wider burrow}
 
 \item{MajBiot, }{Major bioturbator, deep burrowing, mobile and non-shelly; diffusive mixing}
 
 \item{Deep3D, }{Deepest burrowers, long term burrows, very high level of sediment mixing and irrigation}
}

Explanation on its contents is given in its attribute \code{description} (use \code{metadata(Groups)} to derive it)
}

\details{
See Beauchard et al., subm.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}

\examples{

##-----------------------------------------------------
## Show content
##-----------------------------------------------------

NSt <- metadata(Groups)
NSt

##-----------------------------------------------------
## 10 most commonly found organisms per group
##-----------------------------------------------------

MWTL.groups <- merge(MWTL$density, Groups)

for (i in 1:nrow(NSt)){

  NSdta <- subset(MWTL.groups, 
                  subset= (typology == NSt[i,1]))
  
  st    <- sort(table(NSdta$taxon), decreasing=TRUE)
  
  cat(NSt[i,2], " : ", 
     paste(names(st)[1:min(10, length(st))], 
           deparse=", ", sep=""), "\n\n")
}

##-----------------------------------------------------
## Distribution of deep 3D burrowers
##-----------------------------------------------------

# Tsese formulations are all the same:

\dontrun{

# explicit input of data - need to specify the name of the row 
Data <- long2wide(row         = list(station = MWTL.groups$station), 
                  column      = MWTL.groups$typology, 
                  value       = MWTL.groups$density,
                  averageOver = MWTL.groups$year, 
                  subset      = (MWTL.groups$typology=="Deep3D"))
head(Data)

# use of with() to create environment
Data <- with (MWTL.groups, 
              long2wide(row         = list(station = station), 
                        column      = typology, 
                        value       = density,
                        averageOver = year, 
                        subset      = (typology=="Deep3D")))
head(Data)

# use of argument data to create environment
Data <- long2wide(data        = MWTL.groups, 
                  row         = station, 
                  column      = typology, 
                  value       = density,
                  averageOver = year, 
                  subset      = (typology=="Deep3D"))
head(Data)

# piping
Data <- MWTL.groups |> long2wide( 
                  row         = station, 
                  column      = typology, 
                  value       = density,
                  averageOver = year, 
                  subset      = (typology=="Deep3D"))
head(Data)

}


Data <- with (MWTL.groups, 
            long2wide(row         = station, 
                      column      = typology, 
                      value       = density,
                      averageOver = year, 
                      subset      = (typology=="Deep3D")))
head(Data)
Data.st <- merge(MWTL$stations, Data, 
                 by.x = "station", by.y = "descriptor")

with (Data.st, map_key(x = x, y = y, colvar = Deep3D, 
                       main = "deep 3D", 
                       clab = "ind/m2", log = "c"))


}

\references{
The species typology is described in the following manuscript:

Beauchard, O., Murray S.A. Thompson, Kari Elsa Ellingsen, Gerjan Piet, Pascal Laffargue, Karline Soetaert, 2023. Assessing sea floor functional diversity and vulnerability. Marine Ecology Progress Serie v708, p21-43, https://www.int-res.com/abstracts/meps/v708/p21-43
}

\seealso{
  \link{Traits_nioz} for the benthic trait datasets.

  \link{MWTL} for other data sets of the Northsea.

  \link{map_key} for plotting.

  \link{extend_trait} for functions operating on trait data.

  \link{get_density} for functions operating on density data.

  \link{get_trait_density} for functions operating on density and trait data.

  \link{get_Db_index} for extracting bioturbation and bioirrigation indices.
}

\keyword{datasets}
