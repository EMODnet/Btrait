\name{Benthic trait data}
\docType{data}
\alias{Traits_nioz}
\alias{Traits_cefas}
\alias{Traits_Db}
\alias{Traits_irr}
\alias{Taxonomy}
\title{
  Benthic effect and response trait characteristics, and taxonomy.
}

\description{
  Macrofauna species trait information datasets.
  
  \itemize{
   \item{Traits_nioz}: 32 effect and response traits for benthic macrofauna (Beauchard et al, 2021)
   \item{Traits_cefas}: 10 effect and response traist for benthic macrofauna (Clare et al, 2022) 
   \item{Traits_Db}: traits to estimate bioturbation potential index (Queiros et al., 2013)
   \item{Traits_irr}: traits to estimate bioirrigation potential index as in (Wrede et al., 2018)
   \item{Taxonomy}: taxonomic information for all taxa in the trait and density datasets
  }
}

\details{

 \itemize{
 \item{**\code{Traits_nioz}**} contains 15 \code{effect} and 15 \code{response} traits of 281 macrobenthos taxa (as in Beauchard et al., 2021), mainly on species level. Traits are fuzzy coded, i.e. modalities can range inbetween 0 and 1, and a taxon can be attributed a score for several modalities. Within a trait, and for a taxon, the scores sum to 1.  

 \item{**\code{Traits_cefas}**} contains 10 traits for 1025 macrobenthos taxa (as in Beauchard et al., 2021), mainly on genus level. Traits are fuzzy coded.  

 \item{**\code{Traits_Db}**} has the traits for 1094 macrobenthos taxa (as in Queiros et al., 2013), mainly on species level. The original 1032 traits were amended with additional species (as from Ilvo - courtesy Gert van Hoey).
 
  \item{**\code{Traits_irr}**} has the traits for 281 macrobenthos taxa, necessary to estimate the irrigation potential index as in Wrede et al., 2018. This information has been derived from the Traits_nioz dataset (see last example).

 \item{**\code{Taxonomy}**} has the taxonomic information for all taxa, as derived from WoRMS. 
 }
 
}

\format{
Dataset \code{Taxonomy} is a \code{data.frame} with column names \code{taxon}, \code{genus}, \code{family}, \code{order}, \code{class}, \code{phylum} and \code{AphiaID}.

Dataset \code{Traits_nioz} is a \code{data.frame} with column names \code{taxon}, \code{ET1.M1}, \code{ET1.M2}, ... \code{RT10.M3}, \code{RT10.M4}; the latter are the trait modalities, e.g. \code{RT10.M4} is modality 4 for response trait 10.

Dataset \code{Traits_cefas} is a \code{data.frame} with column names \code{taxon}, \code{sr_Less_than_10}, \code{sr_11_to_20}, \code{sr_21_to_100}, \code{sr_101_to_200}, ... \code{b_Downward_conveyer}, \code{b_None}.

Dataset \code{Traits_Db} is a \code{data.frame} with column names \code{taxon}, \code{Ri}, \code{Mi} \code{Fti} and \code{source}. 

All trait datasets have a data.frame called \code{description} in their attributes that explains the meaning of the traits, their modalities and units, and a value and score. 
The \code{score} is a number between 0 and 1, representing the vulnerability of the taxon to perturbations for response traits, or that relate to the complexity or size for effect traits. For motility: a score from 0 (sessile) to 1 (swimmer) is given;  the higher the score, the more mobile the taxon.  
The \code{value} represents an average value for the modality class. For unitless modalities, it is equal to \code{score}. 

\code{metadata(Traits_nioz)}, \code{metadata(Traits_cefas)}, and \code{metadata(Traits_Db)}, extract the desciption for the databases; \code{attributes(Traits_nioz)$origin}, \code{attributes(Traits_cefas)$origin}, and \code{attributes(Traits_Db)$origin}, show the source of the data.


\itemize{
 \item For the *** nioz trait database*** , \code{Traits_nioz} (Beauchard et al., 2021), the following \code{effect} traits are included: 
  \itemize{
   \item{Substratum depth distribution} \code{ET1.Mx}, in cm, with 5 modalities, 0 cm (\code{ET1.M1}), 0-5 cm (\code{ET1.M2}), 5-15 cm (\code{ET1.M3}), 15-30 cm(\code{ET1.M4}), and >30 cm (\code{ET1.M5})

   \item{Biodiffusion} \code{ET2.Mx}, the random dispersion of sediment particles by organisms; 3 modalities: Null, (\code{ET2.M1}), Low (\code{ET2.M2}) and High  (\code{ET2.M3}).     

   \item{Downward conveying} \code{ET3.Mx}, downward displacement of sediment particles (including faeces production); 3 modalities: Null, (\code{ET3.M1}), Low (\code{ET3.M2}) and High  (\code{ET3.M3}).     

   \item{Upward conveying} \code{ET4.Mx}, upward displacement of sediment particles, mostly by burrying; 3 modalities: Null, (\code{ET4.M1}), Low (\code{ET4.M2}) and High  (\code{ET4.M3}).     

   \item{Regeneration} \code{ET5.Mx}, instantaneous up-and downward particle displacement, usually done by crabs; 3 modalities: Null, (\code{ET5.M1}), Low (\code{ET5.M2}) and High  (\code{ET5.M3}).     

   \item{Biodeposition} \code{ET6.Mx}, pelagic-benthic sediment transfer through suspension feeding; 3 modalities: Null, (\code{ET6.M1}), Low (\code{ET6.M2}) and High  (\code{ET6.M3}).     

   \item{Bioerosion} \code{ET7.Mx}, degradation through animal secretion; 3 modalities: Null, (\code{ET7.M1}), Low (\code{ET7.M2}) and High  (\code{ET7.M3}).     

   \item{Biostabilisation} \code{ET8.Mx}, activity that mitigates sediment erodability through hydrological modulation (zoic tubes, mounds reefs); 3 modalities: Null, (\code{ET8.M1}), Low (\code{ET8.M2}) and High  (\code{ET8.M3}).     

   \item{Ventilation/Pumping} \code{ET9.Mx}, generation of water fluxes into the substrate by body or appendix undulation; leads to bioirrigation when combined by burrows and endo-3D structure depth; 3 modalities: Null, (\code{ET9.M1}), Low (\code{ET9.M2}) and High  (\code{ET9.M3}).     

   \item{Burrow width} \code{ET10.Mx}; with 4 modalities: None,  0 mm (\code{ET10.M1}), Narrow <5mm (\code{ET10.M2}), Intermediate 5-10mm (\code{ET10.M3}) and Wide >10mm  (\code{ET10.M4}).     

   \item{Endo-3D structure type} \code{ET11.Mx}; 7 modalities: None (\code{ET11.M1}), Chimney/Funnel  (\code{ET11.M2}), Tube (\code{ET11.M3}), IJ-shaped burrow (\code{ET11.M4}), UY-shaped burrow  (\code{ET11.M5}), Branched burrow (\code{ET11.M6}) and Anastomosed burrow (\code{ET11.M7}). here chimney/funnel has a respiratory function, tube is constructed from organic or mineral material, UY-shaped differs from IJ shaped by enhancing water flux through bioturbation; similar for anastomosed and branched burrows (generally deep structures);

   \item{Endo-3D structure depth} \code{ET12.Mx}, endobenthic vertical effect size (semi- or permanent burrows), in cm; 5 modalities: None (\code{ET12.M1}), 0-5 cm (\code{ET12.M2}), 5-15 cm (\code{ET12.M3}), 15-30 cm (\code{ET12.M4}), >30 cm (\code{ET12.M5})

   \item{Epi-3D structure type} \code{ET13.Mx}, type of epibenthic structure, on top of sediment; 7 modalities: None (\code{ET13.M1}), Mat (\code{ET13.M2}), Mound (\code{ET13.M3}), Tube/Tubular protrusion (\code{ET13.M4}), Shell (\code{ET13.M5}), Stalk/Feather (\code{ET13.M6}), Protuberance (\code{ET13.M7})           

   \item{Epi-3D structure extensions} \code{ET14.xx}, type of epibenthic structure aggregations; 5 modalities: None (\code{ET14.M1}), Isolated/Clumped (\code{ET14.M2}), Mat/Lawn (\code{ET14.M3}), Simple reef (\code{ET14.M4}), Complex reef (\code{ET14.M5})
}  

The **response traits** and their modalities are: 
 
  \itemize{
   \item{Body mass} \code{RT1.Mx}, in gram ash free dry weight (gADWt); 6 modalities: <0.001  gADWt (\code{RT1.M1}), 0.001-0.010 gADWt (\code{RT1.M2}), 0.010-0.100 gADWt (\code{RT1.M3}), 0.100-1.000 gADWt (\code{RT1.M4}), 1.000-10.000 gADWt (\code{RT1.M5}), >10.000 gADWt (\code{RT1.M6}) 

   \item{Body length} \code{RT2.Mx}, in cm; 5 modalities: <1  cm (\code{RT2.M1}), 1-3 cm (\code{RT2.M2}), 3-10 cm (\code{RT2.M3}), 10-20  cm  (\code{RT2.M4}), 20-50  cm (\code{RT2.M5}) 
    
   \item{Body resistance} \code{RT3.Mx}; 5 modalities: Very low (\code{RT3.M1}), Low  (\code{RT3.M2}), Intermediate  (\code{RT3.M3}),      High (\code{RT3.M4}), Very high (\code{RT3.M5})      

   \item{Motility} \code{RT4.Mx}, 4 modalities: Sessile (\code{RT4.M1}), Tubicolous (\code{RT4.M2}), Crawler (\code{RT4.M3}), Swimmer (\code{RT4.M4})        

   \item{Burrowing/Sheltering depth} \code{RT5.Mx}, 4 modalities: 0 cm  (\code{RT5.M1}), 0-5 cm (\code{RT5.M2}), 5-15 cm (\code{RT5.M3}), >15 cm (\code{RT5.M4})
  
   \item{Feeding type} \code{RT6.Mx}; 4 modalities: Deposit feeder (\code{RT6.M1}), Suspension feeder (\code{RT6.M2}), Herbivore/Grazer (\code{RT6.M3}), Carnivore/Scavenger (\code{RT6.M4})  
    
   \item{Life span} \code{RT7.Mx}, time necessary to achieve a life cycle during which at least one reproductive success is ensured, in years; 5 modalities: < 1 yr (\code{RT7.M1}), 1-3 yr (\code{RT7.M2}), 3-10  yr (\code{RT7.M3}), 10-20 yr (\code{RT7.M4}), >20 yr (\code{RT7.M5})

   \item{Age at maturity} \code{RT8.Mx}, time after which reproductive successs can be expected, in years; 3 modalities: <1  yr (\code{RT8.M1}), 1-3 yr (\code{RT8.M2}), 3-5 yr (\code{RT8.M3})

   \item{Offspring type} \code{RT12.Mx}; 3 modalities: Egg  (\code{RT12.M1}), Larva (\code{RT12.M2}), Juvenile (\code{RT12.M3})
    
   \item{Offspring size} \code{RT13.Mx}, in mm; 4 modalities: <0.1 mm (\code{RT13.M1}),0.1-0.5 mm (\code{RT13.M2}), 0.5-1.5 mm (\code{RT13.M3}), 1.5-5 mm (\code{RT13.M4}).
}

 \item In the *** cefas trait database*** , \code{Traits_cefas} (Clare et al., 2022), the following traits are included: 
  \itemize{
   \item{Maximum size}, The maximum size (mm) that the taxon is known to reach during the adult stage, \code{sr_x}, in mm, with 5 modalities: \code{sr_Less_than_10}, \code{sr_11_to_20}, \code{sr_21_to_100}, \code{sr_101_to_200}, \code{sr_201_to_500},     \code{sr_More_than_500}
   \item{Morphology} \code{m_xx} with 6 modalities: \code{m_Soft}, \code{m_Tunic} , \code{m_Exoskeleton}, \code{m_Crustose}, \code{m_Cushion}, \code{m_Stalked} 
   \item{Lifespan}, The maximum reported lifespan (years) of the adult stage, \code{l_xx}, with 4 modalities: \code{l_Less_than_1}, \code{l_1_to_3}, \code{l_3_to_10}, \code{l_More_than_10} 
   \item{Egg development location}, \code{ed_xx}, with 4 modalities: \code{ed_Asexual}, \code{ed_Sexual_pelagic}, \code{ed_Sexual_benthic}, \code{ed_Sexual_brooded}
   \item{Larva development location}, \code{ld_xx}, 3 modalities: \code{ld_Pelagic_planktotrophic}, \code{ld_Pelagic_lecithotrophic}, \code{ld_Benthic_direct}
   \item{Living habit}, \code{lh_xx}, 5 modalities: \code{lh_Burrow_dwelling}, \code{lh_Free_living}, \code{lh_Crevice_hole_under_stones}, \code{lh_Epi_endo_biotic} , \code{lh_Attached_to_substratum} 
   \item{Sediment position}, \code{sp_xx}, 4 modalities: \code{sp_Surface}, \code{sp_Shallow_infauna_0_to_5cm}, \code{sp_Mid_depth_infauna_5_to_10cm}, \code{sp_Deep_infauna_more_than_10cm}
   \item{Feeding mode}, \code{f_xx}, 6 modalities: \code{f_Suspension}, \code{f_Surface_deposit}, \code{f_Subsurface_deposit}, \code{f_Scavenger}, \code{f_Predator}, \code{f_Parasite}
   \item{Mobility}, \code{mob_xx}, 4 modalities: \code{mob_Sessile}, \code{mob_Swim}, \code{mob_Crawl_creep_climb}, \code{mob_Burrower}
   \item{Bioturbation mode}, \code{b_xx}, 5 modalities: \code{b_Diffusive_mixing}, \code{b_Surface_deposition}, \code{b_Upward_conveyor}, \code{b_Downward_conveyer}, \code{b_None} 
 }
 
 \item The *** bioturbation trait database*** , \code{Traits_Db} (Queiros et al., 2013), contains the *Mi* and *Ri* scores, as well as the reworking types (*Fti*) that are necessary to estimate the bioturbation potential index. Their meaning is:

  \itemize{
   \item{Mi}, the mobility scores: \code{1} for organisms that live in fixed tubes; \code{2} indicates limited movement; \code{3} indicates slow, free movement through the sediment matrix; \code{4} indicates free movement, that is, via burrow system. 
   \item{Ri}, the reworking scores: \code{1} for epifauna; \code{2} for surficial modifiers; \code{3} for upward and downward conveyors; \code{4} for biodiffusors; \code{5} for regenerators. 
   \item{Fti}, the reworking types: \code{E} for epifauna; \code{S} for surficial modifiers; \code{UC} for upward conveyors; \code{DC} for downward conveyors; \code{B} for biodiffusors; \code{R} for regenerators.
 }
}
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  Olivier Beauchard
}

\examples{

##-----------------------------------------------------
## Show contents of the trait data sets - NIOZ data
##-----------------------------------------------------

traits.lab <- attributes(Traits_nioz)$description

head  (traits.lab)
unique(traits.lab$trait)

# types of burrows
subset(traits.lab, subset = (trait=="Endo-3D structure type"))

# taxa with U- or Y-shaped burrows (ET9M5):
subset(Traits_nioz, subset = ET11.M5 > 0)$taxon

##-----------------------------------------------------
## Show contents of the trait data sets - CEFAS data
##-----------------------------------------------------

traits_cefas.lab <- attributes(Traits_cefas)$description

head  (traits_cefas.lab)
unique(traits_cefas.lab$trait)

# Living habit
subset(traits_cefas.lab, subset = (trait=="Living habit"))

# taxa that live predominantly in a tube (lh_Tube_dwelling):
subset(Traits_cefas, subset = lh_Tube_dwelling > 0.5)$taxon  # more than 50%

##-----------------------------------------------------
## Feeding type - NIOZ data
##-----------------------------------------------------

FT <- subset(traits.lab, trait=="Feeding type")
FT
head(Traits_nioz[,c("taxon", FT$colname)])

# Select for all taxa the feeding type to which they belong
FTspec <- apply(Traits_nioz[ ,FT$colname], 
                MARGIN = 2,             # apply on columns
                FUN    = function(x) Traits_nioz$taxon[x>0])

names(FTspec) <- FT$modality
FTspec$"Herbivore/Grazer"

## Number of occurrences of feeding types
FTtable <- apply(Traits_nioz[,FT$colname], 
                 MARGIN = 2, 
                 FUN    = function(x) sum(x>0))

data.frame(FT$modality, numberOfOccurrences = FTtable)

##-----------------------------------------------------
## Feeding type - CEFAS data
##-----------------------------------------------------

FT_cefas <- subset(traits_cefas.lab, trait== "Feeding mode")
FT_cefas

FT_cefasSpec <- apply(Traits_cefas[ ,FT_cefas$colname], 
                MARGIN = 2,             # apply on columns
                FUN    = function(x) Traits_cefas$taxon[x>0])
                
names(FT_cefasSpec) <- FT_cefas$modality
FT_cefasSpec$"Parasite"

## number of occurrences
FTtable_cefas <- apply(Traits_cefas[,FT_cefas$colname], 
                 MARGIN = 2, 
                 FUN    = function(x) sum(x>0))

data.frame(FT_cefas$modality, numberOfOccurrences = FTtable_cefas)

##-----------------------------------------------------
## Community mean weighted score.
## Traits from all taxa in the dataset, yearly averaged 
##-----------------------------------------------------

trait.lab <- attributes(Traits_nioz)$description

# traits are extended with taxonomic information (via taxonomy), 
# and converted to crisp values (by passing trait.class, trait.score)

cwm.weight.score <- with (MWTL$density, 
      getTraitDensity(descriptor     = station, 
                      taxon          = taxon, 
                      value          = density, 
                      averageOver    = year, 
                      trait          = Traits_nioz, 
                      taxonomy       = Taxonomy,
                      trait.class    = trait.lab$trait, 
                      trait.score    = trait.lab$score, 
                      scalewithvalue = TRUE))

head(cwm.weight.score, n=c(3,4))  

Stations.traits <- merge(MWTL$stations, cwm.weight.score, 
                         by.x="station", by.y="descriptor")

par(mfrow=c(2,2))

with(Stations.traits, mapKey(x, y, colvar=Biodeposition,
                                main="Biodeposition"))
with(Stations.traits, mapKey(x, y, colvar=Biodiffusion,
                                main="Biodiffusion"))
with(Stations.traits, mapKey(x, y, colvar=Biostabilisation,
                                main="Biostabilisation"))
with(Stations.traits, mapKey(x, y, colvar=Burrow.width,
                                main="Burrow width"))

##-----------------------------------------------------
## Deriving the bio-irrigation traits from the NIOZ dataset
##-----------------------------------------------------

# the burrow types are derived from the nioz traits "Endo-3D structure type", 
# the injection depth from traits "Endo-3D structure depth" 
# the feeding type is derived from the nioz traits "Feeding type". 

# traits and modalities in Traits_nioz usable for bio-irrigation estimates
Irr.meta <- rbind(
  subset(metadata(Traits_nioz), subset = trait == "Endo-3D structure type"),  #BT
  subset(metadata(Traits_nioz), subset = trait == "Endo-3D structure depth"), #ID
  subset(metadata(Traits_nioz), subset = trait == "Feeding type"))            #FT

# Add irrigation scores
Irr.meta$IR.type  <- c(rep("BT", times=7),  rep("ID", times=5), rep("FT", times=4))
Irr.meta$IR.value <- c(1, 1, 1, 2, 2, 3, 3,     1, 2, 3, 4, 4,    3, 2, 1, 1)

# Show the metadata usable together with Traits_nioz
Irr.meta[, c(1:3,8:9)]

cnames <- c("taxon", Irr.meta$colname)  # columns to be extracted from Traits_nioz

Irr.traits <- fuzzy2crisp(trait       = Traits_nioz[, cnames], 
                          trait.class = Irr.meta$IR.type, 
                          trait.score = Irr.meta$IR.value)

# metadata
Irr.traits.meta <- data.frame(
   trait      = c(rep("BT", times=3), rep("ID", times=4),  rep("FT", times=4),
   modality   = c(1:3               , 1:4               , 1:4), 
   description= c("epifauna, internal irrigation (e.g. siphons)", 
                  "open irrigation (e.g. U- or Y- shaped burrows)", 
                  "blind ended irrigation", 
                  "injection pocket 0-2cm", "injection pocket 2-5cm", 
                  "injection pocket 5-10cm", "injection pocket >10 cm depth", 
                  "surface filter feeders", "predators", 
                  "deposit feeders", "sub surface filter feeders")))
attributes(Irr.traits)$description <- Irr.traits.meta
}

\references{

The NIOZ trait database is described in the following paper:

Beauchard O, Brind'Amour A, Schratzberger M, Laffargue P, Hintzen NT, Somerfield PJ, Piet G (2021) A generic approach to develop a trait-based indicator of trawling-induced disturbance. Mar Ecol Prog Ser 675:35-52. https://doi.org/10.3354/meps13840 

The CEFAS trait dataset:

Clare, David S. Bolam, Stefan G. McIlwaine, Paul S.O. Garcia, Clement Murray, Joanna M. Eggleton, Jacqueline D.  (2022). Ten key biological traits of marine benthic invertebrates surveyed in Northwest Europe. Cefas, UK. V2. doi: https://doi.org/10.14466/CefasDataHub.123

The bioturbation data:

Ana M. Queiros, Silvana N. R. Birchenough, Julie Bremner, Jasmin A. Godbold, Ruth E. Parker, Alicia Romero-Ramirez, Henning Reiss, Martin Solan, Paul J. Somerfield, Carl Van Colen, Gert Van Hoey, Stephen Widdicombe, 2013.  A bioturbation classification of European marine infaunal invertebrates. Ecology and Evolution 3 (11), 3958-3985

The bioirrigation traits:

A. Wrede, J.Beermann, J.Dannheim, L.Gutow, T.Brey, 2018. Organism functional traits and ecosystem supporting services - A novel approach to predict bioirrigation. Ecological indicators, 91, 737-743.

}

\seealso{
  \link{Groups} for the benthic typological groups.

  \link{MWTL} for other data sets of the Northsea.

  \link{mapKey} for plotting.

  \link{extendTrait} for functions operating on trait data.

  \link{getDensity} for functions operating on density data.

  \link{getTraitDensity} for functions operating on density and trait data.

  \link{getDbIndex} for extracting bioturbation and bioirrigation indices.
}

\keyword{datasets}

