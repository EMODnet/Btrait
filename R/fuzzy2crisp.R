## ====================================================================
## ====================================================================
## Converts fuzzy-coded taxon x trait data to crisp values 
## (and vice versa)
## ====================================================================
## ====================================================================

fuzzy2crisp <- function(trait,             # species-trait matrix, in WIDE format
                        trait_class,       # indices to trait classes - vector
                        trait_score,       # indices to trait values - vector
                        taxon_column=1,    # nr or name of column(s) with taxa
                        standardize=TRUE){ # rescale to sum to 1 or not

# check input, remove descriptor column if necessary
  Trait <- trait
  trait <- clearRows(trait, taxon_column, "trait")
#  rown  <- row.names(trait)
#  if (attributes(trait)$numericRownames)
#    rown <- as.numeric(rown)
  cn    <- attributes(trait)$cn  # name of the column that has been removed
  rown <- Trait[,cn]
  if (length(trait_class) != length(trait_score))
    stop("'trait_class' and 'trait_score' should have the same length")
  if (ncol(trait) != length(trait_class))
    stop("'trait' should have the same nr of columns or one more, as the number of traits (length 'trait_class')")

  # weighted values
   if (is.numeric(trait_score)){
     SP <- t(t(trait)*trait_score)  # weighted  values

     # a function to convert to crisp values for one row
     FUN <- function(x)tapply(x, INDEX=list(trait_class), FUN=sum)
     
     swm <- apply(SP, MARGIN=1, FUN=FUN)
     
     if (standardize & ! is.vector(swm)){
       Sum <- apply(trait, MARGIN=1, FUN=FUN)  
       swm <- swm/Sum
       swm[is.nan(swm)] <- 0
     }
     
   } else {  # categorical trait scores
     # the number of levels (modalities) for each trait class
     LENlevel <- tapply(trait_class, INDEX=trait_class, FUN=length)
     
     # the maximum level value for each trait class
     FUN <- function(x)tapply(x, INDEX=list(trait_class), FUN=which.max)
     
     Swm <- apply(trait, MARGIN=1, FUN=FUN)
     
     swm <- trait_score[Swm+c(0, LENlevel[-length(LENlevel)])]
     dim(swm) <- dim(Swm)
     if (is.vector(swm)) swm <- matrix(data=swm, nrow=1)
     row.names(swm) <- row.names(Swm)
   }
  
  if (is.vector(swm)) {
    swm <- matrix(data=swm, ncol=1)
    colnames(swm) <-trait_class[1] 
  } else 
    swm <- t(swm)

  if (length(rown)){ 
    swm <- data.frame(rown, swm)
    colnames(swm)[taxon_column] <- cn
    row.names(swm) <- NULL
  }
  as.data.frame(swm)  # this will be alphabetically ordered.
}

## ====================================================================
## Converts crisp-coded taxon x trait data to fuzzy-coded values
## (makes sense for categorical variables only).
## ====================================================================

crisp2fuzzy <- function(trait,
                        taxon_column=1){ 

# check input
  x     <- clearRows(trait, taxon_column, "trait")
  rown  <- row.names(x)

  type <- sapply(x, class)
  cnx  <- colnames(x)
  out  <- as.data.frame(rown)
  cn   <- names(trait)[1]
  
  iscategory <-TRUE
  for (i in 1:ncol(x)){
    CC <- x[,i]
    if (type[i] %in% c("factor", "character", "logical")){
      if (type[i]!="factor") CC <- as.factor(CC)
      out <- data.frame(out, 
                        sapply(levels(CC), 
                               FUN=function(x) as.integer(x == CC)))
      cn <- c(cn, paste(cnx[i],levels(CC), sep="_"))
    } else{
      iscategory <- FALSE
      out <- data.frame(out, CC)
      cn <- c(cn, cnx[i])
    }
  }
  names(out) <- cn
  
  if (iscategory){
     TT <- apply(x, MARGIN=2, FUN=function(X)length(unique(X)))
     trait <- rep(names(TT), times=TT)
     modality <- unlist(apply(x, MARGIN=2, 
                              FUN=function(X)levels(as.factor(X))))
     description <- data.frame(colname=cn[-1], trait=trait, modality=modality)
     attributes(out)$description <- description
   } 
    out
}
