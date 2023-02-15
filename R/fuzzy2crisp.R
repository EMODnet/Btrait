## ====================================================================
## ====================================================================
## Converts fuzzy-coded taxon x trait data to crisp values 
## (and vice versa)
## ====================================================================
## ====================================================================

fuzzy2crisp <- function(trait,             # species-trait matrix, in WIDE format
                        trait.class,       # indices to trait classes - vector
                        trait.score,       # indices to trait values - vector
                        t.column=1,        # nr or name of column(s) with taxa
                        standardize=TRUE){ # rescale to sum to 1 or not

# check input, remove descriptor column if necessary
  Trait <- trait
  trait <- clearRows(trait, t.column, "trait")
#  rown  <- row.names(trait)
#  if (attributes(trait)$numericRownames)
#    rown <- as.numeric(rown)
  cn    <- attributes(trait)$cn  # name of the column that has been removed
  rown <- Trait[,cn]
  if (length(trait.class) != length(trait.score))
    stop("'trait.class' and 'trait.score' should have the same length")
  if (ncol(trait) != length(trait.class))
    stop("'trait' should have the same nr of columns or one more, as the number of traits (length 'trait.class')")

  # weighted values
   if (is.numeric(trait.score)){
     SP <- t(t(trait)*trait.score)  # weighted  values

     # a function to convert to crisp values for one row
     FUN <- function(x)tapply(x, INDEX=list(trait.class), FUN=sum)
     
     swm <- apply(SP, MARGIN=1, FUN=FUN)
     
     if (standardize & ! is.vector(swm)){
       Sum <- apply(trait, MARGIN=1, FUN=FUN)  
       swm <- swm/Sum
       swm[is.nan(swm)] <- 0
     }
     
   } else {  # categorical trait scores
     # the number of levels (modalities) for each trait class
     LENlevel <- tapply(trait.class, INDEX=trait.class, FUN=length)
     
     # the maximum level value for each trait class
     FUN <- function(x)tapply(x, INDEX=list(trait.class), FUN=which.max)
     
     Swm <- apply(trait, MARGIN=1, FUN=FUN)
     
     swm <- trait.score[Swm+c(0, LENlevel[-length(LENlevel)])]
     dim(swm) <- dim(Swm)
     if (is.vector(swm)) swm <- matrix(data=swm, nrow=1)
     row.names(swm) <- row.names(Swm)
   }
  
  if (is.vector(swm)) {
    swm <- matrix(data=swm, ncol=1)
    colnames(swm) <-trait.class[1] 
  } else 
    swm <- t(swm)

  if (length(rown)){ 
    swm <- data.frame(rown, swm)
    colnames(swm)[t.column] <- cn
    row.names(swm) <- NULL
  }
  as.data.frame(swm)  # this will be alphabetically ordered.
}

## ====================================================================
## Converts crisp-coded taxon x trait data to fuzzy-coded values
## (makes sense for categorical variables only).
## ====================================================================

crisp2fuzzy <- function(trait,
                        t.column=1){ 

# check input
  x     <- clearRows(trait, t.column, "trait")
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
