## ====================================================================
## ====================================================================
## Mapping functions for Btrait data
## ====================================================================
## ====================================================================

## ====================================================================
## Mapping with a color key
## ====================================================================

  
mapKey <- function(x=NULL, y=NULL, colvar=NULL, main=NULL, col=NULL, lwd=1,
                   colkey=list(), xlim=NULL, ylim=NULL, clim=NULL, 
                   contours=NULL, draw.levels=FALSE, col.levels="black", 
                   by.levels=1, key.levels=FALSE, lwd.levels=1,
                   axes=TRUE, frame.plot=TRUE, asp=NULL, NApch=43, ...){

  hasContours <- TRUE
  if (is.logical(contours)){
    hasContours <- contours
    contours <- NA # so that its length = 1
  }
  if (is.null(contours)){
    hasContours <- FALSE
    contours <- NA # so that its length = 1
  } 

  if (length(contours)==1) 
    if (is.na(contours)) contours <- list(x=range(x), y=range(y), z=c(0,0))
  
  if (is.null(xlim)) xlim <- range(contours$x, na.rm=TRUE)
  if (is.null(ylim)) ylim <- range(contours$y, na.rm=TRUE)
  
  # open a plot, add contours
  plotContours(x=contours$x, y=contours$y, z=contours$z, 
               xlim=xlim, ylim=ylim, main=main, addCont=hasContours,
               colkey=key.levels, draw.levels=draw.levels, by=by.levels,
               col=col.levels, lwd=lwd.levels, hascv = ! is.null(colvar), 
               axes=axes, frame.plot=frame.plot, asp=asp, clim=clim)
  
  # Now plot the data
  if (!is.null(x)) {
     dot <- list(...)  
     cv  <- colvar
        if (!is.null(dot$log)){                   
          if (is.null(asp) & (length(grep("x", dot$log)) | length(grep("y", dot$log))))
            stop ("cannot have log of 'x' or 'y-axis' when asp=NULL")
          if (length(grep("c", dot$log))) cv <- log(colvar)
        }
     isna  <- which(is.na(cv) | is.infinite(cv))
     if (length(isna)){
       xisna <- x[isna]
       yisna <- y[isna]
       x     <- x[-isna]
       y     <- y[-isna]
       colvar <- colvar[-isna]
       NAcol <- dot$NAcol
       if (is.null(NAcol)) NAcol <- "black"
     }  
       if (length(colvar))
       do.call("points2D",c(alist(x=x, y=y, colvar=colvar, add=TRUE, 
             colkey=colkey, lwd=lwd, col=col), dot))
       if (length(isna))
          points2D(xisna, yisna, pch=NApch, add=TRUE, col=NAcol)
  }
}

plotContours <- function (x=NULL, y=NULL, z=NULL, addCont=TRUE,
                 xlim=NULL, ylim=NULL, main="", cex.main=1.5, 
                 colkey=list(dist=0.08,length=0.3, width=0.5,
                    cex.axis=0.6, cex.clab=par("cex.lab")), 
                 draw.levels=FALSE, by=4, col=NULL, lwd=1, hascv=TRUE,
                 axes=TRUE, frame.plot=TRUE, asp=NULL, clim=NULL, ...) { 

     type <- "n"
     cv   <- z
     if (addCont){
      if (! draw.levels){  # Only 0-contours
        ii     <- which(cv <= 0)
        x      <- x[ii]
        y      <- y[ii]
        cv     <- cv[ii]
        col    <- "black"
        colkey <- FALSE
        type   <- "l"
      }
     
      if (by > 1){
       ii <- sort(unique(c(seq(1, length(x), by=by), which(is.na(x)))))
       x <- x[ii]
       y <- y[ii]
       cv <- cv[ii]
      }
     } else {
       x <- NULL
       y <- NULL
       cv <- NULL
     }
     if (! is.list(colkey)) 
       if (colkey) colkey <- list()
     
     labs <- na.omit(unique(cv))
     dx   <- diff(range(labs))/(length(labs))
     if (length(col) == 1)
       if (is.function(col)) {
         col <- col(length(labs))
         
         if (is.list(colkey)) {
           colkey$at <- seq(from=dx/2, by=dx, length.out=length(labs))
           colkey$labels <- as.character(labs)
           if (is.null(clim)) clim <- range(labs)
         }
       }
     
     pmar <- par("mar")
     if (! is.list(colkey)) {
       if (colkey != FALSE) pm <- par(mar=pmar+c(0,0,0,2)) 
     } else 
       pm <- par(mar=pmar+c(0,0,0,2))
     
     if (is.null(ylim)) ylim <- range(y, na.rm=TRUE)
     if (is.null(xlim)) xlim <- range(x, na.rm=TRUE)
     
     if (is.null(col)) {
         col<- c(grey.colors(length(labs), start=0))
         if (is.list(colkey))  {          
           colkey$at <- seq(from=dx/2, by=dx, length.out=length(labs))
           colkey$labels <- as.character(labs)
           if (is.null(clim)) clim <- range(labs)
         }
     }
     
     if (is.null(asp)) 
       asp <- 1/cos((mean(ylim) * pi)/180)  # aspect ratio for map

     pm <- par(mar=par("mar"))
     if (hascv)  pm <- par(mar=par("mar")+c(0,0,0,2))
     
     ck <- list(plot=FALSE)
     
     points2D(x, y, colvar=cv,
              asp=asp, xlim=xlim, ylim=ylim, clim=clim, type=type, 
              xlab=expression(""^o~E), ylab=expression(""^o~N), 
              col=col, lwd=lwd, main=main, cex.main=cex.main,
              las=1, colkey=ck, axes=axes, frame.plot=frame.plot) 
     
     if (!is.null(x) & draw.levels) 
           lines2D(x, y, colvar=cv, add=TRUE, clim=clim,
                   colkey = colkey, 
                   clab=c("depth","m"), 
                   col=col, lwd=lwd, ...)
     par(mar=pmar)
}     


## ====================================================================
## Mapping with a legend; size of pch ~ value
## ====================================================================

mapLegend <- function (x=NULL, y=NULL, colvar=NULL, main=NULL, col=NULL, lwd=1,
                       legend = list(), xlim=NULL, ylim=NULL, clim=NULL, 
                       scale="simple", pch=18, cex=3, cex.min=cex/20, 
                       contours=NULL, draw.levels=FALSE, col.levels="black", 
                       by.levels=1, key.levels=FALSE, lwd.levels=1, 
                       axes=TRUE, frame.plot=TRUE, asp=NULL, 
                       NApch=43, NAtext="NA", ...){
  
  hasContours <- ! is.null(contours)
  ell <- list(...)
  
  if (is.null(contours)) contours <- NA 
  
  if (length(contours)==1) 
    if (is.na(contours)) contours <- list(x=range(x), y=range(y), z=c(0,0))
  
  if (is.null(xlim)) xlim <- range(contours$x, na.rm=TRUE)
  if (is.null(ylim)) ylim <- range(contours$y, na.rm=TRUE)
  
  if (is.null(x))  x <- contours$x
  if (is.null(y))  y <- contours$y
  if (is.null(colvar))  colvar <- 1:x     

  # check for log transformation
  if (is.null(ell$log))
    islog <- FALSE
  else
    islog <- grepl("c", ell$log)
  
  cv <- colvar
  if (islog) cv <- log10(colvar)
  
  cl <- clim
  if (! is.null(clim)) {
    if (islog) cl <- log10(clim)
    colvar[colvar > max(clim)] <- max(clim)
    colvar[colvar < min(clim)] <- min(clim)
  }
  
  isna <- which(is.na(cv) | is.infinite(cv))
  if (length(isna)){
    xisna <- x[isna]
    yisna <- y[isna]
    x     <- x[-isna]
    y     <- y[-isna]
    colvar <- colvar[-isna]
    cv     <- cv[-isna]
    NAcol <- ell$NAcol
    if (is.null(NAcol)) NAcol <- "black"
  }
  # size range for symbols
  if (cex.min < cex) 
    cex <- cex - cex.min 
  else 
    cex.min <- 0
  
  # size and pch for all colvar
  Pch <- pch[1]
  
  # all same size
  pr <- 0
  
  # for legend
  # Legend
  if (! is.list(legend)) legend <- as.list(legend)
  
  if (is.null(legend$x))
    legend$x <- "center"

  if (! length(cv)){
    legend$pch <- NApch
    legend$legend <- NAtext
    col <- NAcol
    Pch <- NApch
    Cex <- 1
    x     <- xisna
    y     <- yisna
    isna  <- NULL
    colvar <- rep(1, times=length(x))
  }  
  if (length(cv)) {  
   pcv <- pretty(cv)
  
  # scale according to values
   if (scale == "simple") {  
    pr <- range(pcv)  # range
    dr <- diff(pr)
    CexFun <- function(x) 
      cex.min + cex*(x-pr[1])/dr
    
    # scale according to absolute values   
   } else if (scale == "abs"){   
    pr <- range(pcv)  # range
    dr <- diff(pr)
    
    # min value, so that it makes sense for all pos/neg values
    dm <- ifelse (sign(prod(pr)) >= 0, pr[1], max(min(pcv), 0)) 
    
    dr <- diff(c(dm, max(abs(pcv))))
    CexFun <- function(x) cex.min + cex*(abs(x)-dm)/dr
    
    if (length(pch) == 2){  # different pch for pos and neg values
      Pch <- rep(pch[1], length.out=length(colvar))
      Pch[colvar<0] <-pch[2]
    }
   } else CexFun <- function(x) cex   # no scaling
  
   Cex <- CexFun(cv)
  
  
   legend$col <- createKey(pcv, col=col)
   legend$pt.lwd <- lwd
   
   legend$pch <- rep(pch[1], length=length(pcv))
   if (scale == "abs" & length(pch) == 2)  # pch has 2 values: for pos/neg
    legend$pch[pcv < 0] <- pch[2]
  
   legend$pt.cex <- CexFun(pcv)
   legend$pt.bg <- ell$bg
   legend$legend <- pcv
   
   
   }   
   
   legend.side <- legend$side
   if (is.null(legend.side)) legend.side <- 4
   legend.cex  <- legend$cex
   if (is.null(legend.cex)) legend.cex <- 1
   legend.pars <- legend$pars
  
   legend$side <- legend$cex <- legend$pars <- NULL
  
   if (islog) {
     mp  <- min(pcv)
     pcv <- 10^pcv
     legend$legend <- format(pcv, digits=max(1, -mp))
  
    legend$title <- ell$clab
   }
   if ( length(isna)){
     if (! class(NApch) == class(legend$pch)) 
       stop ("cannot merge NApch and pch as one is numeric other character")
     legend$pch <- c(legend$pch, NApch)
     legend$legend <- c(legend$legend,NAtext)
     legend$pt.cex <- c(legend$pt.cex, 1)
     legend$col <- c(legend$col, NAcol)
   }  
  
  asp <- NULL
  if (is.null(asp)) asp <- 1/cos((mean(ylim) * pi)/180)
  sc <- list(method="scatter2D", x=x, y=y, colvar=colvar, colkey=FALSE, 
             main=main, legend=legend, col=col, lwd=lwd, pch=Pch, cex=Cex, 
             xlim=xlim, ylim=ylim, clim=clim,
             frame.plot=frame.plot, axes=axes, legend.side=legend.side, 
             legend.cex=legend.cex, legend.pars=legend.pars)
  sc <- c(alist(...), sc)
  sc$asp <- asp
  do.call("legend.plt", sc)
  #do.call("lplt", sc)
  if (hasContours) {
    x <- contours$x
    y <- contours$y
    z <- contours$z
    if (by.levels > 1){
      ii <- sort(unique(c(seq(1, length(x), by=by.levels), which(is.na(x)))))
      x <- x[ii]
      y <- y[ii]
      z <- z[ii]
    }
    clab <- contours$clab
    if (is.null(clab)) clab <- c("depth","m")
    lines2D(x, y, colvar=z, add=TRUE, #clim=NULL,
            colkey = key.levels, clab = clab, 
            col = col.levels, lwd=lwd.levels)
  }
  if (length(isna))
    points2D(xisna, yisna, pch=NApch, add=TRUE, col="black")
  
}

mapMWTL <- function(x=NULL, y=NULL, colvar=NULL, draw.levels=TRUE, 
                    type="key", ...){
  if (type=="key")
    mapKey(x=x,y=y,colvar=colvar, contours=MWTL$contours, 
           draw.levels=draw.levels, 
           colkey=list(length=0.3, width=0.5, shift=-0.2, dist=-0.2,
                        cex.axis=0.8, cex.clab=par("cex.lab")),
           las=1, ...)
  else
    mapLegend(x=x, y=y, colvar=colvar, contours=MWTL$contours,
              draw.levels=draw.levels, legend=list(x="bottomright", side=0), 
              las=1, ...)
  
}

mapNSBS <- function(x=NULL, y=NULL, colvar=NULL, draw.levels=TRUE, 
                    type="key", ...){
  if (type=="key")
    mapKey(x=x,y=y,colvar=colvar, contours=NSBS$contours, 
           draw.levels=draw.levels, 
           colkey=list(length=0.3, width=0.5,shift=-0.2, dist=-0.2,
                       cex.axis=0.8, cex.clab=par("cex.lab")), 
           las=1, ...)
  else
    mapLegend(x=x, y=y, colvar=colvar, contours=NSBS$contours,
              draw.levels=draw.levels, legend=list(x="bottomright", side=0),
              las=1, ...)
  
}

