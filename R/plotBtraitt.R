## ====================================================================
## ====================================================================
## Mapping functions for Btrait data
## ====================================================================
## ====================================================================

mapBtrait <- function(x=NULL, y=NULL, colvar=NULL, 
                      main=NULL, 
                      col=NULL,
                      colkey=list(length=0.3, width=0.5,
                         cex.axis=0.8, cex.clab=par("cex.lab")), 
                      xlim=NULL, ylim=NULL, 
                      contours=NULL, draw.levels=FALSE, 
                      col.levels=NULL, key.levels=FALSE, lwd.levels=1,
                      axes=TRUE, 
                      frame.plot=TRUE,...){
  if (is.null(contours)){
    contours <- NA # MWTL$contours
  } 
  
  if (length(contours)==1) 
    if (is.na(contours)) contours <- list(x=range(x), y=range(y), z=c(0,0))
  
  if (is.null(xlim)) xlim <- range(contours$x, na.rm=TRUE)
  
  if (is.null(ylim)) ylim <- range(contours$y, na.rm=TRUE)
  
  plotContours(x=contours$x, y=contours$y, z=contours$z, 
               xlim=xlim, ylim=ylim, 
               main=main, 
               cont=draw.levels, 
               colkey=key.levels, col=col.levels, lwd=lwd.levels, 
               hascv = ! is.null(colvar), 
               axes=axes, frame.plot=frame.plot)
  if (!is.null(x)) {
        dot <- list(...)                         
        if (!is.null(dot$log)){                   
          if (length(grep("x", dot$log)) | length(grep("y", dot$log)))
            stop ("cannot have log of 'x' or 'y-axis'")
          if (length(grep("c", dot$log)))
            colvar[colvar <= 0] <- NA
        }
        do.call("points2D",c(alist(x=x, y=y, colvar=colvar, add=TRUE, 
             colkey=colkey, col=col), dot))
  }
}

plotContours <- function (cL=NULL, x=NULL, y=NULL, z=NULL, 
                  main="", cex.main=1.5, 
                 colkey=list(dist=0.08,length=0.3, width=0.5,
                   cex.axis=0.6, cex.clab=par("cex.lab")), 
                 xlim=NULL, ylim=NULL, clim=NULL,
                 cont = FALSE, by=4, col=NULL, lwd=1, hascv=TRUE,
  axes=TRUE, frame.plot=TRUE,
  ...) { 

     if (is.null(x))  x <- cL$x
     if (is.null(y))  y <- cL$y
     if (is.null(z))  z <- cL$z     
     type <- "n"
     cv <- z
     if (! cont){
       ii <- which(cv <= 0)
       x <- x[ii]
       y <- y[ii]
       cv <- cv[ii]
       col <- "black"
       colkey <- FALSE
       type <- "l"
     }
     if (by > 1){
       ii <- sort(unique(c(seq(1, length(x), by=by), which(is.na(x)))))
       x <- x[ii]
       y <- y[ii]
       cv <- cv[ii]
     }
     if (! is.list(colkey)) if (colkey) colkey <- list()
     labs <- na.omit(unique(cv))
     dx <- diff(range(labs))/(length(labs))
     if (length(col) == 1)
       if (is.function(col)) {
         col <- col(length(labs))
         if (is.list(colkey)) {
           colkey$at <- seq(from=dx/2, by=dx, length.out=length(labs))
           colkey$labels <- as.character(labs)
           clim <- range(labs)
         }
       }
     pmar <- par("mar")
     if (! is.list(colkey)) {
       if (colkey != FALSE) pm <- par(mar=pmar+c(0,0,0,2)) 
     } else 
       pm <- par(mar=pmar+c(0,0,0,2))
     
     if (is.null(ylim)) ylim <- range(y,na.rm=TRUE)
     if (is.null(xlim)) xlim <- range(x, na.rm=TRUE)
     if (is.null(col)) {
         col<- c(grey.colors(length(labs), start=0))
         if (is.list(colkey))  {          
           colkey$at <- seq(from=dx/2, by=dx, length.out=length(labs))
           colkey$labels <- as.character(labs)
           clim <- range(labs)
         }
     }
     
     asp <- 1/cos((mean(ylim) * pi)/180)

     pm <- par(mar=par("mar"))
     if (hascv)  pm <- par(mar=par("mar")+c(0,0,0,2))
     ck <- list(plot=FALSE)
     
     points2D(x, y, colvar = cv,
              asp = asp, xlim = xlim, ylim = ylim, 
              type = type, 
              xlab = expression(""^o~E), 
              ylab = expression(""^o~N), col=col, lwd=lwd,
              main = main, cex.main = cex.main,
              las = 1, colkey = ck, 
              axes = axes, frame.plot = frame.plot) 
     if (!is.null(x) & cont) 
           lines2D(x, y, colvar=cv, add=TRUE, clim=clim,
                   colkey = colkey, 
                   clab = c("depth","m"), 
                   col = col, lwd=lwd, ...)
     par(mar=pmar)
}     
