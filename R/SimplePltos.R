## ====================================================================
## ====================================================================
## Simple plotting functions for Btrait data -> NOT YET USED
## ====================================================================
## ====================================================================

binplot <- function(x, y, 
                    nbins = length(x)/10, 
                    bins  = NULL,
                    log   = "",
                    addpts = FALSE,
                    ... ){
  if (length(grep("y", log)))
    y[y<=0] <- NA
  if (length(grep("x", log)))
    x[x<=0] <- NA
  if (is.null(bins))
    bins <- quantile(x, p=seq(0.0, 1., length.out=nbins), na.rm=TRUE)
  
  xn    <- cut(x, breaks=bins)
  
  upper <- as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(xn))) 
  lower <- as.numeric( sub("\\((.+),.*", "\\1", levels(xn)) )
  xmean <- colMeans(rbind(upper, lower))   
  
  xup <- as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", xn)) 
  xlo <- as.numeric( sub("\\((.+),.*", "\\1", xn) )
  xm <- colMeans(rbind(xup, xlo))   
  BPX <- boxplot(y~xn, log=log, at=xmean, xlim=range(x, na.rm=TRUE), 
    boxwex=min(upper-lower)/2, names=xmean, axes=FALSE, frame.plot=TRUE, ...)
  axis(side=1)
  axis(side=2)
  
  if (addpts) {
    xp <- jitter(xm, amount=min(diff(xmean))/10)
    points(xp, y, pch=".", cex=3, col="red")
  }
  abline(v=bins, col="grey", lty=2)
  rug(x, ticksize=0.01)
  BPX$x    <- xmean
  BPX$bins <- bins
  invisible(BPX)
}

