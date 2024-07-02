#' By-chain kernel densities of `jagsUI` object
#' @description By-chain kernel densities of a whole `jagsUI` object, or optional subset of parameter nodes.
#' @param x Posterior `jagsUI` object
#' @param p Parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
#' @param exact Whether `p` should be an exact match (`TRUE`) or just match the
#' beginning of the string (`FALSE`).  Defaults to `FALSE`.
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param lwd Line width for plotting.  Defaults to 1.
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{trace_jags}, \link{chaindens_line}, \link{chaindens_df}
#' @author Matt Tyers
#' @examples
#'
#' chaindens_jags(asdf_jags_out, parmfrow=c(4,2))
#' chaindens_jags(x=asdf_jags_out, p="a", parmfrow=c(3,1))
#' @export
chaindens_jags <- function(x,p=NULL,exact=FALSE,parmfrow=NULL,lwd=1,...) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  suppressWarnings(x_plist <- jags_plist(x, p = p, exact=exact))
  if(length(x_plist)==0) stop("No parameters with matching names")
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }

  nline <- ncol(x_plist[[1]])
  cols <- adjustcolor(rainbow(nline),red.f=.9,blue.f=.9,green.f=.9,alpha.f=.6)
  for(i in 1:length(x_plist)) {
    allbw <- density(as.vector(x_plist[[i]]))$bw
    denses <- apply(x_plist[[i]], 2, density, bw=allbw)
    dx <- sapply(denses, function(x) x$x)
    dy <- sapply(denses, function(x) x$y)
    plot(NA, xlim=range(dx), ylim=range(dy), main=names(x_plist)[i], xlab="",ylab="",...=...)
    for(j in 1:nline) lines(dx[,j], dy[,j], col=cols[j], lwd=lwd)
  }

  # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' Simple by-chain kernel density plot
#' @description By-chain kernel density plot of a single parameter node.
#' @param x Posterior vector
#' @param nline Number of chains
#' @param lwd Line width
#' @param main Plot title
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{chaindens_jags}, \link{chaindens_df}
#' @author Matt Tyers
#' @examples
#' b1 <- jags_df(asdf_jags_out, p="b1")
#'
#' chaindens_line(b1, nline=3, main="b1")
#' @export
chaindens_line <- function(x, nline, lwd=1, main="", ...) {
  # if(is.null(nline)) nline <- length(x)/n
  if(inherits(x,"data.frame")) x <- x[[1]]  #####

  n <- length(x)/nline
  cols <- adjustcolor(rainbow(nline),red.f=.9,blue.f=.9,green.f=.9,alpha.f=.6)
  allbw <- density(x)$bw
  denses <- list()
  for(i in 1:nline) {
    denses[[i]] <- density(x[(n*(i-1)+1):(n*i)], bw=allbw)
  }
  dx <- sapply(denses, function(x) x$x)
  dy <- sapply(denses, function(x) x$y)
  plot(NA,xlim=range(dx),ylim=range(dy),main=main,xlab="",ylab="", ...=...)
  for(i in 1:nline) {
    lines(dx[,i], dy[,i], col=cols[i],lwd=lwd)
  }
}


#' By-chain kernel density of each column of a `data.frame`.
#' @description By-chain kernel density plot of each column of a posterior `data.frame`.
#' @param df Posterior `data.frame`
#' @param nline Number of chains
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments or arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{trace_jags}, \link{trace_line}
#' @author Matt Tyers
#' @examples
#' a <- jags_df(asdf_jags_out, p="a")
#'
#' chaindens_df(a, nline=3, parmfrow=c(3,1))
#' @export
chaindens_df <- function(df, nline, parmfrow=NULL, ...) {
  if(!inherits(df,c("data.frame","matrix"))) stop("Input must be a data.frame")
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }
  for(i in 1:ncol(df)) {
    # chaindens_line(df[,i],main=names(df)[i],nline=nline,...=...)
    chaindens_line(df[,i],main=colnames(df)[i],nline=nline,...=...)
  }
  # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' Compare Density
#' @description Side-by-side kernel density plots for all parameters (or a specified subset) from two `jagsUI`
#' output objects or `data.frame`s.  The intent of this function is easy comparison of inferences from two comparable models.
#'
#' Kernel densities are plotted vertically, either left- or right-facing.  Parameters with the same name are
#'  plotted facing one another.
#' @param x1 Output object returned from `jagsUI`; or alternately a `data.frame`
#' @param x2 Second output object returned from `jagsUI`; or alternately a `data.frame`
#' @param p Optional vector of parameters to subset.  All parameters with names matching the beginning of the
#' string supplied will be returned.  If the default (`NULL`) is accepted, all parameters will be plotted.
#' @param minCI Minimum CI width for plotting.  This is intended as a method for excluding far-outlying MCMC
#' samples when determining the appropriate y-axis limits for plotting.  Defaults to 99%.
#' @param ylim Y-axis limits for plotting.  If the default (`NULL`) is accepted, limits will be automatically determined.
#' @param legendnames Names for optional legend.  If the default `NULL` is accepted, no legend will be drawn.
#' @param legendpos Position for optional legend.  Defaults to `"topleft"`.
#' @param col Colors for kernel density plots.  Defaults to colors `4` and `2` (blue and red).
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{comparecat}, \link{comparepriors}
#' @author Matt Tyers
#' @examples
#' ## This is the same output object twice, but shows functionality.
#' comparedens(x1=asdf_jags_out, x2=asdf_jags_out, p=c("a","b","sig"),
#'             legendnames=c("Model 1", "Model 2"))
#' @export
comparedens <- function(x1,x2, p=NULL, minCI=0.99, ylim=NULL, legendnames=NULL, legendpos="topleft", col=c(4,2), ...) {
  # if(!inherits(x1,"jagsUI") | !inherits(x2,"jagsUI")) {
  #   stop("Inputs must both a output objects returned from jagsUI::jags().")
  # }

  if(!((inherits(x1,"jagsUI") | inherits(x1,"data.frame")) & (inherits(x2,"jagsUI") | inherits(x2,"data.frame")))) {
    stop("Inputs must be data.frames or output objects returned from jagsUI::jags().")
  }

  if(inherits(x1,"jagsUI")) {
    xdf1 <- jags_df(x1)
  } else {
    xdf1 <- x1
  }
  if(inherits(x2,"jagsUI")) {
    xdf2 <- jags_df(x2)
  } else {
    xdf2 <- x2
  }

  if(!is.null(p)) {
    parmx1 <- pull_post(xdf1,p=p)
    parmx2 <- pull_post(xdf2,p=p)
  } else {
    parmx1 <- xdf1
    parmx2 <- xdf2
  }

  # allparms <- sort(unique(c(names(parmx1),names(parmx2))))
  allparms <- unique(c(names(parmx1),names(parmx2)))

  cilo1 <- apply(parmx1, 2, quantile, p=(1-minCI)/2)
  cilo2 <- apply(parmx2, 2, quantile, p=(1-minCI)/2)
  cihi1 <- apply(parmx1, 2, quantile, p=1-((1-minCI)/2))
  cihi2 <- apply(parmx2, 2, quantile, p=1-((1-minCI)/2))
  if(is.null(ylim)) ylim <- range(cilo1,cihi1,cilo2,cihi2,na.rm=T)

  plot(NA,xlim=c(0,length(allparms)+1), ylim=ylim, ylab="",xlab="",xaxt="n",...=...)
  axis(1,at=1:length(allparms),labels=allparms,las=2)
  abline(v=1:length(allparms),lty=3)

  if(length(col)==1) col <- rep(col, 2)
  for(i in 1:length(allparms)) {
    if(allparms[i] %in% names(parmx1)) {
      xxx <- density(parmx1[,which(names(parmx1)==allparms[i])])
      polygon(i-xxx$y/max(xxx$y)*.4, xxx$x, col=adjustcolor(col[1],alpha.f=.5), border=col[1])
    }
    if(allparms[i] %in% names(parmx2)) {
      xxx <- density(parmx2[,which(names(parmx2)==allparms[i])])
      polygon(i+xxx$y/max(xxx$y)*.4, xxx$x, col=adjustcolor(col[2],alpha.f=.5), border=col[2])
    }
  }

  if(!is.null(legendnames)) {
    legend(legendpos,legend=legendnames, fill=adjustcolor(col, alpha.f=.5), border=col)
  }
}





#' Plot kernel densities of single parameter nodes
#' @description Produces a kernel density plot of a single or multiple parameter nodes (overlayed).
#'
#' Input can be of multiple possible formats: either a single or list of output objects
#'  from `jagsUI` with an associated vector of parameter names, or a vector or `data.frame`
#'  of posterior samples.
#' @param df Input object for plotting.  See examples below.
#' @param p Vector of parameter names, if `df` is given as a single or list of output objects
#'  from `jagsUI`
#' @param exact Whether the `p=` argument should match the parameter name exactly.  See
#' \link{jags_df} for details.
#' @param add Whether to add to an existing plot (`TRUE`) or produce a new plot.
#' Defaults to `FALSE`.
#' @param col Vector of colors for plotting.  If the default (`NULL`) is accepted,
#' colors will be automatically selected.
#' @param shade Whether to shade the regions below the kernel density curve(s).
#' Defaults to `TRUE`.
#' @param lwd Line width for kernel density curves.  Defaults to `2`.  Note: setting
#' this to `0` (or `FALSE`) will suppress lines.
#' @param minCI Minimum CI width to include for all density curves.  Defaults to 99%.
#' @param legend Whether to plot a legend.  Defaults to `TRUE`.
#' @param legendpos Position for automatic legend.  Defaults to `"topleft"`.
#' @param legendnames Names for legend
#' @param main Plot title.  Defaults to "".
#' @param xlab X-axis label.  Defaults to "".
#' @param ylab Y-axis label.  Defaults to "Density".
#' @param ... Optional plotting arguments
#' @return `NULL`
#' @seealso \link{comparedens}, \link{comparecat}, \link{comparepriors}
#' @author Matt Tyers
#' @examples
#' ## jagsUI object with a single parameter
#' plotdens(asdf_jags_out, p="b1")
#'
#' ## jagsUI object with multiple nodes of a parameter
#' plotdens(asdf_jags_out, p="a")
#'
#' ## jagsUI object with multiple parameter nodes
#' plotdens(asdf_jags_out, p=c("a[1]","a[2]","a[3]"))
#'
#' ## data.frame with multiple columns
#' plotdens(jags_df(asdf_jags_out, p="a"))
#'
#' ## list of jagsUI objects with a single parameter name
#' plotdens(list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p="b1")
#'
#' ## list of jagsUI objects with a vector of parameter names
#' plotdens(list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p=c("a[1]","a[2]","a[3]"))
#' @export
plotdens <- function(df, p=NULL, exact=FALSE, add=FALSE,
                     col=NULL, shade=TRUE, lwd=2, minCI=0.99,
                     legend=TRUE, legendpos="topleft", legendnames=NULL,
                     main=NULL, xlab="", ylab="Density",...) {  #...
  dflist <- NULL # initial instance to check later

  ## - single density
  # vector
  # df <- jags_df(asdf_jags_out, p="b0")
  # if(inherits(df, "data.frame")) {
  #   if(ncol(df)==1) {
  #     dflist <- list(df[,1])
  #   }
  # }
  # df <- jags_df(asdf_jags_out, p="b0")[,1]
  if(inherits(df,"numeric") & is.null(dim(df))) {
    dflist <- list(df)
  }

  # df <- jags_df(asdf_jags_out, p="a")
  if(inherits(df, "data.frame")) {
    dflist <- as.list(df)
  }

  # df <- as.matrix(jags_df(asdf_jags_out, p="a"))
  if(inherits(df, "matrix")) {
    dflist <- as.list(as.data.frame(df))
  }

  # jags & p
  # df <- asdf_jags_out
  # p <- "a"
  if(inherits(df, "jagsUI")) {
    dflist1 <- jags_df(df, p=p, exact=exact)  # this can be multiple
    dflist <- as.list(dflist1)
  }

  ## - multiple densities
  # df <- list(asdf_jags_out, asdf_jags_out, asdf_jags_out)
  # p <- c("a[1]","a[2]","a[3]")
  # p <- "b1"
  if(inherits(df,"list")) {
    # list of jags, list of p (same length)
    if(length(df)==length(p)) {
      dflist <- list()
      for(i in 1:length(df)) {
        suppressWarnings(dflist[[i]] <- tryCatch(jags_df(x=df[[i]], p=p[i], exact=T)[,1], error=function(e) NA))
        if(all(is.na(dflist[[i]]))) stop("No parameter names are an exact match to p= argument.")
      }
    } else {
      # list of jags, one p
      if(length(p)==1) {
        dflist <- list()
        for(i in 1:length(df)) {
          suppressWarnings(dflist[[i]] <- tryCatch(jags_df(x=df[[i]], p=p, exact=T)[,1], error=function(e) NA))
          if(all(is.na(dflist[[i]]))) stop("No parameter names are an exact match to p= argument.")
        }
      } else {
        stop("p= argument must be of length 1 or length equal to the length of jagsUI object list")
      }
    }

  }

  if(is.null(dflist)) stop("No allowable input detected.  See help(plotdens) for details.")
  # alldims <- sapply(dflist, dim)
  # for(i in 1:length(alldims)) {
  #   if(!is.null(alldims[[i]])) stop("Multiple columns detected for at least one list element")
  # }

  # make xlims
  los <- sapply(dflist, quantile, p=0.5*(1-minCI), na.rm=T)
  his <- sapply(dflist, quantile, p=1-(0.5*(1-minCI)), na.rm=T)
  xlims <- range(los, his, na.rm=T)

  # make denses
  denses <- lapply(dflist, density)
  xs <- sapply(denses, function(x) x$x)
  ys <- sapply(denses, function(x) x$y)

  # plot denses
  if(is.null(main) & length(p)==1) main <- p
  if(is.null(col)) col <- c(4,2,3,rcolors(100))[1:length(dflist)]
  if(!add) {
    plot(NA, xlim=xlims, ylim=range(ys), xlab=xlab, ylab=ylab, main=main,...=...)
  }
  for(i in 1:length(dflist)) {
    if(shade) {
      polygon(x=c(xs[1,i], rev(xs[,i])),
              y=c(ys[1,i], rev(ys[,i])),
              border=NA, col=adjustcolor(col[i],alpha.f=.3))
    }
    lines(denses[[i]], col=col[i], lwd=lwd)
  }

  # legend
  if(!is.null(names(dflist)) & is.null(legendnames)) legendnames <- names(dflist)
  if(legend & !is.null(legendnames)) {
    legend(legendpos, lwd=lwd, col=col, legend=legendnames)
  }
}
# plotdens(df=asdf_jags_out, p="b1")
# plotdens(asdf_jags_out, p="a")
# plotdens(jags_df(asdf_jags_out, p="a"))
# plotdens(df=list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p="a")
# plotdens(df=list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p="a[1]")
# plotdens(list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p="b1",lwd=F)
