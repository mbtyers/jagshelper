#' Trace plot of jagsUI object
#' @description Trace plot of a whole `jagsUI` object, or optional subset of parameter nodes.
#' @param x Posterior `jagsUI` object
#' @param p Parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
#' @param exact Whether `p` should be an exact match (`TRUE`) or just match the
#' beginning of the string (`FALSE`).  Defaults to `FALSE`.
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param lwd Line width for plotting.  Defaults to 1.
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{pairstrace_jags}, \link{trace_df}, \link{trace_line}
#' @author Matt Tyers
#' @examples
#'
#' trace_jags(asdf_jags_out, parmfrow=c(4,2))
#' trace_jags(asdf_jags_out, p="a", parmfrow=c(3,1))
#' @export
trace_jags <- function(x,p=NULL,exact=FALSE,parmfrow=NULL,lwd=1,...) {
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
    plot(NA, xlim=c(1,nrow(x_plist[[i]])), ylim=range(x_plist[[i]]), main=names(x_plist)[i], xlab="iter", ylab="", ...=...)
    for(j in 1:nline) lines(x_plist[[i]][,j], col=cols[j], lwd=lwd)
  }

  # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}



#' Combination of trace plots and by-chain kernel densities of `jagsUI` object
#' @description Combination of trace plots and by-chain kernel densities of a whole `jagsUI` object, or optional subset of parameter nodes.
#' @param x Posterior `jagsUI` object
#' @param p Parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
#' @param exact Whether `p` should be an exact match (`TRUE`) or just match the
#' beginning of the string (`FALSE`).  Defaults to `FALSE`.
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param lwd Line width for plotting.  Defaults to 1.
#' @param shade Whether to add semi-transparent shading to by-chain kernel densities.  Defaults to `TRUE`.
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{trace_jags}, \link{chaindens_jags}, \link{pairstrace_jags}
#' @author Matt Tyers
#' @examples
#'
#' tracedens_jags(asdf_jags_out, parmfrow=c(4,2))
#' tracedens_jags(asdf_jags_out, p="a", parmfrow=c(3,1))
#' @export
tracedens_jags <- function (x, p = NULL, exact=FALSE, parmfrow = NULL, lwd = 1, shade=TRUE, ...)
{
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  suppressWarnings(x_plist <- jags_plist(x, p = p, exact=exact))
  if(length(x_plist)==0) stop("No parameters with matching names")
  if (!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow = parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }
  nline <- ncol(x_plist[[1]])
  cols <- adjustcolor(rainbow(nline), red.f = 0.9, blue.f = 0.9,
                      green.f = 0.9, alpha.f = 0.6)
  for (i in 1:length(x_plist)) {
    allbw <- density(as.vector(x_plist[[i]]))$bw
    denses <- apply(x_plist[[i]], 2, density, bw = allbw)
    dx <- sapply(denses, function(x) x$x)
    dy <- sapply(denses, function(x) x$y)
    dy1 <- (dy/max(dy) * 0.3 + 1) * nrow(x_plist[[i]])
    plot(NA, xlim = c(1, 1.3 * nrow(x_plist[[i]])), ylim = range(x_plist[[i]]),
         main = names(x_plist)[i], xlab = "iter", ylab = "",
         ... = ...)
    for (j in 1:nline) lines(x_plist[[i]][, j], col = cols[j],
                             lwd = lwd)
    for (j in 1:nline) lines(dy1[, j], dx[, j], col = cols[j],
                             lwd = lwd)
    if(shade) {
      for (j in 1:nline) polygon(dy1[, j], dx[, j], border = NA, #cols[j],
                                 col=adjustcolor(cols[j], alpha.f=.2))
    }
  }
  # if (!is.null(parmfrow))
  #   par(mfrow = parmfrow1)
}


#' Simple trace plot
#' @description Trace plot of a single parameter node.
#' @param x Posterior vector
#' @param nline Number of MCMC chains
#' @param lwd Line width
#' @param main Plot title
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{trace_jags}, \link{trace_df}, \link{chaindens_line}
#' @author Matt Tyers
#' @examples
#' b1 <- jags_df(asdf_jags_out, p="b1")
#'
#' trace_line(b1, nline=3, main="b1")
#' @export
trace_line <- function(x, nline, lwd=1, main="", ...) {
  # if(is.null(nline)) nline <- length(x)/n
  if(inherits(x,"data.frame")) x <- x[[1]]  #####

  n <- length(x)/nline
  cols <- adjustcolor(rainbow(nline),red.f=.9,blue.f=.9,green.f=.9,alpha.f=.6)
  plot(NA,xlim=c(0,n),ylim=range(x,na.rm=T),main=main,xlab="iteration",ylab="value", ...=...)
  for(i in 1:nline) {
    lines(1:n, x[(n*(i-1)+1):(n*i)], col=cols[i],lwd=lwd)
  }
}


#' Trace plot of each column of a `data.frame`.
#' @description Trace plot of each column of a posterior 'data.frame'.
#' @param df Posterior data.frame
#' @param nline Number of chains
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments or arguments to `trace_line()`
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{trace_jags}, \link{trace_line}
#' @author Matt Tyers
#' @examples
#' a <- jags_df(asdf_jags_out, p="a")
#'
#' trace_df(a, nline=3, parmfrow=c(3,1))
#' @export
trace_df <- function(df, nline, parmfrow=NULL, ...) {
  if(!inherits(df,c("data.frame","matrix"))) stop("Input must be a data.frame")
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }
  thencol <- ifelse(is.null(ncol(df)), 1, ncol(df))
  if(thencol==1) {
    trace_line(df,nline=nline,...=...)
  } else {
    for(i in 1:ncol(df)) {
      # trace_line(df[,i],main=names(df)[i],nline=nline,...=...)
      trace_line(df[,i],main=colnames(df)[i],nline=nline,...=...)
    }}
  # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' Pairs trace plot
#' @description Two-dimensional trace plots (or alternately, scatter plots or contour plots) of each possible pair of
#' parameters from a possible subset.  May be useful in assessing correlation between parameter nodes, or problematic
#' posterior surfaces.
#' @param x Output object returned from `jagsUI`
#' @param p Optional vector of parameters to subset
#' @param points Whether to plot as scatter plots instead.  Defaults to `FALSE`.
#' @param contour Whether to plot as contour plots instead.  Defaults to `FALSE`.
#' @param lwd Line width for trace plots.  Defaults to 1.
#' @param alpha Opacity of lines (or points, when `points=TRUE`).  Defaults to 0.2.
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments or arguments to `tracedens_jags()`
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{crossplot}
#' @author Matt Tyers
#' @examples
#' pairstrace_jags(SS_out, p="sig", parmfrow=c(2,3), lwd=2)
#' pairstrace_jags(SS_out, p="sig", parmfrow=c(2,3), points=TRUE)
#' pairstrace_jags(SS_out, p="sig", parmfrow=c(2,3), contour=TRUE)
#'
#' pairstrace_jags(asdf_jags_out, parmfrow=c(3,3))
#' pairstrace_jags(asdf_jags_out, parmfrow=c(3,3), points=TRUE)
#' pairstrace_jags(asdf_jags_out, parmfrow=c(3,3), contour=TRUE)
#' @export
pairstrace_jags <- function (x, p = NULL, points=FALSE, contour=FALSE, lwd = 1, alpha=0.2, parmfrow = NULL,...)
{
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")

  x_plist <- jags_plist(x, p = p)
  if (!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow = parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }
  nline <- ncol(x_plist[[1]])
  cols <- adjustcolor(rainbow(nline), red.f = 0.9, blue.f = 0.9,
                      green.f = 0.9, alpha.f = alpha)
  nn <- length(x_plist)
  for (i in 1:(nn-1)) {
    for(j in (i+1):nn) {
      plot(NA, xlim=range(x_plist[[i]], na.rm=T), ylim=range(x_plist[[j]]),
           main="", xlab=names(x_plist)[i], ylab=names(x_plist)[j],...=...)
      if(!contour) {
        for (k in 1:nline) {
          if(!points) {
            lines(x=x_plist[[i]][, k], y=x_plist[[j]][, k], col = cols[k], lwd = lwd)
          } else {
            points(x=x_plist[[i]][, k], y=x_plist[[j]][, k], col = cols[k])
          }
        }
      } else {
        xx <- as.vector(x_plist[[i]])
        yy <- as.vector(x_plist[[j]])
        # if(sum(diff(xx))!=0 & sum(diff(yy))!=0) {
        contour(MASS::kde2d(x=xx, y=yy), main="", xlab=names(x_plist)[i], ylab=names(x_plist)[j], add=T)
        # } else {
        # plot(NA, main="", xlab=names(x_plist)[i], ylab=names(x_plist)[j], xlim=0:1, ylim=0:1)
        # }
      }
    }
  }
  # if (!is.null(parmfrow))
  #   par(mfrow = parmfrow1)
}
