#' Skeleton
#' @description Provides a paste-able skeleton of code with an example JAGS model and output.
#' @param NAME Name to append to JAGS model object, etc.
#' @return `NULL`
#' @note The printed code will use the `cat()` function to write the model code to an
#' external text file.  It may be desirable to use a call to `\link{tempfile}()`
#' instead, to eliminate creation of unneeded files.
#' @author Matt Tyers
#' @importFrom grDevices adjustcolor rainbow
#' @importFrom graphics axis lines par points polygon segments
#' @importFrom stats density quantile
#' @importFrom grDevices rgb
#' @importFrom graphics abline legend
#' @importFrom stats rbeta median
#' @examples
#' skeleton("asdf")
#' @export
skeleton <- function(NAME="NAME") {
cat("
library(jagsUI)

# specify model, which is written to an external file
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0 + b1*x[i] + a[grp[i]]
  }

  for(j in 1:ngrp) {
    a[j] ~ dnorm(0, tau_a)
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  tau_a <- pow(sig_a, -2)
  sig_a ~ dunif(0, 10)
}', file=\"",NAME,"_jags\")


# simulate data to go with the example model
n <- 60
x <- rnorm(n, sd=3)
grp <- sample(1:3, n, replace=T)
y <- rnorm(n, mean=grp-x)

# bundle data to pass into JAGS
",NAME,"_data <- list(x=x,
                 y=y,
                 n=length(x),
                 grp=as.numeric(as.factor(grp)),
                 ngrp=length(unique(grp)))

# JAGS controls
niter <- 10000
ncores <- 3

{
  tstart <- Sys.time()
  print(tstart)
  ",NAME,"_jags_out <- jagsUI::jags(model.file=\"",NAME,"_jags\", data=",NAME,"_data,
                             parameters.to.save=c(\"b0\",\"b1\",\"sig\",\"a\",\"sig_a\"),
                             n.chains=ncores, parallel=T, n.iter=niter,
                             n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}",sep="")
}

#' Extract data.frame
#' @description Extracts the posterior from `jagsUI` output in the form of
#' a `data.frame`.  This simpler construction has a few benefits: operations may
#' be more straightforward, posterior objects will be smaller files, and can be
#' written to an external table or .csv, etc.
#' @param x Output object from `jagsUI::jags()`
#' @return A `data.frame` with a column associated with each parameter and a row
#' associated with each MCMC iteration.
#' @seealso \link{pull_post}
#' @author Matt Tyers
#' @import jagsUI
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#' @export
jags_df <- function(x) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(as.data.frame(as.matrix(x$samples)))
}


#' Subset from posterior data.frame
#' @description Extracts a subset vector or `data.frame` from a larger (likely posterior) `data.frame`, such that column names match a name given in the `p=` argument.
#' @param x Posterior `data.frame`
#' @param p String to begin posterior names.  If `NULL` is used, all parameters will be returned.
#' @param exact Whether name must be an exact match (`TRUE`) or with initial sub-string matching only supplied characters (`FALSE`).  Defaults to `FALSE`.
#' Defaults to `FALSE.`
#' @return A `data.frame` with a column associated with each (subsetted) parameter and a row
#' associated with each MCMC iteration.
#' @seealso \link{jags_df}
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b <- pull_post(out_df, p="b")
#' str(b)
#' a <- pull_post(out_df, p=c("a","sig_a"))
#' str(a)
#' sigs <- pull_post(out_df, p="sig")
#' str(sigs)
#' justsig <- pull_post(out_df, p="sig", exact=TRUE)
#' str(justsig)
#' @export
pull_post <- function(x, p=NULL, exact=FALSE) {
  # x[,substr(names(x),1,nchar(p))==p]
  if(!inherits(x,"data.frame")) stop("Input must be a data.frame")

  if(is.null(p)) {
    these <- rep(T, length(names(x)))
  } else {
    these <- rep(F, length(names(x)))  ## used to be F
    for(i in 1:length(p)) {
      if(exact) {
        these[names(x)==p[i]] <- T
      } else {
        these[substr(names(x),1,nchar(p[i]))==p[i]] <- T
      }
    }
    if(sum(these)==0) warning("No parameters with matching names, returning empty data.frame")
  }
  return(x[these])
}


#' Plist
#' @description Extract a list of nrep X nchain matrices, one for each parameter.
#' @param x `jagsUI` output object
#' @param p String to subset parameter names, if a subset is desired
#' @return A `list` with an element associated with each parameter.  Each element
#' will be a matrix with a column associated with each MCMC chain and a row for
#' each MCMC iteration.
#' @note It is unlikely that a user will need this function; it is included
#' primarily as a helper function used by other functions in this package.
#' @author Matt Tyers
#' @examples
#' out_plist <- jags_plist(asdf_jags_out)
#' str(out_plist)
#'
#' a_plist <- jags_plist(asdf_jags_out, p=c("a","sig_a"))
#' str(a_plist)
#' @export
jags_plist <- function(x, p=NULL) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")

  x_dflist <- lapply(x$samples, as.data.frame)
  x2 <- lapply(1:length(x_dflist[[1]]), function(x) sapply(x_dflist,
                                                           "[[", x))
  names(x2) <- names(x_dflist[[1]])
  these <- rep(F, length(x2))
  if (!is.null(p)) {
    for(i in 1:length(p)) {
      these[substr(names(x2), 1, nchar(p[i])) == p[i]] <- T
    }
    x2 <- x2[these]
  }
  if(length(x2)==0) warning("No parameters with matching names, returning empty list")
  return(x2)
}

#' Trace plot of jagsUI object
#' @description Trace plot of a whole `jagsUI` object, or subset.
#' @param x Posterior `jagsUI` object
#' @param p Parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
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
trace_jags <- function(x,p=NULL,parmfrow=NULL,lwd=1,...) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  suppressWarnings(x_plist <- jags_plist(x, p = p))
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


#' By-chain kernel densities of `jagsUI` object
#' @description By-chain kernel densities of a whole `jagsUI` object, or subset.
#' @param x Posterior `jagsUI` object
#' @param p Parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
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
chaindens_jags <- function(x,p=NULL,parmfrow=NULL,lwd=1,...) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  suppressWarnings(x_plist <- jags_plist(x, p = p))
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


#' Combination of trace plots and by-chain kernel densities of `jagsUI` object
#' @description Combination of trace plots and by-chain kernel densities of a whole `jagsUI` object, or subset.
#' @param x Posterior `jagsUI` object
#' @param p Parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
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
tracedens_jags <- function (x, p = NULL, parmfrow = NULL, lwd = 1, shade=TRUE, ...)
{
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  suppressWarnings(x_plist <- jags_plist(x, p = p))
  if(length(x_plist)==0) stop("No parameters with matching names")
  if (!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow = parmfrow)
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
  if (!is.null(parmfrow))
    par(mfrow = parmfrow1)
}

#' Number of parameters
#' @description Number of individual parameter nodes saved in `jagsUI` output.
#' @param x Output object from `jagsUI::jags()`
#' @return A single numeric value giving the number of parameter nodes.
#' @seealso \link{nbyname}
#' @author Matt Tyers
#' @examples
#' head(jags_df(asdf_jags_out))
#'
#' nparam(asdf_jags_out)
#' @export
nparam <- function(x) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(ncol(jags_df(x)))
}


#' Number of parameter nodes by parameter name
#' @description Returns a summary of the numbers of parameter nodes saved in `jagsUI` output, by parameter name.
#' As a default, what is returned for each list element is a vector of the array dimensions within the JAGS model
#'  (excluding the number of MCMC iterations),
#' or alternately, just the total number of parameters
#' @param x Output object from `jagsUI::jags()`
#' @param justtotal Whether to just report the total number of parameters, as opposed to dimensions.
#' @return A `list` with an element associated with each parameter.  Each element
#' can be interpreted as the vector length or array dimension associated with the
#' given parameter.
#' @seealso \link{nparam}
#' @author Matt Tyers
#' @examples
#' head(jags_df(asdf_jags_out))
#'
#' nbyname(asdf_jags_out)
#'
#' nparam(SS_out)
#' nbyname(SS_out)
#' @export
nbyname <- function(x, justtotal=FALSE) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  # sapply(x$sims.list, function(x) ncol(as.matrix(x)))
  out1 <- sapply(x$sims.list, function(x) dim(x)[-1])
  out1[sapply(out1, is.null)] <- 1
  if(!justtotal) return(out1) #comment
  if(justtotal) return(sapply(out1, prod))
}


#' Quick summary of Rhat values by parameter name
#' @description Returns the mean number of `Rhat` values (by each parameter) that are less than a specified threshold criterion.
#' @param x Output object from `jagsUI::jags()`
#' @param thresh Threshold value (defaults to 1.1)
#' @return Numeric (named) giving the proportion of Rhat values below the given threshold.
#' @seealso \link{check_neff}, \link{traceworstRhat}, \link{plotRhats}
#' @author Matt Tyers
#' @examples
#' check_Rhat(SS_out)
#' @export
check_Rhat <- function(x, thresh=1.1) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(sapply(x$Rhat, function(x) mean(x<thresh, na.rm=TRUE)))
}


#' Quick summary of n.eff values by parameter name
#' @description Returns the mean number of `n.eff` values (by each parameter) that are greater than a specified threshold criterion.
#' @param x Output object from `jagsUI::jags()`
#' @param thresh Threshold value (defaults to 500)
#' @return Numeric (named) giving the proportion of `n.eff` values above the given threshold.
#' @seealso \link{check_Rhat}
#' @author Matt Tyers
#' @examples
#' check_neff(SS_out)
#' @export
check_neff <- function(x, thresh=500) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(sapply(x$n.eff, function(x) mean(x>thresh, na.rm=TRUE)))
}


#' Logit
#' @description Logit log(x/(1-x)
#' @param x Numeric vector
#' @return Numeric vector
#' @seealso \link{expit}
#' @author Matt Tyers
#' @examples
#' logit(0.5)
#' @export
logit <- function(x) log(x/(1-x))


#' Expit, or inverse logit
#' @description Inverse logit, where logit is defined as log(x/(1-x).
#'
#' Expit (inverse logit) is defined as exp(x)/(1+exp(x)).
#' @param x Numeric vector
#' @return Numeric vector
#' @seealso \link{logit}
#' @author Matt Tyers
#' @examples
#' expit(0)
#' @export
expit <- function(x) exp(x)/(1+exp(x))




#' Simple trace plot
#' @description Trace plot of a single parameter.
#' @param x Posterior vector
#' @param nline Number of MCMC chains
#' @param lwd Line width
#' @param main Plot title
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{trace_jags}, \link{trace_df}, \link{chaindens_line}
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
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

#' Simple by-chain kernel density plot
#' @description By-chain kernel density plot of a single parameter.
#' @param x Posterior vector
#' @param nline Number of chains
#' @param lwd Line width
#' @param main Plot title
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{chaindens_jags}, \link{chaindens_df}
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
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

#' Trace plot of each column of a `data.frame`
#' @description Trace plot of each column of a posterior 'data.frame'.
#' @param df Posterior data.frame
#' @param nline Number of chains
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments or arguments to `trace_line()`
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{trace_jags}, \link{trace_line}
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
#'
#' par(mfrow=c(3,1))
#' trace_df(a, nline=3)
#'
#' trace_df(a, nline=3, parmfrow=c(3,1))
#' @export
trace_df <- function(df, nline, parmfrow=NULL, ...) {
  if(!inherits(df,"data.frame")) stop("Input must be a data.frame")
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
    trace_line(df[,i],main=names(df)[i],nline=nline,...=...)
  }}
  # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}

#' By-chain kernel density of each column of a `data.frame`
#' @description By-chain kernel density plot of each column of a posterior `data.frame`.
#' @param df Posterior `data.frame`
#' @param nline Number of chains
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments or arguments
#' @return `NULL`
#' @seealso \link{tracedens_jags}, \link{trace_jags}, \link{trace_line}
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
#'
#' par(mfrow=c(3,1))
#' trace_df(a, nline=3)
#'
#' chaindens_df(a, nline=3, parmfrow=c(3,1))
#' @export
chaindens_df <- function(df, nline, parmfrow=NULL, ...) {
  if(!inherits(df,"data.frame")) stop("Input must be a data.frame")
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }
  for(i in 1:ncol(df)) {
    chaindens_line(df[,i],main=names(df)[i],nline=nline,...=...)
  }
  # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' Envelope plot
#' @description Envelope plot of the posterior densities of a vector of parameter nodes,
#' in which the sequential order of nodes is important, such as a time series.
#' @param df Output object returned from `jagsUI::jags()`; or alternately,
#' two-dimensional `data.frame` or matrix in which parameter node element is
#' given by column and MCMC iteration is given by row.
#' @param p Parameter name, if input to `df` is a `jagsUI` output object.
#' @param x Vector of X-coordinates for plotting.
#' @param row Row to subset, in the case of a 2-d matrix of parameter nodes in-model.
#' @param column Column to subset, in the case of a 2-d matrix of parameter nodes in-model.
#' @param median Whether to include median line
#' @param ci Vector of intervals to overlay.  Defaults to 50 percent and 95 percent.
#' @param col Color for plotting
#' @param add Whether to add to existing plot
#' @param dark Opacity (0-1) for envelopes.  Note that multiple overlapping intervals will darken the envelope.
#' @param outline Whether to just envelope outlines
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param main Plot title.  If the default (`NULL`) is accepted and argument p is used, p will be used for the title.
#' @param ylim Y-axis limits for plotting.  If the default (`NULL`) is accepted, these will be determined automatically.
#' @param ... additional plotting arguments or arguments to `lines()`
#' @return `NULL`
#' @seealso \link{overlayenvelope}, \link{caterpillar}
#' @author Matt Tyers
#' @examples
#' ## usage with input data.frame
#' SS_df <- jags_df(SS_out)
#' trend <- pull_post(SS_df, "trend")
#' envelope(trend, x=SS_data$x)
#'
#' ## usage with jagsUI object
#' envelope(SS_out, p="trend")
#'
#' ## usage with 2-d jagsUI object
#' envelope(SS_out, p="cycle_s", column=1, main="cycle")
#' envelope(SS_out, p="cycle_s", column=2, col=2, add=TRUE)  ## overlay
#' @export
envelope <- function(df,
                     p=NULL,
                     x=NA,
                     row=NULL, column=NULL,
                     median=TRUE,
                     ci=c(0.5,0.95),
                     col=4, add=FALSE, dark=.3, outline=FALSE,
                     xlab="", ylab="", main=NULL,
                     ylim=NULL, ...) {
  # ci <- sort(ci)
  # loq <- apply(df, 2, quantile, p=(1-ci)/2, na.rm=T)
  # hiq <- apply(df, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  # if(length(ci)==1) {
  #   loq <- t(as.matrix(loq))
  #   hiq <- t(as.matrix(hiq))
  # }
  # med <- apply(df, 2, median, na.rm=T)
  # if(all(is.na(x))) x <- 1:ncol(df)
  # if(!add) {
  #   if(is.null(ylim)) ylim <- range(loq,hiq,na.rm=T)
  #   plot(NA, xlim=range(x),ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...=...)
  #   if(median) lines(x, med, col=col)
  # }
  # else
  #   if(median) lines(x, med, type='l', col=col, ...=...)
  # if(outline) {
  #   darks <- rev(1-((1-dark)^(1:length(ci))))
  #   for(i in 1:length(ci)) {
  #     lines(x,loq[i,], col=adjustcolor(col,alpha.f=darks[i]), ...=...)
  #     lines(x,hiq[i,], col=adjustcolor(col,alpha.f=darks[i]), ...=...)
  #   }
  # }
  # else {
  #   for(i in 1:length(ci)) polygon(c(x,rev(x)), c(loq[i,],rev(hiq[i,])), col=adjustcolor(col,alpha.f=dark), border=NA)
  # }

  if(!inherits(df,"jagsUI") & !(inherits(df,"data.frame") | inherits(df,"matrix"))) {
    stop("Input must be a data.frame or output from jagsUI::jags() plus parameter name")
  }
  if(inherits(df,"jagsUI") & length(p)!=1) stop("Need single parameter name in p= argument") ###

  if(inherits(df,"jagsUI") & !is.null(p)) {
  # if(class(df)[1]=="jagsUI" & !is.null(p)) {
    simslist <- df$sims.list
    if(all(names(simslist)!=p)) stop("No parameters with matching names") ###
    theparm <- simslist[names(simslist)==p][[1]]
    if(length(dim(theparm))==2) {
      df <- theparm
    } else {
      if(!is.null(row)) df <- theparm[,row,]
      if(!is.null(column)) df <- theparm[,,column]
    }
    if(is.null(main)) {
      if(is.null(row) & is.null(column)) {
        main <- p
      } else {
        if(!is.null(row)) main <- paste0(p,"[",row,",]")
        if(!is.null(column)) main <- paste0(p,"[,",column,"]")
      }
    }
  }
  if(is.null(main)) main <- ""

  ci <- sort(ci)
  loq <- apply(df, 2, quantile, p=(1-ci)/2, na.rm=T)
  hiq <- apply(df, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  if(length(ci)==1) {
    loq <- t(as.matrix(loq))
    hiq <- t(as.matrix(hiq))
  }
  med <- apply(df, 2, median, na.rm=T)
  if(all(is.na(x))) x <- 1:ncol(df)
  if(!add) {
    if(is.null(ylim)) ylim <- range(loq,hiq,na.rm=T)
    plot(NA, xlim=range(x),ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...=...)
    if(median) lines(x, med, col=col)
  }
  else
    if(median) lines(x, med, type='l', col=col, ...=...)
  if(outline) {
    darks <- rev(1-((1-dark)^(1:length(ci))))
    for(i in 1:length(ci)) {
      lines(x,loq[i,], col=adjustcolor(col,alpha.f=darks[i]), ...=...)
      lines(x,hiq[i,], col=adjustcolor(col,alpha.f=darks[i]), ...=...)
    }
  }
  else {
    for(i in 1:length(ci)) polygon(c(x,rev(x)), c(loq[i,],rev(hiq[i,])), col=adjustcolor(col,alpha.f=dark), border=NA)
  }
}


#' Overlay envelope plots
#' @description Overlays multiple envelope plots of posterior `data.frames`, or outputs returned from `jagsUI`.
#' This would be best suited to a set of posterior `data.frames` or 2-d matrices representing sequential vectors of parameter nodes.
#' @param df Primary input can be specified in a number of ways: either a `list()` of posterior `data.frame`s or matrices,
#' a `list` of output objects returned from `jagsUI::jags()`, a 3-dimensional `array` in which the input matrices to plot
#' are separated according to the 3rd array dimension, or a single output object returned from `jagsUI::jags()` with multiple
#' arguments passed to `p`, following.
#' @param p Parameter name, if input to `df` is a list of `jagsUI` output objects; or a vector of parameter names, if
#' input to `df` is a single `jagsUI` output object.
#' @param x Optional vector of X-coordinates for plotting.
#' @param row Row to subset, in the case of a 2-d matrix of parameter nodes in-model.
#' @param column Column to subset, in the case of a 2-d matrix of parameter nodes in-model.
#' @param median Whether to include median line
#' @param ci Vector of intervals to overlay.  Defaults to 50 percent and 95 percent.
#' @param col Vector of colors for plotting
#' @param add Whether to add to existing plot
#' @param dark Opacity (0-1) for envelopes.  Note that multiple overlapping intervals will darken the envelope.  Defaults to 0.3.
#' @param outline Whether to just envelope outlines
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param main Plot title.  If the default (`NULL`) is accepted and argument `p=` is used, `p` will be used for the title.
#' @param ylim Y-axis limits for plotting.  If the default (`NULL`) is accepted, these will be determined automatically.
#' @param legend Whether to automatically try to add a legend.  Defaults to `TRUE`.
#' @param legendnames Optional vector of names for a legend.
#' @param ... additional plotting arguments or arguments to `lines()`
#' @return `NULL`
#' @seealso \link{envelope}
#' @author Matt Tyers
#' @examples
#' ## usage with list of input data.frames
#' overlayenvelope(df=list(SS_out$sims.list$cycle_s[,,1],
#'                             SS_out$sims.list$cycle_s[,,2]))
#'
#' ## usage with a 3-d input array
#' overlayenvelope(df=SS_out$sims.list$cycle_s)
#'
#' ## usage with a jagsUI output object and parameter name (2-d parameter)
#' overlayenvelope(df=SS_out, p="cycle_s")
#'
#' ## usage with a single jagsUI output object and multiple parameters
#' overlayenvelope(df=SS_out, p=c("trend","rate"))
#' @export
overlayenvelope <- function(df,
                            p=NULL,
                            x=NA,
                            row=NULL, column=NULL,
                            median=TRUE,
                            ci=c(0.5,0.95),
                            col=NULL, add=FALSE, dark=.3, outline=FALSE,
                            xlab="", ylab="", main=NULL,
                            ylim=NULL,
                            legend=TRUE, legendnames=NULL, ...) {
  conditionmet <- F
  ## list of dfs
  # do nothing
  if(inherits(df, "list")) {
    if(inherits(df[[1]],"data.frame") | inherits(df[[1]],"matrix")) {
      conditionmet <- T
    }
  }

  ## list of jagsUIs and parameter name
  if(inherits(df, "list")) {
    if(inherits(df[[1]], "jagsUI")) {
      plist <- list()
      for(i in 1:length(df)) {
        plist[[i]] <- df[[i]]$sims.list[names(df[[i]]$sims.list)==p][[1]]
      }
      df <- plist
      conditionmet <- T
    }
  }



  ## single jagsUI, vector of parameter names
  if(inherits(df, "jagsUI") & length(p)>1) {
    plist <- list()
    for(i in 1:length(p)) {
      plist[[i]] <- df$sims.list[names(df$sims.list)==p[i]][[1]]
    }
    df <- plist
    conditionmet <- T
  }

  ## single jagsUI and a 3d parameter
  if(inherits(df, "jagsUI") & length(p)==1) {
    pp <- df$sims.list[names(df$sims.list)==p][[1]]
    if(length(dim(pp))==3) df <- pp
    conditionmet <- T
  }

  ## 3d array
  if(inherits(df,"array")) {
    if(length(dim(df))==3) {
      df <- apply(df, 3, function(x) x, simplify=F)
      conditionmet <- T
    }
  }

  ## list of jagsUIs and parameter name - and a column name
  ## single jagsUI, vector of parameter names - and a column
  for(i in 1:length(df)) {
    if(length(dim(df[[i]]))>2) {
      if(!is.null(row)) df[[i]] <- df[[i]][,row,]
      if(!is.null(column)) df[[i]] <- df[[i]][,,column]
      conditionmet <- T
    }
  }

  if(!conditionmet) stop("Invalid input.  See ?overlayenvelope for details.")

  #### do stuff
  cilim <- c((1-max(ci))/2, 1-(1-max(ci))/2)
  bounds <- sapply(df, function(x) apply(x,2, quantile, p=cilim, na.rm=T))
  if(is.null(ylim)) ylim <- range(bounds)

  if(length(col) != length(df)) {
    cols <- c(4,2,3,rcolors(100))[1:length(df)]  ## betterize this??
  } else {
    cols <- col
  }

  if(is.null(main)) {
    if(length(p)==1) {
      main <- p
      if(!is.null(row)) main <- paste0(p,"[",row,",]")
      if(!is.null(column)) main <- paste0(p,"[,",column,"]")
    }
  }
  envelope(df[[1]], ci=ci, col=cols[1], add=add, ylim=ylim, main=main, x=x,
           median=median, dark=dark, outline=outline,xlab=xlab,ylab=ylab, ...=...)
  for(i in 2:length(df)) {
    envelope(df[[i]], ci=ci, col=cols[i], add=T, x=x, median=median, dark=dark,
             outline=outline, ...=...)
  }

  ## add legend??
  if(length(p)==length(df) & is.null(legendnames)) {
    legendnames <- p
    if(!is.null(row)) {
      for(i in 1:length(p)) {
        legendnames[i] <- paste0(p[i],"[",row,",]")
      }
    }
    if(!is.null(column)) {
      for(i in 1:length(p)) {
        legendnames[i] <- paste0(p[i],"[,",column,"]")
      }
    }
  }
  if(legend & !is.null(legendnames)) {
    legend("topleft",legend=legendnames, fill=adjustcolor(cols, alpha.f=.5), border=cols)
  }
}



#' Caterpillar plot
#' @description Caterpillar plot of the posterior densities of a vector of parameter nodes,
#' in which the sequential order of nodes might not be important, such as vector of random effects.
#' @param df Output object returned from `jagsUI::jags()`; or alternately,
#' two-dimensional `data.frame` or matrix in which parameter node element is
#' given by column and MCMC iteration is given by row.
#' @param p Parameter name, if input to `df` is a `jagsUI` output object.
#' @param x Vector of X-coordinates for plotting.
#' @param row Row to subset, in the case of a 2-d matrix of parameter nodes in-model.
#' @param column Column to subset, in the case of a 2-d matrix of parameter nodes in-model.
#' @param median Whether to include medians
#' @param mean Whether to include means
#' @param ci Vector of intervals to overlay.  Defaults to 50 percent and 95 percent.
#' @param col Color for plotting
#' @param add Whether to add to existing plot
#' @param lwd Base line width for plotting.  Defaults to 1.
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param main Plot title.  If the default (`NULL`) is accepted and argument `p` is used, `p` will be used for the title.
#' @param xax Vector of possible x-axis tick labels.  Defaults to the `data.frame` column names.
#' @param medlwd Line width of median line
#' @param medwd Relative width of median line.  Defaults to 1, perhaps smaller numbers will look better?
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{envelope}
#' @author Matt Tyers
#' @examples
#' ## usage with input data.frame
#' data(asdf_jags_out)
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
#'
#' caterpillar(a)
#' caterpillar(a, ci=seq(.1,.9,by=.1))
#' caterpillar(a, lwd=2)
#' caterpillar(a, xax=c("effect 1", "effect 2", "effect 3"))
#'
#'
#' ## usage with input as jagsUI object
#' caterpillar(asdf_jags_out, p="a")
#' caterpillar(SS_out, p="rate")
#'
#' ## usage with a 2-d parameter matrix
#' caterpillar(SS_out, p="cycle_s", column=1)
#' caterpillar(SS_out, p="cycle_s", column=2)
#' @export
caterpillar <- function(df,
                        p=NULL,
                        x=NA,
                        row=NULL, column=NULL,
                        median=TRUE, mean=FALSE,
                        ci=c(0.5,0.95),
                        lwd=1, col=4, add=FALSE,
                        xlab="", ylab="", main=NULL,
                        xax=NA,
                        medlwd=lwd, medwd=1,...) {
  # ci <- rev(sort(ci))
  # loq <- apply(df, 2, quantile, p=(1-ci)/2, na.rm=T)
  # hiq <- apply(df, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  # if(length(ci)==1) {
  #   loq <- t(as.matrix(loq))
  #   hiq <- t(as.matrix(hiq))
  # }
  # med <- apply(df, 2, median, na.rm=T)
  # if(all(is.na(x))) x <- 1:ncol(df)
  # d <- diff(x[1:2])
  # nn <- ncol(df)
  # if(all(is.na(xax))) xax<-names(df)
  # lwds <- (1+2*(1:length(ci)-1))*lwd
  # if(!add) {
  #   plot(NA, type='l', ylim=range(loq,hiq,na.rm=T), xlim=range(x-(.2*d),x+(.2*d)), xlab=xlab, ylab=ylab, main=main, xaxt="n", ...=...)
  #   axis(1,x,labels=xax)
  # }
  # if(median) {
  #   segments(x0=x-.2*d*medwd,x1=x+.2*d*medwd,y0=med,y1=med,col=col,lwd=medlwd, lend=1)
  # }
  # if(mean) points(x,colMeans(df, na.rm=T), pch=16, col=col)
  # for(i in 1:length(ci)) segments(x0=x,x1=x,y0=loq[i,],y1=hiq[i,],col=col,lwd=lwds[i],lend=1)

  if(!inherits(df,"jagsUI") & !inherits(df,"data.frame")) {
    stop("Input must be a data.frame or output from jagsUI::jags() plus parameter name")
  }
  if(inherits(df,"jagsUI") & length(p)!=1) stop("Need single parameter name in p= argument") ###

  if(class(df)[1]=="jagsUI" & !is.null(p)) {
    simslist <- df$sims.list
    if(all(names(simslist)!=p)) stop("No parameters with matching names") ###
    theparm <- simslist[names(simslist)==p][[1]]
    if(is.null(row) & is.null(column)) {
      df <- theparm
    } else {
      theparm <- simslist[names(simslist)==p][[1]]
      if(!is.null(row)) df <- theparm[,row,]
      if(!is.null(column)) df <- theparm[,,column]
    }
    if(is.null(main) & length(p)==1) {
      if(is.null(row) & is.null(column)) {
        main <- p
      } else {
        if(!is.null(row)) main <- paste0(p,"[",row,",]")
        if(!is.null(column)) main <- paste0(p,"[,",column,"]")
      }
    }
  }
  if(is.null(main)) main <- ""

  ci <- rev(sort(ci))
  loq <- apply(df, 2, quantile, p=(1-ci)/2, na.rm=T)
  hiq <- apply(df, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  if(length(ci)==1) {
    loq <- t(as.matrix(loq))
    hiq <- t(as.matrix(hiq))
  }
  med <- apply(df, 2, median, na.rm=T)
  if(all(is.na(x))) x <- 1:ncol(df)
  d <- diff(x[1:2])
  nn <- ncol(df)
  if(all(is.na(xax))) xax<-names(df)
  lwds <- (1+2*(1:length(ci)-1))*lwd
  if(!add) {
    plot(NA, type='l', ylim=range(loq,hiq,na.rm=T), xlim=range(x-(.2*d),x+(.2*d)), xlab=xlab, ylab=ylab, main=main, xaxt="n", ...=...)
    axis(1,x,labels=xax)
  }
  if(median) {
    segments(x0=x-.2*d*medwd,x1=x+.2*d*medwd,y0=med,y1=med,col=col,lwd=medlwd, lend=1)
  }
  if(mean) points(x,colMeans(df, na.rm=T), pch=16, col=col)
  for(i in 1:length(ci)) segments(x0=x,x1=x,y0=loq[i,],y1=hiq[i,],col=col,lwd=lwds[i],lend=1)
}


#' Example data: asdf jags out
#'
#' A simple model, equivalent to that produced by the output produced by `\link{skeleton}`.
#'
"asdf_jags_out"

#' Example data: SS JAGS out
#'
#' A time series model with multiple observations of a single time series, and with two stochastic cycle components.
#'
#' This model is included partly to show a model with vectors or 2-dimensional
#' matrices of parameter nodes, and also to give an example of poor model convergence.
#'
"SS_out"

#' Example data: Time series associated with SS JAGS out
#'
#' The time series and time measurements associated with the time series model `\link{SS_out}`.
#'
"SS_data"



#' Trace plots corresponding to the worst values of Rhat
#' @description Trace plots with kernel densities will be created for parameters with the largest (worst) associated values of `Rhat`.
#' This function is primarily intended for parameters with a vector (or array) of values.
#' @param x Output object returned from `jagsUI`
#' @param p Optional vector of parameters to subset
#' @param n.eff Whether to plot parameters with the smallest associated values of `n.eff` instead.  Defaults to `FALSE`.
#' @param margin In the case of a 2+ dimensional array associated with a given parameter, this will have the effect
#' of plotting the worst `Rhat` corresponding to each margin specified.  For example, specifying `margin=2` (column)
#' will plot the parameter with the worst `Rhat` value from each column.  In contrast, specifying `margin=NULL` (the default)
#' will cause the function to plot the single array element with the largest Rhat value.
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments or arguments to `tracedens_jags()`
#' @return `NULL`
#' @seealso \link{plotRhats}, \link{check_Rhat}
#' @author Matt Tyers
#' @examples
#' ## plotting everything
#' traceworstRhat(SS_out, parmfrow=c(3,2))
#' SS_out$Rhat  # the associated values
#'
#' traceworstRhat(SS_out, parmfrow=c(3,2), n.eff=TRUE)
#' SS_out$n.eff  # the associated values
#'
#' ## in the case of a 2-D array, setting margin=2 gives the max Rhat
#' ## associated with each column, rather than the global max
#' traceworstRhat(x=SS_out, p="cycle_s", margin=2, parmfrow=c(2,2))
#' SS_out$Rhat
#' traceworstRhat(x=SS_out, p="cycle_s", margin=2, parmfrow=c(2,2), n.eff=TRUE)
#' SS_out$n.eff
#' @export
traceworstRhat <- function(x,p=NULL,n.eff=FALSE,margin=NULL,parmfrow=NULL,...) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")

  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }
  if(!n.eff) {
    if(is.null(p)) {
      rhatlist <- x$Rhat
      p <- names(rhatlist)
    } else {
      rhatlist <- x$Rhat[names(x$Rhat) %in% p]    # do the thing like plist, that allows sig??
    }
  } else {
    if(is.null(p)) {
      rhatlist <- x$n.eff
      p <- names(rhatlist)
    } else {
      rhatlist <- x$n.eff[names(x$n.eff) %in% p]    # do the thing like plist, that allows sig??
    }
  }
  if(length(rhatlist)==0) stop("No parameters with matching names")
  for(i in 1:length(rhatlist)) {
    pp <- p[i]
    rhats <- rhatlist[[i]]
    if(!is.null(dim(rhats))) {
      if(!is.null(margin)) {
        if(max(margin)<=length(dim(rhats))) {
          if(!n.eff) maxes <- apply(rhats, MARGIN=margin, max, na.rm=T)
          if(n.eff) maxes <- apply(rhats, MARGIN=margin, min, na.rm=T)
          thenames <- NULL
          for(imax in 1:length(maxes)) {
            whichone1 <- which(rhats==maxes[imax], arr.ind=T)
            whichone <- whichone1[which.max(whichone1[,2]==imax),]
            thenames <- c(thenames, paste0(pp,"[",paste(whichone,collapse=","),"]"))
          }
        } else {
          stop("invalid margin argument")    ##### add functionality of using margin when applicable??  SS_out doesn't work
        }

      } else {
        if(length(dim(rhats))>1) {  ##### this might not work as intended for >2d
          if(!n.eff) whichones <- which(rhats==max(rhats,na.rm=T), arr.ind=T)[1,]
          if(n.eff) whichones <- which(rhats==min(rhats,na.rm=T), arr.ind=T)[1,]
        } else {
          if(!n.eff) whichones <- which(rhats==max(rhats,na.rm=T), arr.ind=T)[1]
          if(n.eff) whichones <- which(rhats==min(rhats,na.rm=T), arr.ind=T)[1]
        }
        thenames <- paste0(pp,"[",paste(whichones,collapse=","),"]")
      }
    } else {
      thenames <- pp
    }
    tracedens_jags(x,p=thenames,...=...)
  }
  # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' Random Colors
#' @description Creates a vector of randomly-generated colors.
#' @param n Vector length
#' @return A vector of colors
#' @author Matt Tyers
#' @examples
#' n <- 1000
#' cols <- rcolors(n)
#' x <- runif(n)
#' y <- runif(n)
#' plot(x,y, col=cols, pch=16)
#' @export
rcolors <- function(n) {
  bparm <- .8
  hi <- .9
  lo <- .1
  rr <- rbeta(n,bparm,bparm)*(hi-lo)+lo#runif(n=n, min=0.3, max=0.8)
  gg <- rbeta(n,bparm,bparm)*(hi-lo)+lo#runif(n=n, min=0.3, max=0.8)
  bb <- rbeta(n,bparm,bparm)*(hi-lo)+lo#runif(n=n, min=0.3, max=0.8)
  cols <- rgb(red=rr, green=gg, blue=bb)
  return(cols)
}


#' Plotting all Rhat values
#' @description Plotting all values of `Rhat` (or alternately `n.eff`) from an output object returned by `jagsUI`, or perhaps a subset of parameters.
#' This function is intended as a quick graphical check of which parameters have adequately converged.
#' @param x Output object returned from `jagsUI`
#' @param p Optional vector of parameters to subset
#' @param n.eff Optionally, whether to plot `n.eff` instead of `Rhat`.  Defaults to `FALSE`.
#' @param fence Value of horizontal lines to overlay as reference.  Accepting the default value (`NULL`) will give `fence` values
#' of 1.1 (a commonly used value) and 1.01 for Rhat, or 100 and 500 for `n.eff.`
#' @param plotsequence Whether to plot parameter vectors (or matrices) in a sequence, running left to right, which may
#' be useful for time series models, etc.  If the default (`FALSE`) is used, a vertical cluster will be plotted
#' for each parameter, resulting in a simpler plot if there are many parameters.  Note that the `Rhat` values will still be
#' plotted in sequence if the default (`FALSE`) is used.
#' @param splitarr Whether to split 2+ dimensional parameter arrays by a given dimension, rather than plotting the full
#' array in one vertical cluster or continuous sequence.  This may be recommended in the case of large arrays.  Defaults to `FALSE`.
#' @param margin If `splitarr=` is set to `TRUE`, which array margin to split by.  In the case of a 2-dimensional array, setting
#' `margin=2` will separate the array by column.  If the default (`NULL`) is accepted, the function will split by the smallest dimension,
#' therefore splitting into the fewest groups.
#' @return `NULL`
#' @seealso \link{traceworstRhat}, \link{check_Rhat}
#' @param ... additional plotting arguments
#' @author Matt Tyers
#' @examples
#' ## plotting everything
#' plotRhats(SS_out)
#' str(SS_out$Rhat)  # the associated values
#'
#' plotRhats(SS_out, n.eff=TRUE)
#' str(SS_out$n.eff)  # the associated values
#'
#' ## behavior of splitarr and margin are shown
#' plotRhats(SS_out)
#' plotRhats(SS_out, splitarr=TRUE)
#' str(SS_out$Rhat) # the associated values
#'
#' ## plotsequence may be useful in the case of a sequence of values
#' plotRhats(SS_out, p=c("trend", "cycle_s"), splitarr=TRUE, plotsequence=TRUE)
#' @export
plotRhats <- function(x,
                      p=NULL,
                      n.eff=FALSE,
                      fence=NULL,
                      plotsequence=FALSE,
                      splitarr=FALSE,
                      margin=NULL,
                      ...) {   # sepbymargin??
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")

  ## make a plotneffs( )

  if(!n.eff) {
    main <- "R hat"
    if(is.null(fence)) fence <- c(1.01,1.1)
    if(is.null(p)) {
      rhats <- x$Rhat
      p <- names(x$Rhat)
    } else {
      rhats <- x$Rhat[names(x$Rhat) %in% p]
    }
  } else {
    main <- "n eff"
    if(is.null(fence)) fence <- c(500,100)
    if(is.null(p)) {
      rhats <- x$n.eff
      p <- names(x$n.eff)
    } else {
      rhats <- x$n.eff[names(x$n.eff) %in% p]
    }
  }

  if(length(rhats)==0) stop("No parameters with matching names")  ##
  ylims <- range(fence,unlist(rhats),na.rm=T)  ## do something with fence?

  if(splitarr) {
    rhats1 <- list()
    ilist <- 1
    for(irhat in 1:length(rhats)) {
      if(!is.null(dim(rhats[[irhat]]))) {
        if(length(dim(rhats[[irhat]]))>1) {
          thedims <- dim(rhats[[irhat]])
          if(is.null(margin)) {
            margin1 <- which.min(thedims)
          } else {
            margin1 <- margin
          }
          newRhats <- apply(rhats[[irhat]], margin1, as.vector)
          for(inew in 1:ncol(newRhats)) {
            rhats1[[ilist]] <- newRhats[,inew]
            suffix <- paste0("[",rep(",",margin1-1),inew,rep(",",length(thedims)-margin1),"]")
            names(rhats1)[ilist] <- paste0(names(rhats)[irhat],suffix)
            ilist <- ilist+1
          }
        } else {
          rhats1[[ilist]] <- rhats[[irhat]]
          names(rhats1)[ilist] <- names(rhats)[irhat]
          ilist <- ilist+1
        }
      } else {
        rhats1[[ilist]] <- rhats[[irhat]]
        names(rhats1)[ilist] <- names(rhats)[irhat]
        ilist <- ilist+1
      }
    }
  } else {
    rhats1 <- lapply(rhats, as.vector)
  }

  if(plotsequence) {
    lengths <- sapply(rhats1,length)
    xlims <- c(-0.1*max(lengths), max(lengths))
    # cols <- sample(adjustcolor(rainbow(length(rhats1)), red.f=.8, blue.f=.8, green.f=.8))
    cols <- rcolors(length(rhats1))
    plot(NA, xlim=xlims, ylim=ylims, main=main, ylab="", xlab="",log='y',...=...)
    legend("topleft",col=cols, legend=names(rhats1),pch=16)
    for(i in 1:length(rhats1)) {
      points(x=1:length(rhats1[[i]]), y=rhats1[[i]], col=cols[i],pch=16)
      lines(x=1:length(rhats1[[i]]), y=rhats1[[i]], col=cols[i])
    }
  } else {
    lengths <- sapply(rhats1,length)
    xlims <- c(0.7,length(rhats1)+0.3)
    # cols <- sample(adjustcolor(rainbow(length(rhats1)), red.f=.8, blue.f=.8, green.f=.8))
    cols <- rcolors(length(rhats1))
    plot(NA, xlim=xlims, ylim=ylims, main=main, ylab="", xlab="", xaxt='n',log='y',...=...)
    axis(side=1, at=1:length(names(rhats1)), labels=names(rhats1), las=2)
    for(i in 1:length(rhats1)) {
      # points(x=jitter(rep(i, length(rhats1[[i]]))), y=rhats1[[i]])
      dx <- 1/max(lengths)
      # xx <- seq(from=i-.3, to=i+.3, length.out=length(rhats1[[i]]))
      xx <- i+seq(from=-.3*dx*lengths[i], to=.3*dx*lengths[i], length.out=lengths[i])
      points(x=xx, y=rhats1[[i]], col=cols[i])
    }
  }
  if(length(fence)==2) {
    ltylty <- 3:2
  } else {
    ltylty <- 3
  }
  abline(h=fence,lty=ltylty)
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
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{comparecat}
#' @author Matt Tyers
#' @examples
#' ## This is the same output object twice, but shows functionality.
#' comparedens(x1=asdf_jags_out, x2=asdf_jags_out, p=c("a","b","sig"))
#' @export
comparedens <- function(x1,x2, p=NULL,...) {
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

  allparms <- sort(unique(c(names(parmx1),names(parmx2))))

  plot(NA,xlim=c(0,length(allparms)+1), ylim=range(parmx1,parmx2,na.rm=T), ylab="",xlab="",xaxt="n",...=...)
  axis(1,at=1:length(allparms),labels=allparms,las=2)
  abline(v=1:length(allparms),lty=3)

  for(i in 1:length(allparms)) {
    if(allparms[i] %in% names(parmx1)) {
      xxx <- density(parmx1[,which(names(parmx1)==allparms[i])])
      polygon(i-xxx$y/max(xxx$y)*.4, xxx$x, col=adjustcolor(2,alpha.f=.5), border=2)
    }
    if(allparms[i] %in% names(parmx2)) {
      xxx <- density(parmx2[,which(names(parmx2)==allparms[i])])
      polygon(i+xxx$y/max(xxx$y)*.4, xxx$x, col=adjustcolor(4,alpha.f=.5), border=4)
    }
  }
}


#' Compare Caterpillar Plots
#' @description Interleaved caterpillar plots for all parameters (or a specified subset) from a list of `jagsUI`
#' output objects or `data.frame`s.  The intent of this function is easy comparison of inferences from multiple comparable models.
#' @param x List of output objects returned from `jagsUI` or `data.frame`s
#' @param p Optional vector of parameters to subset.  All parameters with names matching the beginning of the
#' string supplied will be returned.  If the default (`NULL`) is accepted, all parameters will be plotted.
#' @param ci Credible intervals widths to plot.  Defaults to 50% and 95%.
#' @param ylim Y-axis limits for plotting
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{caterpillar}, \link{comparedens}
#' @author Matt Tyers
#' @examples
#' ## This is the same output object three times, but shows functionality.
#' comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),
#'            p=c("a","b","sig"))
#' @export
comparecat <- function(x,p=NULL,ci=c(0.5,0.95),ylim=NULL,...) {
  if(!inherits(x,"list")) stop("Input must be a (single) list of outputs from jagsUI::jags() or data.frames.")
  xdf <- list()
  for(i in 1:length(x)) {
    if(!inherits(x[[i]], "jagsUI") & !inherits(x[[i]], "data.frame")) {
      stop("Each input element must be an output from jagsUI::jags() or data.frame.")
    }
    if(inherits(x[[i]],"jagsUI")) {
      xdf[[i]] <- jags_df(x[[i]])
    } else {
      xdf[[i]] <- x[[i]]
    }
  }
  # xdf <- lapply(x, jags_df)
  parmx <- list()
  # for(i in 1:length(x)) parmx[[i]] <- cbind(pull_post(xdf[[i]],"sig"), phi=pull_post(xdf[[i]],"phi"))
  for(i in 1:length(x)) parmx[[i]] <- pull_post(xdf[[i]],p=p)

  allparms <- sort(unique(unlist(lapply(parmx,names))))

  cilo <- sort(1-ci)/2
  cihi <- 1-cilo
  lwds <- seq(from=1,by=2,length.out=length(ci))

  if(!is.null(ylim)) {
    ylims <- ylim
  } else {
    ylims <- c(min(sapply(parmx, function(x) apply(x,2,quantile,p=min(cilo))),na.rm=T),
               max(sapply(parmx, function(x) apply(x,2,quantile,p=max(cihi))),na.rm=T))
  }

  plot(NA,xlim=c(0,length(allparms)+1), ylim=ylims, ylab="",xlab="",xaxt="n",...=...)
  axis(1,at=1:length(allparms),labels=allparms,las=2)
  # abline(v=1:length(allparms),lty=3)
  for(i in 1:length(allparms)) {
    for(ii in 1:length(x)) {
      if(allparms[i] %in% names(parmx[[ii]])) {
        xplot <- i+(seq(-.3,.3,length.out=length(x)))[ii]
        vec <- parmx[[ii]][,which(names(parmx[[ii]])==allparms[i])]
        segments(x0=rep(xplot,2), y0=quantile(vec,p=cilo), y1=quantile(vec,p=cihi),
                 lwd=lwds, lend=1, col=ii+1)
        points(xplot,median(vec),pch=16,col=ii+1)
      }
    }
  }
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
#' @seealso \link{tracedens_jags}
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
  if (!is.null(parmfrow))
    par(mfrow = parmfrow1)
}



niggle <- function() print("He was the sort of painter who can paint leaves better than trees. ")



