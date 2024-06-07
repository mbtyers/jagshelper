#' Skeleton
#' @description Prints an example 'JAGS' model and associated 'jagsUI' code to
#' the console, along with code to simulate a corresponding dataset.  This is
#' intended to serve as a template that can be altered as needed by the user.
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
#' @importFrom graphics abline legend rect text
#' @importFrom stats rbeta median cor
#' @examples
#' skeleton("asdf")
#' @export
skeleton <- function(NAME="NAME") {
cat("
library(jagsUI)

# specify model, which is written to a temporary file
",NAME,"_jags <- tempfile()
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
}', file=",NAME,"_jags)


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
# ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  ",NAME,"_jags_out <- jagsUI::jags(model.file=",NAME,"_jags, data=",NAME,"_data,
                             parameters.to.save=c(\"b0\",\"b1\",\"sig\",\"a\",\"sig_a\"),
                             n.chains=ncores, parallel=T, n.iter=niter,
                             n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(",NAME,"_jags_out)
plotRhats(",NAME,"_jags_out)
traceworstRhat(",NAME,"_jags_out, parmfrow = c(3, 3))",sep="")
}

#' Extract data.frame
#' @description Extracts the posterior samples from `jagsUI` output in the form of
#' a `data.frame`.  This simpler construction has a few benefits: operations may
#' be more straightforward, and posterior objects will be smaller files and can be
#' written to an external table or .csv, etc.
#' @param x Output object from `jagsUI::jags()`
#' @param p Optional string to begin posterior names.  If `NULL` is used, all parameters will be returned.
#' @param exact Whether name must be an exact match (`TRUE`) or with initial sub-string matching only supplied characters (`FALSE`).
#' Defaults to `FALSE.`
#' @return A `data.frame` with a column associated with each parameter and a row
#' associated with each MCMC iteration.
#' @seealso \link{pull_post}
#' @author Matt Tyers
#' @import jagsUI
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#' @export
jags_df <- function(x, p=NULL, exact=FALSE) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  xdf <- as.data.frame(as.matrix(x$samples))
  if(!is.null(p)) xdf <- pull_post(xdf, p=p, exact=exact)
  return(xdf)
}


#' Subset from posterior data.frame
#' @description Extracts a subset vector or `data.frame` from a `data.frame` consisting of more columns,
#' such that column names match a name given in the `p=` argument.  This may be useful
#' in creating smaller objects consisting of MCMC samples.
#' @param x Posterior `data.frame`
#' @param p String to begin posterior names.  If `NULL` is used, all parameters will be returned.
#' @param exact Whether name must be an exact match (`TRUE`) or with initial sub-string matching only supplied characters (`FALSE`).
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
  # # x[,substr(names(x),1,nchar(p))==p]
  # if(!inherits(x,"data.frame")) stop("Input must be a data.frame")
  #
  # if(is.null(p)) {
  #   these <- rep(T, length(names(x)))
  # } else {
  #   these <- rep(F, length(names(x)))  ## used to be F
  #   for(i in 1:length(p)) {
  #     if(exact) {
  #       these[names(x)==p[i]] <- T
  #     } else {
  #       these[substr(names(x),1,nchar(p[i]))==p[i]] <- T
  #     }
  #   }
  #   if(sum(these)==0) warning("No parameters with matching names, returning empty data.frame")
  # }
  # return(x[these])
  if(!inherits(x,"data.frame")) stop("Input must be a data.frame")

  if(is.null(p)) {
    these <- rep(T, length(names(x)))
    return(x)
  } else {
    these <- rep(F, length(names(x)))  ## used to be F
    asalist <- list()
    for(i in 1:length(p)) {
      if(exact) {
        these[names(x)==p[i]] <- T
        asalist[[i]] <- x[names(x)==p[i]]
      } else {
        these[substr(names(x),1,nchar(p[i]))==p[i]] <- T
        asalist[[i]] <- x[substr(names(x),1,nchar(p[i]))==p[i]]
      }
    }
    if(sum(these)==0) warning("No parameters with matching names, returning empty data.frame")
  }
  # return(x[these])
  return(do.call(cbind, asalist))
}

# pull_post1 <- function(x, p=NULL, exact=FALSE) {
#   # x[,substr(names(x),1,nchar(p))==p]
#   if(!inherits(x,"data.frame")) stop("Input must be a data.frame")
#
#   if(is.null(p)) {
#     these <- rep(T, length(names(x)))
#   } else {
#     these <- rep(F, length(names(x)))  ## used to be F
#     asalist <- list()
#     for(i in 1:length(p)) {
#       if(exact) {
#         these[names(x)==p[i]] <- T
#         asalist[[i]] <- x[names(x)==p[i]]
#       } else {
#         these[substr(names(x),1,nchar(p[i]))==p[i]] <- T
#         asalist[[i]] <- x[substr(names(x),1,nchar(p[i]))==p[i]]
#       }
#     }
#     if(sum(these)==0) warning("No parameters with matching names, returning empty data.frame")
#   }
#   # return(x[these])
#   return(do.call(cbind, asalist))
# }
# xx <- pull_post1(out_df, p=c("a","sig","b"))
# xx <- pull_post1(out_df, p=c("a","sig","b"), exact=T)
# xx <- pull_post1(out_df, p=c("a","sig","b","bob"), exact=T)
# xx <- pull_post1(out_df, p=c("bob"), exact=T)


#' Plist
#' @description Extracts a list of matrices, one for each saved parameter node.  Each
#' list element will be all posterior samples from that parameter node, arranged in
#' a matrix with a column associated with each MCMC chain and a row for
#' each MCMC iteration.
#' @param x `jagsUI` output object
#' @param p String to subset parameter names, if a subset is desired
#' @param exact Whether `p` should be an exact match (`TRUE`) or just match the
#' beginning of the string (`FALSE`).  Defaults to `FALSE`.
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
jags_plist <- function(x, p=NULL, exact=FALSE) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")

  x_dflist <- lapply(x$samples, as.data.frame)
  x2 <- lapply(1:length(x_dflist[[1]]), function(x) sapply(x_dflist,
                                                           "[[", x))
  names(x2) <- names(x_dflist[[1]])
  these <- rep(F, length(x2))
  if (!is.null(p)) {
    for(i in 1:length(p)) {
      if(!exact) {
        these[substr(names(x2), 1, nchar(p[i])) == p[i]] <- T
      } else {
        these[names(x2) == p[i]] <- T
      }
      # these[sapply(strsplit(names(x2), split="\\["), FUN="[", 1)==p[i]] <- T   # this is a weird hack
    }
    x2 <- x2[these]
  }
  if(length(x2)==0) warning("No parameters with matching names, returning empty list")
  return(x2)
}

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

#' Number of parameters
#' @description Total number of individual parameter nodes saved in `jagsUI` output.
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
#' @description Returns a list of the numbers of parameter nodes saved in `jagsUI` output, by parameter name.
#' As a default, what is returned for each list element is a vector of the array dimensions within the JAGS model
#'  (that is, excluding the dimension associated with the number of MCMC samples for each parameter node),
#' or alternately, just the total number of parameter nodes.
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
#' @description Returns the mean number of `Rhat` values for each parameter (by each parameter)
#' that are less than a specified threshold criterion.
#'
#' `Rhat` (Gelman-Rubin Convergence Diagnostic, or Potential Scale Reduction Factor)
#' is calculated within 'JAGS', and is
#' commonly used as a measure of convergence for a given parameter node.  Values close
#' to 1 are seen as evidence of adequate convergence.
#' @param x Output object from `jagsUI::jags()`
#' @param thresh Threshold value (defaults to 1.1)
#' @return Numeric (named) giving the proportion of Rhat values below the given threshold.
#' @seealso \link{check_neff}, \link{traceworstRhat}, \link{plotRhats}, \link{qq_postpred}, \link{ts_postpred}
#' @author Matt Tyers
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from Iterative Simulation
#' Using Multiple Sequences. *Statistical Science, 7*(4), 457–472. http://www.jstor.org/stable/2246093
#' @examples
#' check_Rhat(SS_out)
#' @export
check_Rhat <- function(x, thresh=1.1) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(sapply(x$Rhat, function(x) mean(x<thresh, na.rm=TRUE)))
}


#' Quick summary of n.eff values by parameter name
#' @description Returns the mean number of `n.eff` values (by each parameter) that are greater than a specified threshold criterion.
#'
#' `n.eff` is calculated within 'JAGS', and may be interpreted as a crude measure of
#' effective sample size for a given parameter node.
#' @param x Output object from `jagsUI::jags()`
#' @param thresh Threshold value (defaults to 500)
#' @return Numeric (named) giving the proportion of `n.eff` values above the given threshold.
#' @seealso \link{check_Rhat}, \link{traceworstRhat}, \link{plotRhats}, \link{qq_postpred}, \link{ts_postpred}
#' @author Matt Tyers
#' @examples
#' check_neff(SS_out)
#' @export
check_neff <- function(x, thresh=500) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(sapply(x$n.eff, function(x) mean(x>thresh, na.rm=TRUE)))
}


#' Logit
#' @description Logit log(x/(1-x))
#' @param x Numeric vector
#' @return Numeric vector
#' @seealso \link{expit}
#' @author Matt Tyers
#' @examples
#' logit(0.5)
#' @export
logit <- function(x) log(x/(1-x))


#' Expit, or inverse logit
#' @description Inverse logit, where logit is defined as log(x/(1-x)).
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


#' Envelope plot
#' @description Envelope plot of the posterior densities of a vector of parameter nodes,
#' in which the sequential order of nodes is important, such as a time series.
#'
#' This produces a plot of overlayed shaded strips, each corresponding to a given
#' interval width (defaults to 50 percent and 95 percent), with an overlayed
#' median line.
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
#' @param transform Should the y-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transform="exp"`is used, consider
#' adding additional plotting argument `log="y"`.
#' @param ... additional plotting arguments or arguments to `lines()`
#' @return `NULL`
#' @seealso \link{overlayenvelope}, \link{caterpillar}
#' @author Matt Tyers
#' @examples
#' ## usage with input data.frame
#' trend <- jags_df(SS_out, p="trend")
#' envelope(trend, x=SS_data$x)
#'
#' ## usage with jagsUI object
#' envelope(SS_out, p="trend")
#'
#' ## usage with 2-d jagsUI object
#' envelope(SS_out, p="cycle_s", column=1, main="cycle")
#' envelope(SS_out, p="cycle_s", column=2, col=2, add=TRUE)  ## overlay
#'
#' ## scale transformation
#' envelope(SS_out, p="trend", transform="exp", ylab="exp transform")
#' envelope(SS_out, p="trend", transform="exp", ylab="exp transform", log="y")
#' @export
envelope <- function(df,
                     p=NULL,
                     x=NA,
                     row=NULL, column=NULL,
                     median=TRUE,
                     ci=c(0.5,0.95),
                     col=4, add=FALSE, dark=.3, outline=FALSE,
                     xlab="", ylab="", main=NULL,
                     ylim=NULL,
                     transform=c("none", "exp", "expit"), ...) {


  # if(!inherits(df,"jagsUI") & !(inherits(df,"data.frame") | inherits(df,"matrix"))) {
  #   stop("Input must be a data.frame or output from jagsUI::jags() plus parameter name")
  # }
  # if(inherits(df,"jagsUI") & length(p)!=1) stop("Need single parameter name in p= argument") ###
  #
  # if(inherits(df,"jagsUI") & !is.null(p)) {
  # # if(class(df)[1]=="jagsUI" & !is.null(p)) {
  #   simslist <- df$sims.list
  #   if(all(names(simslist)!=p)) stop("No parameters with matching names") ###
  #   theparm <- simslist[names(simslist)==p][[1]]
  #   if(length(dim(theparm))==2) {
  #     df <- theparm
  #   } else {
  #     if(!is.null(row)) df <- theparm[,row,]
  #     if(!is.null(column)) df <- theparm[,,column]
  #   }
  #   if(is.null(main)) {
  #     if(is.null(row) & is.null(column)) {
  #       main <- p
  #     } else {
  #       if(!is.null(row)) main <- paste0(p,"[",row,",]")
  #       if(!is.null(column)) main <- paste0(p,"[,",column,"]")
  #     }
  #   }
  # }
  # if(is.null(main)) main <- ""
  #
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

  if(all(is.na(x))) x <- 1:ncol(df)
  df <- df[, order(x)]  # reordering for plotting if x is not in order
  x <- x[order(x)]

  transform <- match.arg(transform)
  if(transform == "exp") df <- exp(df)
  if(transform == "expit") df <- expit(df)

  ci <- sort(ci)
  loq <- apply(df, 2, quantile, p=(1-ci)/2, na.rm=T)
  hiq <- apply(df, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  if(length(ci)==1) {
    loq <- t(as.matrix(loq))
    hiq <- t(as.matrix(hiq))
  }
  med <- apply(df, 2, median, na.rm=T)

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
    # whichnotna <- unname(apply(loq,2,function(x) !all(is.na(x))))
    # firsts <- lasts <- NA
    # ilist <- whichnotna[1]
    # if(whichnotna[1]) firsts[1] <- 1
    # for(ii in 2:(length(whichnotna)-1)) {
    #   if(whichnotna[ii] & !whichnotna[ii-1]) {
    #     ilist <- ilist+1
    #     firsts[ilist] <- ii
    #   }
    #   if(whichnotna[ii] & !whichnotna[ii+1]) {
    #     lasts[ilist] <- ii
    #   }
    # }
    # if(whichnotna[length(whichnotna)]) lasts[ilist] <- length(whichnotna)
    #
    # for(ii in 1:length(firsts)) {
    #   for(i in 1:length(ci)) {
    #     polygon(c(x[firsts[ii]:lasts[ii]],rev(x[firsts[ii]:lasts[ii]])),
    #             c(loq[i,firsts[ii]:lasts[ii]],rev(hiq[i,firsts[ii]:lasts[ii]])),
    #             col=adjustcolor(col,alpha.f=dark), border=NA)
    #   }
    # }

    whichnotna <- which(unname(apply(loq,2,function(x) !all(is.na(x)))))
    # firsts <- lasts <- NA
    # ilist <- whichnotna[1]
    # if(whichnotna[1]) firsts[1] <- 1
    # for(ii in 2:(length(whichnotna)-1)) {
    #   if(whichnotna[ii] & !whichnotna[ii-1]) {
    #     ilist <- ilist+1
    #     firsts[ilist] <- ii
    #   }
    #   if(whichnotna[ii] & !whichnotna[ii+1]) {
    #     lasts[ilist] <- ii
    #   }
    # }
    # if(whichnotna[length(whichnotna)]) lasts[ilist] <- length(whichnotna)

    # for(ii in 1:length(firsts)) {
      for(i in 1:length(ci)) {
        polygon(c(x[whichnotna],rev(x[whichnotna])),
                c(loq[i,whichnotna],rev(hiq[i,whichnotna])),
                col=adjustcolor(col,alpha.f=dark), border=NA)
      }
    # }
  }
}


#' Overlay envelope plots
#' @description Overlays multiple envelope plots of posterior `data.frames`, or outputs returned from `jagsUI`.
#' This would be best suited to a set of posterior `data.frames` or 2-d matrices representing sequential vectors of parameter nodes.
#'
#' Here a single \link{envelope} plot is defined as a set of overlayed shaded strips, each corresponding to a given
#' interval width (defaults to 50 percent and 95 percent), with an overlayed
#' median line.
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
#' @param legendpos Position for optional legend.  Defaults to `"topleft"`.
#' @param ... additional plotting arguments or arguments to `lines()`
#' @param transform Should the y-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transform="exp"`is used, consider
#' adding additional plotting argument `log="y"`.
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
#'
#' ## exponential transformation
#' overlayenvelope(df=SS_out, p="cycle_s", transform="exp",
#'                 ylab="exp transform")
#' overlayenvelope(df=SS_out, p="cycle_s", transform="exp",
#'                 ylab="exp transform", log="y")
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
                            legend=TRUE, legendnames=NULL, legendpos="topleft",
                            transform=c("none", "exp", "expit"), ...) {
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

  transform <- match.arg(transform)
  if(transform == "exp") ylim <- exp(ylim)
  if(transform == "expit") ylim <- expit(ylim)

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
           median=median, dark=dark, outline=outline,xlab=xlab,ylab=ylab,
           transform=transform, ...=...)
  for(i in 2:length(df)) {
    envelope(df[[i]], ci=ci, col=cols[i], add=T, x=x, median=median, dark=dark,
             outline=outline,
             transform=transform)#, ...=...)
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
    legend(legendpos,legend=legendnames, fill=adjustcolor(cols, alpha.f=.5), border=cols)
  }
}



#' Caterpillar plot
#' @description Caterpillar plot of the posterior densities of a vector of parameter nodes,
#' in which the sequential order of nodes might not be important, such as vector of random effects.
#'
#' This produces a set of overlayed interval bars (default values are 50 percent and 95 percent),
#' with overlayed median markings, for each of a vector of parameter nodes.
#' @param df Output object returned from `jagsUI::jags()`; or alternately,
#' two-dimensional `data.frame` or matrix in which parameter node element is
#' given by column and MCMC iteration is given by row.  A vector may also be used,
#' that expresses MCMC iterations of a single parameter node.
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
#' @param ylim Y-axis limits.  If the default (`NULL`) is accepted, the limits will be determined automatically.
#' @param xax Vector of possible x-axis tick labels.  Defaults to the `data.frame` column names.
#' @param transform Should the y-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transform="exp"`is used, consider
#' adding additional plotting argument `log="y"`.
#' @param medlwd Line width of median line
#' @param medwd Relative width of median line.  Defaults to 1, perhaps smaller numbers will look better?
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{envelope}
#' @author Matt Tyers
#' @examples
#' ## usage with input data.frame
#' a <- jags_df(asdf_jags_out, p="a")
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
#'
#' ## usage with an exponential transformation
#' caterpillar(SS_out, p="trend", transform="exp", ylab="exp transform")
#' caterpillar(SS_out, p="trend", transform="exp", ylab="exp transform", log="y")
#' caterpillar(SS_out, p="trend", transform="expit", ylab="expit (inv logit) transform")
#' @export
caterpillar <- function(df,
                        p=NULL,
                        x=NA,
                        row=NULL, column=NULL,
                        median=TRUE, mean=FALSE,
                        ci=c(0.5,0.95),
                        lwd=1, col=4, add=FALSE,
                        xlab="", ylab="", main=NULL, ylim=NULL,
                        xax=NA,
                        transform=c("none", "exp", "expit"),
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

  if(!inherits(df,"jagsUI") & !inherits(df,c("matrix","data.frame","numeric","integer"))) {  ## put matrix in here
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

  df <- as.matrix(df)  ################

  transform <- match.arg(transform)
  if(transform == "exp") df <- exp(df)
  if(transform == "expit") df <- expit(df)

  ci <- rev(sort(ci))
  loq <- apply(df, 2, quantile, p=(1-ci)/2, na.rm=T)
  hiq <- apply(df, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  if(length(ci)==1) {
    loq <- t(as.matrix(loq))
    hiq <- t(as.matrix(hiq))
  }
  med <- apply(df, 2, median, na.rm=T)
  if(all(is.na(x))) x <- 1:ncol(df)
  d <- ifelse(length(x)>1, diff(x[1:2]), 1)  #########
  nn <- ncol(df)
  if(all(is.na(xax))) xax<-names(df)
  lwds <- (1+2*(1:length(ci)-1))*lwd
  if(!add) {
    if(is.null(ylim)) ylim <- range(loq,hiq,na.rm=T)
    plot(NA, type='l', xlim=range(x-(.2*d),x+(.2*d)), xlab=xlab, ylab=ylab, main=main, ylim=ylim, xaxt="n", ...=...)
    axis(1,x,labels=xax, las=list(...)$las)
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

#' Example data: asdf prior jags out
#'
#' A simple model, equivalent to that produced by the output produced by `\link{skeleton}`,
#' with the addition of prior samples for all parameters.
#'
"asdf_prior_jags_out"

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
#'
#' `Rhat` (Gelman-Rubin Convergence Diagnostic, or Potential Scale Reduction Factor)
#' is calculated within 'JAGS', and is
#' commonly used as a measure of convergence for a given parameter node.  Values close
#' to 1 are seen as evidence of adequate convergence.  `n.eff` is also calculated within 'JAGS', and may be interpreted as a crude measure of
#' effective sample size for a given parameter node.
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
#' @seealso \link{plotRhats}, \link{check_Rhat}, \link{qq_postpred}, \link{ts_postpred}
#' @author Matt Tyers
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from Iterative Simulation
#' Using Multiple Sequences. *Statistical Science, 7*(4), 457–472. http://www.jstor.org/stable/2246093
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
    if(all(is.na(rhats))) rhats[1] <- 1  # this is a weird hack to eliminate NA bug
    if(!is.null(dim(rhats))) {
      if(!is.null(margin)) {
        domargin <- (max(margin)<=length(dim(rhats)))
      } else {
        domargin <- F
      }
      if(domargin) {
        if(max(margin)<=length(dim(rhats))) {
          if(!n.eff) maxes <- apply(rhats, MARGIN=margin, max, na.rm=T)
          if(n.eff) maxes <- apply(rhats, MARGIN=margin, min, na.rm=T)
          thenames <- NULL
          for(imax in 1:length(maxes)) {
            whichone1 <- which(rhats==maxes[imax], arr.ind=T)
            whichone <- whichone1[which.max(whichone1[,2]==imax),]
            thenames <- c(thenames, paste0(pp,"[",paste(whichone,collapse=","),"]"))
          }
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
    tracedens_jags(x,p=thenames, exact=T)#,...=...)
  }

  # if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  #
  # if(!is.null(parmfrow)) {
  #   parmfrow1 <- par("mfrow")
  #   par(mfrow=parmfrow)
  #   on.exit(par(mfrow=parmfrow1))
  # }
  # if(!n.eff) {
  #   if(is.null(p)) {
  #     rhatlist <- x$Rhat
  #     p <- names(rhatlist)
  #   } else {
  #     rhatlist <- x$Rhat[names(x$Rhat) %in% p]    # do the thing like plist, that allows sig??
  #   }
  # } else {
  #   if(is.null(p)) {
  #     rhatlist <- x$n.eff
  #     p <- names(rhatlist)
  #   } else {
  #     rhatlist <- x$n.eff[names(x$n.eff) %in% p]    # do the thing like plist, that allows sig??
  #   }
  # }
  # if(length(rhatlist)==0) stop("No parameters with matching names")
  # for(i in 1:length(rhatlist)) {
  #   pp <- p[i]
  #   rhats <- rhatlist[[i]]
  #   if(!is.null(dim(rhats))) {
  #     if(!is.null(margin)) {
  #       if(max(margin)<=length(dim(rhats))) {
  #         if(!n.eff) maxes <- apply(rhats, MARGIN=margin, max, na.rm=T)
  #         if(n.eff) maxes <- apply(rhats, MARGIN=margin, min, na.rm=T)
  #         thenames <- NULL
  #         for(imax in 1:length(maxes)) {
  #           whichone1 <- which(rhats==maxes[imax], arr.ind=T)
  #           whichone <- whichone1[which.max(whichone1[,2]==imax),]
  #           thenames <- c(thenames, paste0(pp,"[",paste(whichone,collapse=","),"]"))
  #         }
  #       } else {
  #         stop("invalid margin argument")    ##### add functionality of using margin when applicable??  SS_out doesn't work
  #       }
  #
  #     } else {
  #       if(length(dim(rhats))>1) {  ##### this might not work as intended for >2d
  #         if(!n.eff) whichones <- which(rhats==max(rhats,na.rm=T), arr.ind=T)[1,]
  #         if(n.eff) whichones <- which(rhats==min(rhats,na.rm=T), arr.ind=T)[1,]
  #       } else {
  #         if(!n.eff) whichones <- which(rhats==max(rhats,na.rm=T), arr.ind=T)[1]
  #         if(n.eff) whichones <- which(rhats==min(rhats,na.rm=T), arr.ind=T)[1]
  #       }
  #       thenames <- paste0(pp,"[",paste(whichones,collapse=","),"]")
  #     }
  #   } else {
  #     thenames <- pp
  #   }
  #   tracedens_jags(x,p=thenames,...=...)
  # }
  # # if(!is.null(parmfrow)) par(mfrow=parmfrow1)
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
#'
#' `Rhat` (Gelman-Rubin Convergence Diagnostic, or Potential Scale Reduction Factor)
#' is calculated within 'JAGS', and is
#' commonly used as a measure of convergence for a given parameter node.  Values close
#' to 1 are seen as evidence of adequate convergence.  `n.eff` is also calculated within 'JAGS', and may be interpreted as a crude measure of
#' effective sample size for a given parameter node.
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
#' @seealso \link{traceworstRhat}, \link{check_Rhat}, \link{qq_postpred}, \link{ts_postpred}
#' @param ... additional plotting arguments
#' @author Matt Tyers
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from Iterative Simulation
#' Using Multiple Sequences. *Statistical Science, 7*(4), 457–472. http://www.jstor.org/stable/2246093
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


#' Compare Caterpillar Plots
#' @description Interleaved caterpillar plots for all parameters (or a specified subset) from a list of `jagsUI`
#' output objects or `data.frame`s.  The intent of this function is easy comparison of inferences from multiple comparable models.
#'
#' Here a \link{caterpillar} plot is defined as a set of overlayed interval bars (default values are 50 percent and 95 percent),
#' with overlayed median markings, for each of a vector of parameter nodes.
#' @param x List of output objects returned from `jagsUI` or `data.frame`s
#' @param p Optional vector of parameters to subset.  All parameters with names matching the beginning of the
#' string supplied will be returned.  If the default (`NULL`) is accepted, all parameters will be plotted.
#' @param ci Credible intervals widths to plot.  Defaults to 50% and 95%.
#' @param ylim Y-axis limits for plotting.  If the default (`NULL`) is accepted, limits will be automatically determined.
#' @param col Vector of colors for plotting.  If the default (`NULL`) is accepted, colors will be automatically drawn.
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param transform Should the y-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transform="exp"`is used, consider
#' adding additional plotting argument `log="y"`.
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{caterpillar}, \link{comparedens}, \link{comparepriors}
#' @author Matt Tyers
#' @examples
#' ## This is the same output object three times, but shows functionality.
#' comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),
#'            p=c("a","b","sig"))
#'
#' ## Transformed
#' comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),
#'            p=c("sig"), transform="exp")
#' comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),
#'            p=c("sig"), transform="exp", log="y")
#' @export
comparecat <- function(x, p=NULL, ci=c(0.5,0.95), ylim=NULL, col=NULL, xlab="", ylab="",
                       transform=c("none", "exp", "expit"), ...) {
  if(!inherits(x,"list")) stop("Input must be a (single) list of outputs from jagsUI::jags() or data.frames.")
  xdf <- list()
  for(i in 1:length(x)) {
    if(!inherits(x[[i]], "jagsUI") & !inherits(x[[i]], "data.frame") & !inherits(x[[i]], "matrix")) {
      stop("Each input element must be an output from jagsUI::jags() or data.frame.")
    }
    if(inherits(x[[i]],"jagsUI")) {
      xdf[[i]] <- jags_df(x[[i]])
    }
    if(inherits(x[[i]],"matrix")) {
      xdf[[i]] <- as.data.frame(x[[i]])
    }
    if(inherits(x[[i]],"data.frame")) {
      xdf[[i]] <- x[[i]]
    }
  }
  # xdf <- lapply(x, jags_df)
  parmx <- list()
  # for(i in 1:length(x)) parmx[[i]] <- cbind(pull_post(xdf[[i]],"sig"), phi=pull_post(xdf[[i]],"phi"))

  transform <- match.arg(transform)
  for(i in 1:length(x)) {
    parmx[[i]] <- pull_post(xdf[[i]], p=p)
    # parmx[[i]] <- as.matrix(parmx[[i]])

    if(transform == "exp") {
      for(j in seq_along(parmx[[i]])) {    # hopefully seq_along() goes along columns
        parmx[[i]][,j] <- exp(parmx[[i]][,j])
      }
    }
    if(transform == "expit") {
      for(j in seq_along(parmx[[i]])) {    # hopefully seq_along() goes along columns
        parmx[[i]][,j] <- expit(parmx[[i]][,j])
      }
    }
  }

  # # allparms <- sort(unique(unlist(lapply(parmx,names))))
  # allparms <- unique(unlist(lapply(parmx,names)))
  allparms <- unique(unlist(lapply(parmx, colnames)))

  cilo <- sort(1-ci)/2
  cihi <- 1-cilo
  lwds <- seq(from=1,by=2,length.out=length(ci))

  if(!is.null(ylim)) {
    ylims <- ylim
  } else {
    ylims <- c(min(unlist(sapply(parmx, function(x) apply(x,2,quantile,p=min(cilo),na.rm=T)))),
               max(unlist(sapply(parmx, function(x) apply(x,2,quantile,p=max(cihi),na.rm=T)))))
  }

  plot(NA,xlim=c(0,length(allparms)+1), ylim=ylims, ylab=ylab, xlab=xlab, xaxt="n", ...=...)
  axis(1,at=1:length(allparms),labels=allparms,las=2)
  # abline(v=1:length(allparms),lty=3)

  if(is.null(col)) {
    cols <- c(4,2,3,rcolors(100))[1:length(x)]
  } else {
    cols <- rep(col, 100)
  }

  for(i in 1:length(allparms)) {
    for(ii in 1:length(x)) {
      if(allparms[i] %in% names(parmx[[ii]])) {
        xplot <- i+(seq(-.3,.3,length.out=length(x)))[ii]
        vec <- parmx[[ii]][,which(names(parmx[[ii]])==allparms[i])]
        segments(x0=rep(xplot,2), y0=quantile(vec,p=cilo), y1=quantile(vec,p=cihi),
                 lwd=lwds, lend=1, col=cols[ii])  #col=ii+1
        points(xplot,median(vec),pch=16,col=cols[ii])  #col=ii+1
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



niggle <- function() print("He was the sort of painter who can paint leaves better than trees. ")


#' Correlation matrix from a JAGS object
#' @description Computes a correlation matrix of all MCMC samples from an object
#' returned by 'jagsUI', or an optional subset of parameter nodes.
#' @param x Output object returned from `jagsUI`
#' @param p Optional string to begin posterior names.  If `NULL` is used, all parameters will be used
#' @param exact Whether name must be an exact match (`TRUE`) or with initial sub-string matching only supplied characters (`FALSE`).
#' Defaults to `FALSE.`
#' @return A 2-dimensional correlation matrix (n X n, where n is the number of parameter nodes)
#' @seealso \link{plotcor_jags}
#' @author Matt Tyers
#' @examples
#' cor_jags(asdf_jags_out)
#' @export
cor_jags <- function(x, p=NULL, exact=FALSE) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  df <- jags_df(x, p=p, exact=exact)   ### make this jags_df
  suppressWarnings(dfcor <- cor(df))
  return(dfcor)
}


#' Plot a correlation matrix from a JAGS object
#' @description Plots a correlation matrix of all MCMC samples from an object
#' returned by 'jagsUI', or an optional subset of parameter nodes.  Correlation is
#' plotted as shades of red (positive) or blue (negative).
#'
#' In the case of vectors or arrays of nodes for each parameter name, a single axis
#' tick will be used for all nodes with a single name.  This has the effect of
#' giving greater visual weight to single parameters, and reducing plot clutter.
#'
#' Values of correlation are overlayed for all parameters with few nodes, with
#' character size scaled according to the absolute correlation.
#' @param x Output object returned from `jagsUI`
#' @param p Optional string to begin posterior names.  If `NULL` is used, all parameters will be used
#' @param exact Whether name must be an exact match (`TRUE`) or with initial sub-string matching only supplied characters (`FALSE`).
#' Defaults to `FALSE.`
#' @param mincor Minimum (absolute) correlation to use for text labels.  Defaults to 0 (all will be plotted)
#' @param maxn Maximum number of nodes per parameter name for text labels, to prevent plot clutter.  Defaults to 4.
#' @param maxcex Maximum character expansion factor for text labels.  Defaults to 1.
#' @param legend Whether to produce a plot legend.  Defaults to `TRUE`.
#' @param ... Optional plotting arguments
#' @return `NULL`
#' @seealso \link{plotcor_jags}
#' @author Matt Tyers
#' @examples
#' plotcor_jags(asdf_jags_out, maxcex=0.7)
#'
#' plotcor_jags(SS_out, p=c("trend","rate","sig"))
#' @export
plotcor_jags <- function(x, p=NULL, exact=FALSE, mincor=0, maxn=4, maxcex=1, legend=TRUE, ...) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  dfcor <- cor_jags(x=x, p=p, exact=exact)
  if(all(dim(dfcor)==0)) stop("No parameters with matching names")

  dfnames <- dimnames(dfcor)[[1]] #names(df)
  dfwhich <- sapply(strsplit(dfnames,split="[",fixed=T),FUN="[",1)
  dfhowmany <- rep(NA,length(dfnames))
  for(i in 1:length(dfnames)) dfhowmany[i] <- sum(dfwhich==dfwhich[i])
  dfdim <- cumsum(1/dfhowmany)
  dfdim1 <- c(0, dfdim[-length(dfdim)])

  xmat <- matrix(dfdim, nrow=length(dfdim), ncol=length(dfdim))
  ymat <- matrix(dfdim, nrow=length(dfdim), ncol=length(dfdim), byrow=T)
  xmat1 <- matrix(dfdim1, nrow=length(dfdim), ncol=length(dfdim))
  ymat1 <- matrix(dfdim1, nrow=length(dfdim), ncol=length(dfdim), byrow=T)

  cols <- 0*xmat
  for(i in 1:nrow(xmat)) {
    for(j in 1:i) {
      if(!is.na(dfcor[i,j])) {
        if(dfcor[i,j] > 0) {
          cols[i,j] <- cols[j,i] <- adjustcolor(2,alpha.f=dfcor[i,j])
        }
        if(dfcor[i,j] < 0) {
          cols[i,j] <- cols[j,i] <- adjustcolor(4,alpha.f=-dfcor[i,j])
        }
      }
      if(i==j) cols[i,j] <- 1
    }
  }

  plot(NA, xlim=(1+.1*legend)*range(0,dfdim), ylim=rev(range(0,dfdim)),
       yaxt="n", xaxt="n", ylab="", xlab="", bty='n',...=...)#,yaxs="i", xaxs="i"
  dfwhichunique <- unique(dfwhich)
  axis(side=1, at=1:length(dfwhichunique)-.5, labels=dfwhichunique, las=2)
  axis(side=2, at=1:length(dfwhichunique)-.5, labels=dfwhichunique, las=2)

  rect(xleft=xmat1, xright=xmat, ybottom=ymat1, ytop=ymat, border=cols, col=cols)
  for(i in 1:nrow(xmat)) {
    for(j in 1:i) {
      if((dfhowmany[i]<=maxn) & (dfhowmany[j]<=maxn) & (abs(dfcor[i,j])>=mincor)) {
        text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
        text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
      }
    }
  }
  # abline(v=0:length(dfwhichunique))
  segments(x0=rep(0, length(dfwhichunique)+1),
           x1=rep(length(dfwhichunique), length(dfwhichunique)+1),
           y0=0:length(dfwhichunique))
  segments(y0=rep(0, length(dfwhichunique)+1),
           y1=rep(length(dfwhichunique), length(dfwhichunique)+1),
           x0=0:length(dfwhichunique))

  if(legend) {
    legendby <- .25
    legendn <- 2/legendby+1

    legendl <- rep(1.05*max(dfdim),legendn)
    legendr <- rep(1.1*max(dfdim),legendn)
    legendb <- seq(from=.5*max(dfdim), to=0, length.out=legendn+1)[-legendn-1]
    legendt <- seq(from=.5*max(dfdim), to=0, length.out=legendn+1)[-1]

    legendcols <- rep(0, legendn)
    legendcors <- seq(-1,1,length.out=legendn)
    for(i in 1:(1/legendby)) {
      legendcols[i] <- adjustcolor(4, alpha.f=-legendcors[i])
      legendcols[legendn+1-i] <- adjustcolor(2, alpha.f=-legendcors[i])
    }

    rect(xleft=legendl, xright=legendr, ytop=legendt, ybottom=legendb, col=legendcols, border=NA)
    text(x=.5*(legendl+legendr), y=.5*(legendt+legendb), labels=legendcors, cex=.7)
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


#' Quantile-quantile plot from posterior predictive distribution
#' @description Produces a quantile-quantile plot, calculated from the quantiles of
#' a vector of data (most likely a time series), with respect to the matrix of associated posterior
#' predictive distributions.
#'
#' While not an omnibus posterior predictive check, this plot can be useful
#' for detecting an overparameterized model, or else improper specification
#' of observation error.  Like a traditional Q-Q plot, a well-specified model
#' will have points that lie close to the x=y line.  In the case of this
#' function, an overparametrized model will typically produce a plot with a
#' much shallower slope, possibly with many associated posterior predictive quantiles close
#' to 0.5.
#'
#' It should be noted that this function will only produce meaningful results
#' with a vector of data, as opposed to a single value.
#'
#' The posterior predictive distribution can be specified in two possible ways:
#' either a single output object from `jagsUI` with an associated parameter
#' name, or as a matrix or `data.frame` of posterior samples.
#' @param ypp Either a matrix or `data.frame` of posterior samples, or an
#' output object returned from `jagsUI` and a supplied parameter name
#' @param y The associated data vector
#' @param p A character name, if a `jagsUI` object is passed to `ypp`
#' @param add Whether to add the plot to an existing plot.  Defaults to `FALSE`.
#' @param ... Optional plotting arguments
#' @return `NULL`
#' @note This function assumes the existence of a matrix of posterior predictive
#' samples corresponding to a data vector, the construction of which must be
#' left to the user.  This can be accomplished within JAGS, or using appropriate
#' simulation from the posterior samples.
#' @seealso \link{ts_postpred}, \link{check_Rhat}, \link{check_neff}, \link{traceworstRhat}, \link{plotRhats}
#' @author Matt Tyers
#' @examples
#' # first, a quick look at the example data...
#' str(SS_data)
#' str(SS_out$sims.list$ypp)
#'
#' # plotting the example posterior predictive distribution with the data
#' # points overlayed.  Note the overdispersion in the posterior predictive.
#' caterpillar(SS_out, p="ypp")
#' points(SS_data$y)
#'
#' # using a jagsUI object as ypp input
#' qq_postpred(ypp=SS_out, p="ypp", y=SS_data$y)
#'
#' # using a matrix as ypp input
#' qq_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y)
#' @export
qq_postpred <- function(ypp, y, p=NULL, add=FALSE, ...) { # ypp is a matrix, y is a vector
  if(!inherits(ypp, c("matrix","data.frame")) & !inherits(ypp, "jagsUI")) stop("Argument ypp must be a posterior matrix or jagsUI object.")
  if(inherits(ypp, "jagsUI") & is.null(p)) stop("Parameter name must be supplied to p= argument if jagsUI object is used in argument ypp")
  if(inherits(ypp, "jagsUI") & !is.null(p)) {
    ypp <- ypp$sims.list[names(ypp$sims.list)==p][[1]]   # rework this with jags_df?
  }
  if(length(y)<=1) stop("Data (argument y) must be a vector for meaningful diagnostics")
  if(ncol(ypp) > length(y)) {
    stop("Posterior matrix ypp has more columns than length of data matrix y")
  }
  if(ncol(ypp) < length(y)) {
    warning("Posterior matrix ypp has fewer columns than length of data matrix y")
    # define ymat somehow differently
    # ymat <- matrix(NA, nrow=nrow(ypp), ncol=ncol(ypp))
    ymat <- matrix(y, nrow=nrow(ypp), ncol=length(y), byrow=T)   # actually this should work for both
    ypp1 <- matrix(NA, nrow=nrow(ypp), ncol=length(y))
    ypp1[,1:ncol(ypp)] <- ypp
    ypp <- ypp1
  }
  if(ncol(ypp) == length(y)){
    ymat <- matrix(y, nrow=nrow(ypp), ncol=ncol(ypp), byrow=T)
  }
  qpp <- sort(colMeans(ymat>=ypp))
  qtheo <- (1:length(qpp))/length(qpp)
  if(!add) {
    plot(qtheo, qpp, xlim=0:1, ylim=0:1, xlab="Theoretical quantile", ylab="Posterior Predictive quantile", ...=...)
    abline(0,1, lty=2)
  } else {
    points(qtheo, qpp, ...=...)
  }
}


#' Time series plot of centered posterior predictive distribution
#' @description Produces a plot of centered posterior predictive distributions
#' associated with a vector of data (most likely a time series),
#' defined as the difference between posterior predictive and posterior predictive
#' median.
#'
#' Also overlays the posterior predictive residuals, defined as the differences
#' between data values and their respective posterior predictive medians.
#'
#' While not an omnibus posterior predictive check, this plot can be useful
#' for detecting an overparameterized model, or else improper specification
#' of observation error.
#'
#' It should be noted that this function will only produce meaningful results
#' with a vector of data, as opposed to a single value.
#'
#' The posterior predictive distribution can be specified in two possible ways:
#' either a single output object from `jagsUI` with an associated parameter
#' name, or as a matrix or `data.frame` of posterior samples.
#' @param ypp Either a matrix or `data.frame` of posterior samples, or an
#' output object returned from `jagsUI` and a supplied parameter name
#' @param y The associated data vector
#' @param p A character name, if a `jagsUI` object is passed to `ypp`
#' @param x The time measurements associated with time series `y`.  If the default
#' `NULL` is accepted, equally-spaced integer values will be used.
#' @param lines Whether to add a line linking data time series points.  Defaults to `FALSE`.
#' @param transform Should the y-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transform="exp"`is used, consider
#' adding additional plotting argument `log="y"`.
#' @param ... Additional arguments to \link{envelope}
#' @return `NULL`
#' @note This function assumes the existence of a matrix of posterior predictive
#' samples corresponding to a data vector, the construction of which must be
#' left to the user.  This can be accomplished within JAGS, or using appropriate
#' simulation from the posterior samples.
#' @seealso \link{qq_postpred}, \link{check_Rhat}, \link{check_neff}, \link{traceworstRhat}, \link{plotRhats}
#' @author Matt Tyers
#' @examples
#' # first, a quick look at the example data...
#' str(SS_data)
#' str(SS_out$sims.list$ypp)
#'
#' # plotting the example posterior predictive distribution with the data
#' # points overlayed.  Note the overdispersion in the posterior predictive.
#' caterpillar(SS_out, p="ypp")
#' points(SS_data$y)
#'
#' # using a jagsUI object as ypp input
#' ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y)
#'
#' # using a matrix as ypp input
#' ts_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y)
#'
#' # exp transformation
#' ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y, transform="exp")
#' ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y, transform="exp", log="y")
#' @export
ts_postpred <- function(ypp, y, p=NULL, x=NULL, lines=FALSE,
                        transform=c("none", "exp", "expit"), ...) { #p=NULL  ?? style it after qq_postpred
  if(!inherits(ypp, c("matrix","data.frame")) & !inherits(ypp, "jagsUI")) stop("Argument ypp must be a posterior matrix or jagsUI object.")
  if(inherits(ypp, "jagsUI") & is.null(p)) stop("Parameter name must be supplied to p= argument if jagsUI object is used in argument ypp")
  if(inherits(ypp, "jagsUI") & !is.null(p)) {
    ypp <- ypp$sims.list[names(ypp$sims.list)==p][[1]]   # rework this with jags_df?
  }
  if(length(y)<=1) stop("Data (argument y) must be a vector for meaningful diagnostics")
  if(ncol(ypp)!=length(y)) stop("Posterior matrix ypp must have the same number of columns as length of data matrix y")
  meds <- apply(ypp, 2, median, na.rm=T)
  ypp_resid <- ypp - matrix(meds, byrow=TRUE, nrow=nrow(ypp), ncol=ncol(ypp))
  transform <- match.arg(transform)
  yplot <- y-meds

  ## finding y limits for plotting
  dots <- list(...)
  if(!is.null(dots$ci)) {
    ci <- max(dots$ci)
  } else {
    ci <- 0.95
  }
  ylim1 <- apply(ypp_resid, 2, quantile, p=(1-ci)/2)
  ylim2 <- apply(ypp_resid, 2, quantile, p=1-((1-ci)/2))

  ## transforming if needed
  if(transform == "exp") {
    yplot <- exp(yplot)
    ylim1 <- exp(ylim1)
    ylim2 <- exp(ylim2)
  }
  if(transform == "expit") {
    yplot <- expit(yplot)
    ylim1 <- expit(ylim1)
    ylim2 <- expit(ylim2)
  }

  envelope(ypp_resid, x=x, ylab="Diff from post pred median",
           transform=transform, ylim=range(ylim1, ylim2, yplot, na.rm=TRUE),
           ...=...)
  if(is.null(x)) x <- seq_along(y)
  points(x=x, y=yplot)
  if(lines) lines(x=x, y=yplot)
}


#' @export
plot_postpred <- function(ypp, y, p=NULL, x=NULL, lines=FALSE,
                          transform=c("none", "exp", "expit"), ...) {   # include a parmfrow

  if(!inherits(ypp, c("matrix","data.frame")) & !inherits(ypp, "jagsUI")) stop("Argument ypp must be a posterior matrix or jagsUI object.")
  if(inherits(ypp, "jagsUI") & is.null(p)) stop("Parameter name must be supplied to p= argument if jagsUI object is used in argument ypp")
  if(inherits(ypp, "jagsUI") & !is.null(p)) {
    ypp <- ypp$sims.list[names(ypp$sims.list)==p][[1]]   # rework this with jags_df?
  }
  if(length(y)<=1) stop("Data (argument y) must be a vector for meaningful diagnostics")
  if(ncol(ypp)!=length(y)) stop("Posterior matrix ypp must have the same number of columns as length of data matrix y")

  # plot 1
  qq_postpred(ypp=ypp, y=y)

  # plot 2
  if(is.null(x)) x <- seq_along(y)
  ts_postpred(ypp=ypp, y=y, x=x, lines=lines, transform=transform)

  # plot 3
  ymeds <- apply(ypp, 2, median, na.rm=TRUE)
  ts_postpred(ypp=ypp, y=y, x=ymeds, lines=lines, transform=transform, xlab="Post pred median")

  # plot 4
  thecat <- cut(ymeds, breaks=floor(sqrt(length(ymeds))))  ## this is a throwaway breaks=, make it smarter please

  # nperbin <- 5

  xplot <- tapply(y, thecat, mean, na.rm=TRUE)
  thesd <- tapply(y, thecat, sd, na.rm=TRUE)

  # ylims <- c(min(thesd, na.rm=TRUE)-0.5*diff(range(thesd, na.rm=TRUE)), max(thesd, na.rm=TRUE))
  ylims <- range(thesd, na.rm=TRUE)

  plot(xplot, thesd, type="b",
       xlim=range(ymeds, na.rm=TRUE), ylim=ylims,
       xlab="Post pred median", ylab="PP residual SD (binned)")

  # points(x=ymeds, y=0*ymeds+ylims[1]) # jitter this somehow

  # segments(x0=ymeds,
  #          y0=rep(0, length(ymeds)), y1=rep(ylims[1], length(ymeds)),
  #          col=as.numeric(as.factor(thecat))+1)
}
# par(mfrow=c(2,2))
# plot_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y, x=SS_data$x)

#' Compare Priors
#' @description Side-by-side kernel density plots for all parameters with parameter
#' names ending in `"_prior"`, and corresponding parameters without.  It should
#' be noted that these parameters must be specified in JAGS as well as the
#' corresponding parameters, and this is left to the user.
#'
#' This function is a wrapper of \link{comparedens}.
#'
#' Kernel densities are plotted vertically, either left- or right-facing.  Parameters with the same name are
#'  plotted facing one another.
#' @param x Output object returned from jagsUI::jags()
#' @param parmfrow Optional call to `par(mfrow)` for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional arguments to \link{comparedens}
#' @return `NULL`
#' @seealso \link{comparecat}, \link{comparedens}, \link{plotdens}
#' @author Matt Tyers
#' @examples
#' ## a look at what parameters exist in the input object
#' nbyname(asdf_prior_jags_out)
#'
#' ## then, showing the function usage
#' comparepriors(asdf_prior_jags_out, parmfrow=c(2, 3))
#' @export
comparepriors <- function(x, parmfrow=NULL,...) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")

  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
    on.exit(par(mfrow=parmfrow1))
  }

  # get names
  thenames <- names(x$sims.list)

  # find names ending in "_prior"
  thepriors <- NULL
  for(i in 1:length(thenames)) {
    thesplit <- strsplit(thenames[i], split="_")[[1]]
    if(thesplit[length(thesplit)] == "prior") thepriors <- c(thepriors, i)
  }
  if(is.null(thepriors)) warning('No parameter names ending in "_prior"')

  for(i_prior in thepriors) {
    # find a posterior with matching name
    thepriorname <- thenames[i_prior]
    thepostname <- NULL
    for(i_post in 1:length(thenames)) {
      if(thenames[i_post] == substr(thepriorname, 1, nchar(thepriorname)-6)) {
        thepostname <- thenames[i_post]
      }
    }
    if(!is.null(thepostname)) {
      priordf <- as.data.frame(x$sims.list[thepriorname])
      postdf <- as.data.frame(x$sims.list[thepostname])
      names(priordf) <- names(postdf)
      comparedens(x1=priordf, x2=postdf, main=thepostname, legendnames=c("prior","post"),...=...)
    }
  }
}
