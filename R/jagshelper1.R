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







niggle <- function() print("He was the sort of painter who can paint leaves better than trees. ")
