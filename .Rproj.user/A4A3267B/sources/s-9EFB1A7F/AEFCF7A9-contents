#' Skeleton
#' @description Provides a paste-able skeleton of code with an example JAGS model and output.
#' @param NAME Name to append to JAGS model object, etc.
#' @author Matt Tyers
#' @examples
#' skeleton("asdf")
#' @export
skeleton <- function(NAME="NAME") {
cat("
library(jagsUI)

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


niter <- 10000
ncores <- 3

",NAME,"_data <- list(x=x,
                 y=y,
                 n=length(x),
                 grp=as.numeric(as.factor(grp)),
                 ngrp=length(unique(grp)))

tstart <- Sys.time()
print(tstart)
",NAME,"_jags_out <- jagsUI::jags(model.file=\"",NAME,"_jags\", data=",NAME,"_data,
                             parameters.to.save=c(\"b0\",\"b1\",\"sig\",\"a\",\"sig_a\"),
                             n.chains=ncores, parallel=T, n.iter=niter,
                             n.burnin=niter/2, n.thin=niter/2000)
print(Sys.time() - tstart)
",sep="")
}

#' Extract df
#' @description Extract df from jagsUI output.
#' @param x Output object from jagsUI::jags()
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#' @export
jags_df <- function(x) as.data.frame(as.matrix(x$samples))

#' Subset from posterior df
#' @description Extract a subset vector or df from posterior df, with parameter names starting with a given string.
#' @param x Posterior df
#' @param p String to begin posterior names
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
#' @export
pull_post <- function(x, p) x[,substr(names(x),1,nchar(p))==p]

# skeleton()
#
# n <- 30
# x <- rnorm(n)
# grp <- sample(letters[1:3],n,replace=T)
# a <- rnorm(3)
# y <- 2*x + a[as.numeric(as.factor(grp))] + rnorm(n)
#
# skeleton("asdf")
#
# library(jagsUI)
#
# library(jagsUI)
#
# cat('model {
#   for(i in 1:n) {
#     y[i] ~ dnorm(mu[i], tau)
#     mu[i] <- b0 + b1*x[i] + a[grp[i]]
#   }
#
#   for(j in 1:ngrp) {
#     a[j] ~ dnorm(0, tau_a)
#   }
#
#   tau <- pow(sig, -2)
#   sig ~ dunif(0, 10)
#   b0 ~ dnorm(0, 0.001)
#   b1 ~ dnorm(0, 0.001)
#
#   tau_a <- pow(sig_a, -2)
#   sig_a ~ dunif(0, 10)
# }', file="asdf_jags")
#
#
# niter <- 10000
# ncores <- 3
#
# asdf_data <- list(x=x,
#                   y=y,
#                   n=length(x),
#                   grp=as.numeric(as.factor(grp)),
#                   ngrp=length(unique(grp)))
#
# tstart <- Sys.time()
# print(tstart)
# asdf_jags_out <- jagsUI::jags(model.file="asdf_jags", data=asdf_data,
#                               parameters.to.save=c("b0","b1","sig","a","sig_a"),
#                               n.chains=ncores, parallel=T, n.iter=niter,
#                               n.burnin=niter/2, n.thin=niter/2000)
# print(Sys.time() - tstart)
#
# out_df <- jags_df(asdf_jags_out)
# b1 <- pull_post(out_df,"b1")
# a <- pull_post(out_df,"a")

#' Simple traceplot
#' @description Traceplot of a single parameter.
#' @param x Posterior vector
#' @param nline Number of chains
#' @param lwd Line width
#' @param main Plot title
#' @param ... additional plotting arguments
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
  n <- length(x)/nline
  cols <- adjustcolor(rainbow(nline),red.f=.9,blue.f=.9,green.f=.9,alpha.f=.6)
  plot(NA,xlim=c(0,n),ylim=range(x,na.rm=T),main=main,xlab="iteration",ylab="value", ...=...)
  for(i in 1:nline) {
    lines(1:n, x[(n*(i-1)+1):(n*i)], col=cols[i],lwd=lwd)
  }
}

#' Traceplot of each column of a df
#' @description Traceplot of each column of a posterior df.
#' @param df Posterior df
#' @param p nline Number of chains
#' @param ... additional plotting arguments or arguments to trace_line
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
#'
#' par(mfrow=c(3,1))
#' trace_df(a, nline=3)
#' @export
trace_df <- function(df, nline, ...) {
  for(i in 1:ncol(df)) {
    trace_line(df[,i],main=names(df)[i],nline=nline,...=...)
  }
}

# trace_line(b1, nline=3)
# trace_df(b1)
# trace_df(a)
