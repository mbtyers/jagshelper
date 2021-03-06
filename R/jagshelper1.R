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


#' Plist
#' @description Extract a list of nrep X nchain matrices, one for each parameter.
#' @param x jagsUI object
#' @param p String to subset parameter names, if a subset is desired
#' @author Matt Tyers
#' @examples
#' out_plist <- jags_plist(asdf_jags_out)
#' str(out_plist)
#'
#' a_plist <- jags_plist(asdf_jags_out, p="a")
#' str(a_plist)
#' @export
jags_plist <- function(x, p=NULL) {
  x_dflist <- lapply(x$samples, as.data.frame)
  x2 <- lapply(1:length(x_dflist[[1]]), function(x) sapply(x_dflist, "[[",x))
  names(x2) <- names(x_dflist[[1]])
  if(!is.null(p)) x2 <- x2[substr(names(x2),1,nchar(p))==p]
  return(x2)
}

#' Traceplot of jagsUI object
#' @description Smarter traceplot of a whole jagsUI object, or subset.
#' @param x Posterior jagsUI object
#' @param p parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
#' @param parmfrow Optional call to par(mfrow) for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments
#' @author Matt Tyers
#' @examples
#'
#' trace_jags(asdf_jags_out, parmfrow=c(4,2))
#' trace_jags(asdf_jags_out, p="a", parmfrow=c(3,1))
#' @export
trace_jags <- function(x,p=NULL,parmfrow=NULL,lwd=1,...) {
  x_plist <- jags_plist(x,p=p)
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
  }

  nline <- ncol(x_plist[[1]])
  cols <- adjustcolor(rainbow(nline),red.f=.9,blue.f=.9,green.f=.9,alpha.f=.6)
  for(i in 1:length(x_plist)) {
    plot(NA, xlim=c(1,nrow(x_plist[[i]])), ylim=range(x_plist[[i]]), main=names(x_plist)[i], xlab="iter", ylab="", ...=...)
    for(j in 1:nline) lines(x_plist[[i]][,j], col=cols[j], lwd=lwd)
  }

  if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' By-chain kernel densities of jagsUI object
#' @description By-chain kernel densities of a whole jagsUI object, or subset.
#' @param x Posterior jagsUI object
#' @param p parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
#' @param parmfrow Optional call to par(mfrow) for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments
#' @author Matt Tyers
#' @examples
#'
#' chaindens_jags(asdf_jags_out, parmfrow=c(4,2))
#' chaindens_jags(x=asdf_jags_out, p="a", parmfrow=c(3,1))
#' @export
chaindens_jags <- function(x,p=NULL,parmfrow=NULL,lwd=1,...) {
  x_plist <- jags_plist(x,p=p)
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
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

  if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' Combination of trace plots and by-chain kernel densities of jagsUI object
#' @description Combination of trace plots and by-chain kernel densities of a whole jagsUI object, or subset.
#' @param x Posterior jagsUI object
#' @param p parameter name for subsetting: if this is specified, only parameters with names beginning with this string will be plotted.
#' @param parmfrow Optional call to par(mfrow) for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
#' @param ... additional plotting arguments
#' @author Matt Tyers
#' @examples
#'
#' tracedens_jags(asdf_jags_out, parmfrow=c(4,2))
#' tracedens_jags(x=asdf_jags_out, p="a", parmfrow=c(3,1))
#' @export
tracedens_jags <- function(x,p=NULL,parmfrow=NULL,lwd=1,...) {
  x_plist <- jags_plist(x,p=p)
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
  }

  nline <- ncol(x_plist[[1]])
  cols <- adjustcolor(rainbow(nline),red.f=.9,blue.f=.9,green.f=.9,alpha.f=.6)
  for(i in 1:length(x_plist)) {
    allbw <- density(as.vector(x_plist[[i]]))$bw
    denses <- apply(x_plist[[i]], 2, density, bw=allbw)
    dx <- sapply(denses, function(x) x$x)
    dy <- sapply(denses, function(x) x$y)
    # abline(v=nrow(x_plist[[i]]))
    dy1 <- (dy/max(dy)*0.3 + 1)*nrow(x_plist[[i]])
    plot(NA, xlim=c(1, 1.3*nrow(x_plist[[i]])), ylim=range(x_plist[[i]]), main=names(x_plist)[i], xlab="iter", ylab="", ...=...)
    for(j in 1:nline) lines(x_plist[[i]][,j], col=cols[j], lwd=lwd)
    # plot(NA, xlim=range(dx), ylim=range(dy), main=names(x_plist)[i], xlab="",ylab="",...=...)
    for(j in 1:nline) lines(dy1[,j], dx[,j], col=cols[j], lwd=lwd)
  }

  if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}

#' Number of parameters
#' @description Number of parameters saved in jagsUI output.
#' @param x Output object from jagsUI::jags()
#' @author Matt Tyers
#' @examples
#' head(jags_df(asdf_jags_out))
#'
#' nparam(asdf_jags_out)
#' @export
nparam <- function(x) ncol(jags_df(x))


#' Number of parameters by parameter name
#' @description Returns a summary of the numbers of parameters saved in jagsUI output by parameter name.
#' As a default, what is actually returned is a vector of the array dimensions (excluding the number of mcmc iterations),
#' or alternately, just the total number of parameters
#' @param x Output object from jagsUI::jags()
#' @param justtotal Whether to just report the total number of parameters, as opposed to dims.
#' @author Matt Tyers
#' @examples
#' head(jags_df(asdf_jags_out))
#'
#' nbyname(asdf_jags_out)
#'
#' nparam(ts_jags_out)
#' nbyname(ts_jags_out)
#' @export
nbyname <- function(x, justtotal=FALSE) {
  # sapply(x$sims.list, function(x) ncol(as.matrix(x)))
  out1 <- sapply(x$sims.list, function(x) dim(x)[-1])
  out1[sapply(out1, is.null)] <- 1
  if(!justtotal) return(out1)
  if(justtotal) return(sapply(out1, prod))
}


#' Quick summary of Rhat values by parameter name
#' @description Returns the mean number of Rhat values (by each parameter) that are less than a specified threshold criterion.
#' @param x Output object from jagsUI::jags()
#' @param thresh Threshold value (defaults to 1.1)
#' @author Matt Tyers
#' @examples
#' check_Rhat(ts_jags_out)
#' @export
check_Rhat <- function(x, thresh=1.1) sapply(x$Rhat, function(x) mean(x<thresh, na.rm=T))


#' Quick summary of n.eff values by parameter name
#' @description Returns the mean number of n.eff values (by each parameter) that are greater than a specified threshold criterion.
#' @param x Output object from jagsUI::jags()
#' @param thresh Threshold value (defaults to 500)
#' @author Matt Tyers
#' @examples
#' check_neff(ts_jags_out)
#' @export
check_neff <- function(x, thresh=500) sapply(x$n.eff, function(x) mean(x>thresh, na.rm=T))


#' Logit
#' @description Logit log(x/(1-x)
#' @param x A number
#' @author Matt Tyers
#' @examples
#' logit(0.5)
#' @export
logit <- function(x) log(x/(1-x))


#' Expit, or inverse logit
#' @description Inverse logit, where logit is log(x/(1-x).  Expit is exp(x)/(1+exp(x))
#' @param x A number
#' @author Matt Tyers
#' @examples
#' expit(0)
#' @export
expit <- function(x) exp(x)/(1+exp(x))




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

#' Simple by-chain kernel density plot
#' @description By-chain kernel density plot of a single parameter.
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
#' chaindens_line(b1, nline=3, main="b1")
#' @export
chaindens_line <- function(x, nline, lwd=1, main="", ...) {
  # if(is.null(nline)) nline <- length(x)/n
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

#' Traceplot of each column of a df
#' @description Traceplot of each column of a posterior df.
#' @param df Posterior df
#' @param nline Number of chains
#' @param parmfrow Optional call to par(mfrow) for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
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
#'
#' trace_df(a, nline=3, parmfrow=c(3,1))
#' @export
trace_df <- function(df, nline, parmfrow=NULL, ...) {
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
  }
  thencol <- ifelse(is.null(ncol(df)), 1, ncol(df))
  if(thencol==1) {
    trace_line(df,nline=nline,...=...)
  } else {
  for(i in 1:ncol(df)) {
    trace_line(df[,i],main=names(df)[i],nline=nline,...=...)
  }}
  if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}

#' By-chain kernel density of each column of a df
#' @description By-chain kernel density plot of each column of a posterior df.
#' @param df Posterior df
#' @param nline Number of chains
#' @param parmfrow Optional call to par(mfrow) for the number of rows & columns of plot window.  Returns the graphics device to previous state afterward.
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
#'
#' chaindens_df(a, nline=3, parmfrow=c(3,1))
#' @export
chaindens_df <- function(df, nline, parmfrow=NULL, ...) {
  if(!is.null(parmfrow)) {
    parmfrow1 <- par("mfrow")
    par(mfrow=parmfrow)
  }
  for(i in 1:ncol(df)) {
    chaindens_line(df[,i],main=names(df)[i],nline=nline,...=...)
  }
  if(!is.null(parmfrow)) par(mfrow=parmfrow1)
}


#' Envelope plot
#' @description Envelope plot of posterior df.  This would be best suited to a single posterior df representing a time series.  I've written and re-written this function many times!
#' @param df Posterior df
#' @param x Vector of X-coordinates for plotting.
#' @param median Whether to include median line
#' @param ci Vector of intervals to overlay.  Defaults to 50 percent and 95 percent.
#' @param col Color for plotting
#' @param add Whether to add to existing plot
#' @param dark Opacity (0-1) for envelopes.  Note that multiple overlapping intervals will darken the envelope.
#' @param outline Whether to just envelope outlines
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param main Plot title
#' @param ... additional plotting arguments or arguments to trace_line
#' @note Calling this function with argument ci=0 (or F) is effectively a tricksy way of producing a plot with just a median line.
#' @author Matt Tyers
#' @examples
#' out_df <- jags_df(asdf_jags_out)
#'
#' b1 <- pull_post(out_df,"b1")
#' a <- pull_post(out_df,"a")
#'
#' envelope(a)
#' envelope(a, ci=seq(.1,.9,by=.1), dark=.15, lwd=3)
#' envelope(a, ci=seq(.1,.9,by=.1), dark=.2, outline=T)
#' envelope(a, ci=F, add=T, lwd=3)
#'
#' ts_out_df <- jags_df(ts_jags_out)
#' mu <- pull_post(ts_out_df, "mu")
#' ypred <- pull_post(ts_out_df, "ypred")
#' envelope(x=2010:2020, mu)
#' envelope(x=2010:2020, ypred, outline=T, median=F, add=T, ci=.95, col=1, lty=2)
#' @export
envelope <- function(df, x=NA, median=T, ci=c(0.5,0.95), col=4, add=F, dark=.3, outline=F, xlab="", ylab="", main="", ...) {
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
    plot(NA, xlim=range(x),ylim=range(loq,hiq,na.rm=T), xlab=xlab, ylab=ylab, main=main, ...=...)
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



#' Caterpillar plot
#' @description Caterpillar plot of posterior df.  This would be best suited to a single posterior df representing random effects or group means.
#' @param df Posterior df
#' @param x Vector of X-coordinates for plotting.
#' @param median Whether to include medians
#' @param mean Whether to include means
#' @param ci Vector of intervals to overlay.  Defaults to 50 percent and 95 percent.
#' @param col Color for plotting
#' @param add Whether to add to existing plot
#' @param lwd Base line width for plotting.  Defaults to 1.
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param main Plot title
#' @param xax vector of possible x-axis tick labels.  Defaults to the df column names
#' @param medlwd Line width of median line
#' @param medwd Relative width of median line.  Defaults to 1, perhaps smaller numbers will look better?
#' @param ... additional plotting arguments or arguments to trace_line
#' @note Calling this function with argument ci=0 (or F) is effectively a tricksy way of producing a plot with just a median line.
#' @author Matt Tyers
#' @examples
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
#' ts_out_df <- jags_df(ts_jags_out)
#' mu <- pull_post(ts_out_df, "mu")
#' ypred <- pull_post(ts_out_df, "ypred")
#' caterpillar(x=2010:2020, mu, xax=2010:2020)
#' @export
caterpillar <- function(df, x=NA, median=T, mean=F, ci=c(0.5,0.95), lwd=1, col=4, add=F, xlab="", ylab="", main="", xax=NA, medlwd=lwd, medwd=1,...) {
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
