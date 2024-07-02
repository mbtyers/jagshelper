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
#' @param x Output object returned from `jagsUI`, or a data.frame with MCMC output
#' @param p Optional string to begin posterior names.  If `NULL` is used, all parameters will be used
#' @param exact Whether name must be an exact match (`TRUE`) or with initial sub-string matching only supplied characters (`FALSE`).
#' Defaults to `FALSE.`
#' @param mincor Minimum (absolute) correlation to use for text labels.  Defaults to 0 (all will be plotted)
#' @param maxn Maximum number of nodes per parameter name for text labels, to prevent plot clutter.  Defaults to 4.
#' @param maxcex Maximum character expansion factor for text labels.  Defaults to 1.
#' @param legend Whether to produce a plot legend.  Defaults to `TRUE`.
#' @param ... Optional plotting arguments
#' @return `NULL`
#' @seealso \link{cor_jags}
#' @author Matt Tyers
#' @examples
#' plotcor_jags(asdf_jags_out, maxcex=0.7)
#'
#' plotcor_jags(SS_out, p=c("trend","rate","sig"))
#' @export
plotcor_jags <- function(x, p=NULL, exact=FALSE, mincor=0, maxn=4, maxcex=1, legend=TRUE, ...) {
  if(!inherits(x, c("data.frame","matrix","jagsUI"))) {
    stop("Input must be an output object returned from jagsUI::jags(), or data.frame with MCMC output.")
  }

  if(inherits(x,"jagsUI")) {
    dfcor <- cor_jags(x=x, p=p, exact=exact)
    if(all(dim(dfcor)==0)) stop("No parameters with matching names")
  }
  if(inherits(x, c("data.frame","matrix"))) {
    dfcor <- cor(x=as.data.frame(x))
  }

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

  # rect(xleft=xmat1, xright=xmat, ybottom=ymat1, ytop=ymat, border=cols, col=cols)
  # for(i in 1:nrow(xmat)) {
  #   for(j in 1:i) {
  #     if((dfhowmany[i]<=maxn) & (dfhowmany[j]<=maxn) & (abs(dfcor[i,j])>=mincor)) {
  #       text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
  #       text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
  #     }
  #   }
  # }

  rect(xleft=xmat1, xright=xmat, ybottom=ymat1, ytop=ymat, border=cols, col=cols)
  dfcor4cex <- dfcor
  dfcor4cex[is.na(dfcor4cex)] <- 0.5
  for(i in 1:nrow(xmat)) {
    for(j in 1:i) {
      # if((dfhowmany[i]<=maxn) & (dfhowmany[j]<=maxn) & (abs(dfcor[i,j])>=mincor)) {
      if((dfhowmany[i]<=maxn) & (dfhowmany[j]<=maxn) & (abs(dfcor4cex[i,j])>=mincor)) {
        # text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
        # text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
        text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor4cex[i,j])^.3)
        text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor4cex[i,j])^.3)
        if(is.na(dfcor[i,j])) {
          text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels="NA", cex=maxcex*abs(dfcor4cex[i,j])^.3)
          text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels="NA", cex=maxcex*abs(dfcor4cex[i,j])^.3)
        }
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
  # thecat <- cut(ymeds, breaks=floor(sqrt(length(ymeds))))  ## this is a throwaway breaks=, make it smarter please
  thecat <- cut(rank(ymeds),
                breaks=seq(from=1, to=length(ymeds),
                           length.out=floor(sqrt(length(ymeds)))),
                include.lowest = TRUE)

  # nperbin <- 5

  xplot <- tapply(y, thecat, mean, na.rm=TRUE)
  thesd <- tapply(y, thecat, sd, na.rm=TRUE)

  # ylims <- c(min(thesd, na.rm=TRUE)-0.5*diff(range(thesd, na.rm=TRUE)), max(thesd, na.rm=TRUE))
  ylims <- range(c(0, thesd), na.rm=TRUE)

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

# to do:
# - xlab for plot 2
# - rug or SOMETHING for plot 4
# - maybe overlay rolling / binned SD on plots 2 & 3?
# - decide if these are even the plots!
# - add parmfrow=c(2,2)?? maybe don't supply default actually
# - get inspiration from plot.lm
# - visually differentiate plots 2 and 3?  maybe just add titles
# - test drive this with another project: maybe lake trout stuff
# -- see if I can find the bottleneck!!
# - jags_df??  also for qq_ and ts_
# - add documentation / tests / NEWS of course
# - maybe this is worthy of adding to README and Vignette
