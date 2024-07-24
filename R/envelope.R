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

  ###
  df <- df[, !is.na(x)] # taking out NA in x (not sure what it would do)
  x <- x[!is.na(x)]
  ###

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
    ### vvv this is new
    plot(NA, xlim=range(x[!is.na(med)], na.rm=TRUE),ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...=...)
    if(median) lines(x[!is.na(x) & !is.na(med)], med[!is.na(x) & !is.na(med)], col=col)
  }
  else
    if(median) lines(x[!is.na(x) & !is.na(med)], med[!is.na(x) & !is.na(med)], col=col, ...=...)
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
#' @seealso \link{envelope}, \link{crossplot}
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
  if(is.null(ylim)) {
    bounds <- sapply(df, function(x) apply(x,2, quantile, p=cilim, na.rm=T))
    ylim <- range(bounds)
  }

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
