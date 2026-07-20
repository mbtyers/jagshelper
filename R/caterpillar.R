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
#' @param na.rm Whether to exclude NA columns.  Defaults to TRUE.
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
#' @param xlim X-axis limits.  If the default (`NULL`) is accepted, the limits will be determined automatically.
#' @param ylim Y-axis limits.  If the default (`NULL`) is accepted, the limits will be determined automatically.
#' @param xlim_add Additional X-axis limit.  Defaults to `NULL`.  May be useful in having a axis start at 0.
#' @param ylim_add Additional Y-axis limit.  Defaults to `NULL`.  May be useful in having a axis start at 0.
#' @param xax Vector of possible x-axis tick labels.  Defaults to the `data.frame` column names.
#' @param xorder Whether to draw interval bars in order of posterior median.  Defaults to FALSE.
#' @param horizontal Whether to produce a horizontal plot, that is, with intervals on the x-axis.  Defaults to FALSE.
#' @param transform Should the y-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transform="exp"`is used, consider
#' adding additional plotting argument `log="y"`.
#' @param medlwd Line width of median line
#' @param medwd Relative width of median line.  Defaults to 1, and acts as a multiplicative adjustment.
#' @param padwd Relative width of padding beyond the median line.  Defaults to 1, and acts as a multiplicative adjustment.
#' @param all_ticks Whether to force plotting of all axis tick labels.  Defaults to `FALSE`.
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{envelope}, \link{crossplot}
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
#' caterpillar(a, xax=c("effect 1", "effect 2", "effect 3"),
#'             horizontal = TRUE, las=1)
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
                        na.rm=TRUE,   # think about why this would't be
                        row=NULL, column=NULL,
                        median=TRUE, mean=FALSE,
                        ci=c(0.5,0.95),
                        lwd=1, col=4, add=FALSE,
                        xlab="", ylab="", main=NULL,
                        xlim=NULL, ylim=NULL,
                        xlim_add=NULL, ylim_add=NULL,
                        xax=NA,
                        xorder = FALSE, horizontal=FALSE,
                        transform=c("none", "exp", "expit"),
                        medlwd=lwd, medwd=1, padwd=1,
                        all_ticks=FALSE, ...) {
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
  if(na.rm) {
    keep <- which(!is.na(colMeans(df)))
    # df <- as.matrix(df[, keep])  ################
    df <- df[, keep, drop=FALSE]
    if(!all(is.na(x))) {
      x <- x[keep]
    }
    if(!all(is.na(xax))) {
      xax <- xax[keep]
    }
    if(length(col) > 1) {
      col <- col[keep]
    }
  }

  ### new
  dots <- list(...)  # probably a more efficient way to do this
  xaxlog <- FALSE
  if(!is.null(dots$log)) {
    if(!horizontal & dots$log %in% c("x","xy")) {
      xaxlog <- TRUE
    }
    if(horizontal & dots$log %in% c("y","xy")) {
      xaxlog <- TRUE
    }
  }
  if(add & par("xlog")) {
    xaxlog <- TRUE
  }
  ### new

  ### new
  if(any(is.na(x)) & !all(is.na(x))) {
    # df <- df[, !is.na(x)]
    df <- df[, which(!is.na(x) & seq_along(x) %in% 1:ncol(df)),
             drop=FALSE]
    if(any(!is.na(xax))) {
      xax <- xax[!is.na(x)]
    }
    if(length(col) > 1) {
      col <- col[!is.na(x)]
    }
    x <- x[!is.na(x)]
  }
  ### new

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
  # if(all(is.na(x))) x <- 1:ncol(df)

  ## new
  xorig <- x
  xaxorig <- xax
  if (all(is.na(x))) {
    x <- 1:ncol(df)
  }
  if(xorder) {
    x <- rank(med, na.last=TRUE,ties.method="first")
  }
  ## new

  ### this was smarterized, make sure it doesn't break other code!!!
  # d <- ifelse(length(x)>1, diff(x[1:2]), 1)  #########
  # d <- ifelse(length(x)>1, diff(range(x, na.rm=TRUE))/length(x), 1)  #########

  ## new
  maxl <- min(length(x), 20)
  if(xaxlog) {
    d <- ifelse(length(x) > 1, diff(range(log(x), na.rm = TRUE))/maxl, 1)
  } else {
    d <- ifelse(length(x) > 1, diff(range(x, na.rm = TRUE))/maxl, 1)
  }
  # halfwidth <- 0.2*d  ## is the old way
  # mextra <- 0.3*d  ## is the old way

  halfwidth <- 0.2*d*(1 - 0.05*(length(x)^(-2)))
  mextra <- padwd*halfwidth#d*length(x)^(-2)
  ## new

  nn <- ncol(df)
  if(all(is.na(xax))) xax<-names(df)
  lwds <- (1+2*(1:length(ci)-1))*lwd
  if(!add) {
    # if(is.null(ylim)) ylim <- range(loq,hiq,na.rm=T)
    ci_lim <- range(loq, hiq, na.rm=TRUE)

    # xlims <- range(x-(.5*d*(1+(length(x)==1))),x+(.5*d*(1+(length(x)==1))))

    ## new
    if(xaxlog) {
      # xlims <- exp(range(log(x) - ((halfwidth+mextra)* (1 + (length(x) == 1))),
      #                    log(x) + ((halfwidth+mextra) * (1 + (length(x) == 1)))))
      ### was xlims
      cat_lim <- exp(c(min(log(x)) - ((halfwidth+mextra)* (1 + (length(x) == 1))),
                       max(log(x)) + ((halfwidth+mextra) * (1 + (length(x) == 1)))))
    } else {
      # xlims <- range(x - ((halfwidth+mextra) * (1 + (length(x) == 1))),
      #                x + ((halfwidth+mextra) * (1 + (length(x) == 1))))
      cat_lim <- c(min(x) - ((halfwidth+mextra) * (1 + (length(x) == 1))),
                   max(x) + ((halfwidth+mextra) * (1 + (length(x) == 1))))
    }
    if(!horizontal) {
      if(!is.null(xlim)) {
        cat_lim <- xlim
      }
      if(!is.null(ylim)) {
        ci_lim <- ylim
      }
      cat_lim <- range(cat_lim, xlim_add)
      ci_lim <- range(ci_lim, ylim_add)
    } else {
      if(!is.null(xlim)) {
        ci_lim <- xlim
      }
      if(!is.null(ylim)) {
        cat_lim <- ylim
      }
      cat_lim <- range(cat_lim, ylim_add)
      ci_lim <- range(ci_lim, xlim_add)
    }
    ## new

    # plot(NA, type='l', xlim=xlims, xlab=xlab, ylab=ylab, main=main, ylim=ylim, xaxt="n", ...=...)
    # axis(1,x,labels=xax, las=list(...)$las)

    ## new
    if(!horizontal) {
      if(all(is.na(xorig)) | !all(is.na(xaxorig))) {
        plot(NA, type = "l", xlim = cat_lim, xlab = xlab, ylab = ylab,
             main = main, ylim = ci_lim, xaxt = "n", ... = ...)
        smashaxis(side=1, at=x,
             labels = xax[seq_along(x)],
             las = list(...)$las,
             cex.axis = list(...)$cex.axis,
             smashall = all_ticks)
      } else {
        plot(NA, type = "l", xlim = cat_lim, xlab = xlab, ylab = ylab,
             main = main, ylim = ci_lim, ... = ...)
      }
    } else {
      if(all(is.na(xorig)) | !all(is.na(xaxorig))) {
        plot(NA, type = "l", ylim = rev(cat_lim), xlab = xlab, ylab = ylab,
             main = main, xlim = ci_lim, yaxt = "n", ... = ...)
        smashaxis(side=2, at=x,
             labels = xax[seq_along(x)],
             las = list(...)$las,
             cex.axis = list(...)$cex.axis,
             smashall = all_ticks)
      } else {
        plot(NA, type = "l", ylim = cat_lim, xlab = xlab, ylab = ylab,
             main = main, xlim = ci_lim, ... = ...) # ylim = rev(xlims)
      }
    }
    ## new
  }
  if(median) {
    # segments(x0=x-.2*d*medwd,x1=x+.2*d*medwd,y0=med,y1=med,col=col,lwd=medlwd, lend=1)

    ## new
    if(!horizontal) {
      if(xaxlog) {
        segments(x0 = exp(log(x) - halfwidth * medwd),
                 x1 = exp(log(x) + halfwidth * medwd),
                 y0 = med, y1 = med, col = col, lwd = medlwd,
                 lend = 1)
      } else {
        segments(x0 = x - halfwidth * medwd,
                 x1 = x + halfwidth * medwd,
                 y0 = med, y1 = med, col = col, lwd = medlwd,
                 lend = 1)
      }
    } else {
      if(xaxlog) {
        segments(y0 = exp(log(x) - halfwidth * medwd),
                 y1 = exp(log(x) + halfwidth * medwd),
                 x0 = med, x1 = med, col = col, lwd = medlwd,
                 lend = 1)
      } else {
        segments(y0 = x - halfwidth * medwd,
                 y1 = x + halfwidth * medwd,
                 x0 = med, x1 = med, col = col, lwd = medlwd,
                 lend = 1)
      }
    }
    ## new
  }
  # if(mean) points(x,colMeans(df, na.rm=T), pch=16, col=col)
  # for(i in 1:length(ci)) segments(x0=x,x1=x,y0=loq[i,],y1=hiq[i,],col=col,lwd=lwds[i],lend=1)

  ## new
  if(!horizontal) {
    if (mean) {
      points(x=x, y=colMeans(df, na.rm = T), pch = 16, col = col)
    }
    for (i in 1:length(ci)) {
      segments(x0 = x, x1 = x,
               y0 = loq[i, ], y1 = hiq[i, ],
               col = col, lwd = lwds[i], lend = 1)
    }
  } else {
    if (mean) {
      points(y=x, x=colMeans(df, na.rm = T), pch = 16, col = col)
    }
    for (i in 1:length(ci)) {
      segments(y0 = x, y1 = x,
               x0 = loq[i, ], x1 = hiq[i, ],
               col = col, lwd = lwds[i], lend = 1)
    }
  }
  ## new
}

# a helper function to smash all axis labels together and
# force all to print
smashaxis <- function(side, at=NULL, labels=TRUE,
                      smash=FALSE, smashall=FALSE, ...) {
  if(!is.null(at) & length(labels) <= 1) {
    labels <- at
  }
  if(!is.null(at) & (length(at) == length(labels)) & smash) {
    at1 <- at[seq(from=1, by=2, to=length(at))]
    at2 <- at[seq(from=2, by=2, to=length(at))]
    labels1 <- labels[seq(from=1, by=2, to=length(at))]
    labels2 <- labels[seq(from=2, by=2, to=length(at))]

    axis(side=side, at=at1, labels=labels1, ...=...)
    axis(side=side, at=at2, labels=labels2, ...=...)
  }
  if(!is.null(at) & (length(at) == length(labels)) & smashall) {
    for(iax in seq_along(at)) {
      axis(side=side, at=at[iax], labels=labels[iax], ...=...)
    }
  }
  if(!smash & !smashall) {
    axis(side=side, at=at, labels=labels, ...=...)
  }

  # axis(side=side, at=at, labels=labels, ...=...)
}
# plot(1:10, yaxt='n')
# axis(2, cex.axis=.8)
# smashaxis(side=2, at=1:10, labels=letters[1:10],
#           las=1, cex.axis=1.8)

# SS_df <- jags_df(SS_out)
# trend <- pull_post(SS_df, "trend")
# rate <- pull_post(SS_df, "rate")
# caterpillar(rate)
# caterpillar(rate, all_ticks=TRUE)
# caterpillar(rate, all_ticks=TRUE, las=2)
# caterpillar(rate, all_ticks=TRUE, las=2, horizontal=TRUE)
# caterpillar(rate, all_ticks=TRUE, las=2, horizontal=TRUE,
#             cex.axis=.2)
# caterpillar(rate, all_ticks=TRUE, las=2, horizontal=FALSE,
#             cex.axis=.2)
# caterpillar(rate, all_ticks=TRUE, las=2, horizontal=FALSE,
#             cex.axis=.2, xlab="xlab", ylab="ylab")
# caterpillar(rate, all_ticks=TRUE, las=2, horizontal=FALSE,
#             cex.axis=.2, xlab="xlab", ylab="ylab", cex.lab=.4)


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
#' @seealso \link{caterpillar}, \link{crossplot}, \link{comparedens}, \link{comparepriors}
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


