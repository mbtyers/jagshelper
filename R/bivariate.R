#' Bivariate Plot of Posterior Densities
#' @description Bivariate plot of the posterior densities of corresponding vectors
#' of parameter nodes.  Three plotting methods are provided, that may be overlayed
#' if desired.
#'
#' * If `drawcross == TRUE`, [caterpillar]-like plots will be produced, with quantile
#' intervals in the x- and y- directions.
#'
#' * If `drawx == TRUE`, [caterpillar]-like plots will be produced, but rotated
#' along the standardized principal component axes.  This may be useful to draw if correlation
#' is present.
#'
#' * If `drawblob == TRUE`, smoothed polygons will be produced, each containing
#' approximately `ci=` x100% of the associated MCMC samples.
#'
#' All methods can overlay multiple bars or polygons, depending on the length of `ci=`.
#' @param dfx Output object returned from `jagsUI::jags()`; or alternately,
#' two-dimensional `data.frame` or matrix in which parameter node element is
#' given by column and MCMC iteration is given by row.  A vector may also be used,
#' that expresses MCMC iterations of a single parameter node.  If used with `dfy=`,
#' this will be plotted in the x-direction.
#' @param dfy Optionally, a
#' two-dimensional `data.frame` or matrix in which parameter node element is
#' given by column and MCMC iteration is given by row.  A vector may also be used,
#' that expresses MCMC iterations of a single parameter node.  If used, this will
#' be plotted in the y-direction.
#' @param p Vector of parameter names, if input to `dfx` is a `jagsUI` output object.
#' If used, this must be of length 2.
#' @param col Color for plotting, or recyclable vector of colors.  Defaults to `4`.
#' If `col == "random"`, [rcolors] will be used to generate a random vector of colors.
#' @param drawcross Whether to draw quantile bars in the x- and y-directions.
#' Defaults to `TRUE`.
#' @param drawx Whether to draw quantile bars along the standardized principal component axes.
#' Defaults to `FALSE`.
#' @param drawblob Whether to draw smoothed quantile polygons.
#' Defaults to `FALSE`.
#' @param blobres Optional tuning parameter for drawing quantile polygons, and
#' corresponds to the number of polygon vertices.  If the default `NULL` is accepted,
#' the function will supply a value based on the number of MCMC samples.
#' @param blobsmooth Optional tuning parameter for drawing quantile polygons, and
#' corresponds to half the number of polygon vertices used for local smoothing.
#' If the default `NULL` is accepted,
#' the function will supply a value based on the number of MCMC samples and the
#' number of vertices.
#' @param outline Whether to draw quantile polygons as lines rather than filled regions.  Defaults to `FALSE`.
#' @param ci Vector of intervals to overlay.  Defaults to 50 percent and 95 percent.
#' @param lwd Base line width for plotting.  Defaults to 1.
#' @param mean Whether to include points for means.  Defaults to `FALSE`.
#' @param link Whether to link medians in sequence.  Defaults to `FALSE`.
#' @param linklwd Line width to use for linking.  Defaults to `1`.
#' @param labels Whether to add labels, or a vector of labels to add.  Defaults to `FALSE`.
#' @param labelpos Optionally, an argument to `pos` in \link[graphics]{text} for labels.  Defaults to `NULL`.
#' @param labelcex Optional character expansion for labels.  Defaults to `0.7`.
#' @param whichx Element to subset for x, if only one element of a  vector of parameter nodes is desired for plotting.
#' @param rowx Row to subset for x, in the case of a 2-d matrix of parameter nodes in-model.
#' @param columnx Column to subset for x, in the case of a 2-d matrix of parameter nodes in-model.
#' @param whichy Element to subset for x, if only one element of a  vector of parameter nodes is desired for plotting.
#' @param rowy Row to subset for y, in the case of a 2-d matrix of parameter nodes in-model.
#' @param columny Column to subset for y, in the case of a 2-d matrix of parameter nodes in-model.
#' @param xlab X-axis label.  If the default `NULL` is accepted, this will be drawn automatically.
#' @param ylab Y-axis label.  If the default `NULL` is accepted, this will be drawn automatically.
#' @param main Plot title.
#' @param xlim X-axis limits.  If the default (`NULL`) is accepted, the limits will be determined automatically.
#' @param ylim Y-axis limits.  If the default (`NULL`) is accepted, the limits will be determined automatically.
#' @param transformx Should the x-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transformx="exp"`is used, consider
#' adding additional plotting argument `log="x"` or `log="xy"`.
#' @param transformy Should the y-axis be (back)transformed?  Options are `"exp"`,
#' indicating exponential, or `"expit"`, indicating inverse-logit. Defaults to
#' `"none"`, indicating no transformation.  Note: if `transformy="exp"`is used, consider
#' adding additional plotting argument `log="y"` or `log="xy"`.
#' @param add Whether to add to existing plot
#' @param ... additional plotting arguments
#' @return `NULL`
#' @seealso \link{caterpillar}, \link{pairstrace_jags}
#' @author Matt Tyers
#' @importFrom stats princomp qnorm sd
#' @examples
#' ## basic functionality with cross geometry
#' crossplot(SS_out, p=c("trend","rate"))
#'
#' ## default labels
#' crossplot(SS_out, p=c("trend","cycle"), labels=TRUE)
#'
#' ## showing:
#' ## - link lines
#' ## - blob geometry (smoothed confidence polygons)
#' ## - random colors with col="random"
#' crossplot(SS_out, p=c("trend","cycle"),
#'           labels=SS_data$x, labelpos=1, link=TRUE, drawblob=TRUE,
#'           col="random")
#'
#' ## adding x geometry and showing usage with a single vector element (41)
#' crossplot(SS_out, p=c("trend","cycle"),
#'           whichx=41, whichy=41,
#'           drawblob=TRUE, drawx=TRUE)
#'
#' ## single vectors (or data.frames or 2d matrices) can also be used
#' xx <- SS_out$sims.list$trend[,41]
#' yy <- SS_out$sims.list$cycle[,41]
#'
#' par(mfrow = c(2, 2))
#' plot(xx, yy, col=adjustcolor(1, alpha.f=.1), pch=16, main="cross geometry")
#' crossplot(xx, yy, add=TRUE, col=1)
#' plot(xx, yy, col=adjustcolor(1, alpha.f=.1), pch=16, main="x geometry")
#' crossplot(xx, yy, add=TRUE, col=1,
#'           drawcross=FALSE, drawx=TRUE)
#' plot(xx, yy, col=adjustcolor(1, alpha.f=.1), pch=16, main="blob geometry")
#' crossplot(xx, yy, add=TRUE, col=1,
#'           drawcross=FALSE, drawblob=TRUE)
#' plot(xx, yy, col=adjustcolor(1, alpha.f=.1), pch=16, main="blob outlines")
#' crossplot(xx, yy, add=TRUE, col=1,
#'           drawcross=FALSE, drawblob=TRUE, outline=TRUE)
#' @export
crossplot <- function(dfx,
                      dfy=NULL,
                      p=NULL,

                      col=4,

                      drawcross=TRUE, drawx=FALSE, drawblob=FALSE,
                      blobres=NULL, blobsmooth=NULL,
                      outline=FALSE,

                      ci=c(0.5,0.95),
                      lwd=1,
                      mean=FALSE,

                      link=FALSE, linklwd=1,
                      labels=FALSE, labelpos=NULL, labelcex=0.7,

                      whichx=NULL, rowx=NULL, columnx=NULL,
                      whichy=NULL, rowy=NULL, columny=NULL,

                      xlab=NULL, ylab=NULL, main=NULL, xlim=NULL, ylim=NULL,

                      transformx=c("none", "exp", "expit"),
                      transformy=c("none", "exp", "expit"),

                      add=FALSE,   ...) {


  if(!inherits(dfx,"jagsUI")) {
    if(!inherits(dfx,c("matrix","data.frame","numeric","integer")) &
       !inherits(dfy,c("matrix","data.frame","numeric","integer"))) {
      stop("Input must be a data.frames or output from jagsUI::jags() plus two parameter names")
    }
  }
  if(inherits(dfx,"jagsUI") & length(p)!=2) stop("Need two parameter names in p= argument") ###

  if(inherits(dfx,"jagsUI") & length(p)==2) {
    simslist <- dfx$sims.list
    if(any(!p %in% names(simslist))) stop("No parameters with matching names") ###

    theparmx <- simslist[names(simslist)==p[1]][[1]]
    if(is.null(rowx) & is.null(columnx) & is.null(whichx)) {
      dfx <- theparmx
    } else {
      theparmx <- simslist[names(simslist)==p[1]][[1]]
      if(!is.null(whichx)) dfx <- theparmx[,whichx]
      if(!is.null(rowx)) dfx <- theparmx[,rowx,]
      if(!is.null(columnx)) dfx <- theparmx[,,columnx]
    }
    if(is.null(ylab)) {
      if(is.null(rowx) & is.null(columnx) & is.null(whichx)) {
        xlab <- p[1]
      } else {
        if(!is.null(whichx)) {
          if(length(whichx)==1) {
            xlab <- paste0(p[1],"[",whichx,"]")
          } else {
            xlab <- p[1]
          }
        }
        if(!is.null(rowx)) xlab <- paste0(p[1],"[",rowx,",]")
        if(!is.null(columnx)) xlab <- paste0(p[1],"[,",columnx,"]")
      }
    }

    theparmy <- simslist[names(simslist)==p[2]][[1]]
    if(is.null(rowy) & is.null(columny) & is.null(whichy)) {
      dfy <- theparmy
    } else {
      theparmy <- simslist[names(simslist)==p[2]][[1]]
      if(!is.null(whichy)) dfy <- theparmy[,whichy]
      if(!is.null(rowy)) dfy <- theparmy[,rowy,]
      if(!is.null(columny)) dfy <- theparmy[,,columny]
    }
    if(is.null(ylab)) {
      if(is.null(rowy) & is.null(columny) & is.null(whichy)) {
        ylab <- p[2]
      } else {
        if(!is.null(whichy)) {
          if(length(whichy)==1) {
            ylab <- paste0(p[2],"[",whichy,"]")
          } else {
            ylab <- p[2]
          }
        }
        if(!is.null(rowy)) ylab <- paste0(p[2],"[",rowy,",]")
        if(!is.null(columny)) ylab <- paste0(p[2],"[,",columny,"]")
      }
    }
  }

  if(!isTRUE(all.equal(dim(dfx), dim(dfy)))) stop("Dimension or length mismatch between X and Y")

  dots <- list(...)   ####### this does not do anything at present

  if(is.null(main)) main <- ""
  if(is.null(xlab)) xlab <- ""
  if(is.null(ylab)) ylab <- ""

  dfx <- as.matrix(dfx)
  dfy <- as.matrix(dfy)

  transformx <- match.arg(transformx)
  if(transformx == "exp") dfx <- exp(dfx)
  if(transformx == "expit") dfx <- expit(dfx)
  transformy <- match.arg(transformy)
  if(transformy == "exp") dfy <- exp(dfy)
  if(transformy == "expit") dfy <- expit(dfy)

  ci <- rev(sort(ci))

  loqx <- apply(dfx, 2, quantile, p=(1-ci)/2, na.rm=T)
  hiqx <- apply(dfx, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  loqy <- apply(dfy, 2, quantile, p=(1-ci)/2, na.rm=T)
  hiqy <- apply(dfy, 2, quantile, p=1-(1-ci)/2, na.rm=T)
  if(length(ci)==1) {
    loqx <- t(as.matrix(loqx))
    hiqx <- t(as.matrix(hiqx))
    loqy <- t(as.matrix(loqy))
    hiqy <- t(as.matrix(hiqy))
  }
  medx <- apply(dfx, 2, median, na.rm=T)
  medy <- apply(dfy, 2, median, na.rm=T)

  nn <- ncol(dfx)
  # if(all(is.na(xax))) xax<-names(df)
  lwds <- (1+2*(1:length(ci)-1))*lwd

  # determine which columns can be used
  # whichcolumns <- which(colSums(!is.na(dfx))>0 & colSums(!is.na(dfy))>0)
  whichcolumns <- which(apply(dfx, 2, sd, na.rm=TRUE)>0 & apply(dfy, 2, sd, na.rm=TRUE)>0)

  # computing coords from x and blob to help set xlim & ylim
  x_coordlist <- NULL
  if(drawx) {
    x_coordlist <- list()
    i_x <- 1
    for(i_data in whichcolumns) {
      for(i_ci in seq_along(ci)) {
        x_coordlist[[i_x]] <- addx(xx=dfx[,i_data], yy=dfy[,i_data], ci=ci[i_ci], lwd=lwds[i_ci], col=cols[i_data],
                                    return_coords = TRUE)#
        i_x <- i_x + 1
      }
    }
  }

  blob_coordlist <- NULL
  if(drawblob) {
    blob_coordlist <- list()
    i_blob <- 1
    for(i_data in whichcolumns) {
      # for(i_data in which(!is.na(colSums(dfx, na.rm=T)) & !is.na(colSums(dfy, na.rm=T)))) {
      for(i_ci in seq_along(ci)) {
        blob_coordlist[[i_blob]] <- confidenceblob(x=dfx[,i_data], y=dfy[,i_data], ci=ci[i_ci], lwd=lwds[i_ci], col=cols[i_data],
                                              return_coords = TRUE,
                                              blobres=blobres, blobsmooth=blobsmooth)#
        i_blob <- i_blob + 1
      }
    }
  }

  if(!add) {
    if(is.null(xlim)) xlim <- range(loqx, hiqx,
                                    sapply(x_coordlist, function(x) x$x),
                                    sapply(blob_coordlist, function(x) x$x),
                                    na.rm=TRUE)
    if(is.null(ylim)) ylim <- range(loqy, hiqy,
                                    sapply(x_coordlist, function(x) x$y),
                                    sapply(blob_coordlist, function(x) x$y),
                                    na.rm=TRUE)
    plot(NA, type='l', xlim=xlim, xlab=xlab, ylab=ylab, main=main, ylim=ylim, ...=...)# , xaxt="n"
    # axis(1,x,labels=xax, las=list(...)$las)
  }

  if(mean) points(colMeans(dfx, na.rm=T), colMeans(dfy, na.rm=T), pch=16, col=1)

  if(isTRUE(col=="random")) {
    cols <- rcolors(ncol(dfx))
  } else {
    cols <- rep(col, length.out=ncol(dfx))  # this makes colors recycle
  }
  if(drawcross) {
    for(i in 1:length(ci)) {
      segments(x0=medx, x1=medx, y0=loqy[i,], y1=hiqy[i,], col=cols, lwd=lwds[i], lend=1)
      segments(y0=medy, y1=medy, x0=loqx[i,], x1=hiqx[i,], col=cols, lwd=lwds[i], lend=1)
    }
  }

  if(drawx) {
    i_x <- 1
    for(i_data in whichcolumns) {
      for(i_ci in seq_along(ci)) {
        # addx3(xx=dfx[,i_data], yy=dfy[,i_data], ci=ci[i_ci], lwd=lwds[i_ci], col=cols[i_data])#
        # if(!return_coords) segments(x0=corners[c(1,3),1], y0=corners[c(1,3),2],
        #                             x1=corners[c(2,4),1], y1=corners[c(2,4),2],
        #                             lwd=lwd, col=col, lty=lty)
        # if(return_coords) return(data.frame(x=corners[,1], y=corners[,2]))
        segments(x0=x_coordlist[[i_x]]$x[c(1,3)], y0=x_coordlist[[i_x]]$y[c(1,3)],
                 x1=x_coordlist[[i_x]]$x[c(2,4)], y1=x_coordlist[[i_x]]$y[c(2,4)],
                 lwd=lwds[i_ci], col=cols[i_data])# , lty=lty
        i_x <- i_x + 1
      }
    }
  }

  if(drawblob) {
    i_blob <- 1
    for(i_data in whichcolumns) {
    # for(i_data in which(!is.na(colSums(dfx, na.rm=T)) & !is.na(colSums(dfy, na.rm=T)))) {
      for(i_ci in seq_along(ci)) {
        # confidenceblob2(x=dfx[,i_data], y=dfy[,i_data], ci=ci[i_ci], lwd=lwds[i_ci], col=cols[i_data])#
        if(outline) {
          fillcol <- NA
          linecol <- adjustcolor(cols[i_data],
                                 alpha.f=seq(to=1, length.out=length(ci), by=1/length(ci))[i_ci])
        } else {
          fillcol <- adjustcolor(cols[i_data], alpha.f=.3)
          linecol <- NA
        }
        # polygon(blob_coordlist[[i_blob]]$x, blob_coordlist[[i_blob]]$y,
        #         col=adjustcolor(cols[i_data], alpha.f=.3), border=NA)
        polygon(blob_coordlist[[i_blob]]$x, blob_coordlist[[i_blob]]$y,
                col=fillcol, border=linecol, lwd=lwd)
        i_blob <- i_blob + 1
      }
    }
  }

  if(link) lines(x=medx, y=medy, col=1, lwd=linklwd)
  if(isTRUE(labels)) labels <- seq_along(medx)
  if(length(labels) > 1) {
    text(x=medx, y=medy, labels=labels, pos=labelpos, cex=labelcex)
  }
}

# addx <- function(xx, yy, ci, lwd=1, col=1, lty=1) {
#   covwt <- cov.wt(na.omit(cbind(xx, yy)))
#   eigenthings <- eigen(covwt$cov)
#   lengths <- sqrt(2*outer(eigenthings$values, qf(ci, 2, length(xx)-1)))
#   corners <- matrix(c(
#     covwt$center + outer(lengths[1] , eigenthings$vectors[,1]),
#     covwt$center - outer(lengths[1] , eigenthings$vectors[,1]),
#     covwt$center + outer(lengths[2] , eigenthings$vectors[,2]),
#     covwt$center - outer(lengths[2] , eigenthings$vectors[,2])), nrow=2)
#   segments(x0=rep(covwt$center[1], ncol(corners)), y0=rep(covwt$center[2],  ncol(corners)),
#            x1=corners[1,], y1=corners[2,],
#            lwd=lwd, col=col, lty=lty)
# }
#
# addx2 <- function(xx, yy, ci, lwd=1, col=1, lty=1) {
#   xmn <- mean(xx, na.rm=TRUE)
#   ymn <- mean(yy, na.rm=TRUE)
#   xsd <- sd(xx, na.rm=TRUE)
#   ysd <- sd(yy, na.rm=TRUE)
#   x1 <- (xx-xmn)/xsd
#   y1 <- (yy-ymn)/ysd
#
#   covwt <- cov.wt(na.omit(cbind(x1, y1)))
#   eigenthings <- eigen(covwt$cov)
#   lengths <- sqrt(2*outer(eigenthings$values, qf(ci, 2, length(xx)-1)))
#   corners1 <- matrix(c(
#     covwt$center + outer(lengths[1] , eigenthings$vectors[,1]),
#     covwt$center - outer(lengths[1] , eigenthings$vectors[,1]),
#     covwt$center + outer(lengths[2] , eigenthings$vectors[,2]),
#     covwt$center - outer(lengths[2] , eigenthings$vectors[,2])), nrow=2)
#   corners <- NA*corners1
#   corners[1,] <- corners1[1,]*xsd + xmn
#   corners[2,] <- corners1[2,]*ysd + ymn
#
#   segments(x0=rep(xmn, ncol(corners)), y0=rep(ymn,  ncol(corners)),
#            x1=corners[1,], y1=corners[2,],
#            lwd=lwd, col=col, lty=lty)
# }

addx <- function(xx, yy, ci, lwd=1, col=1, lty=1, return_coords = FALSE) {
  mnx <- mean(xx, na.rm=T)
  mny <- mean(yy, na.rm=T)
  sdx <- sd(xx, na.rm=T)
  sdy <- sd(yy, na.rm=T)

  x1 <- (xx-mnx)/sdx
  y1 <- (yy-mny)/sdy
  pc1 <- princomp(cbind(x1,y1))

  ci2 <- c((1-ci)/2, 1-((1-ci)/2))

  corners1 <- data.frame(x=c(quantile(pc1$scores[,1], ci2), rep(0,2)),
                         y=c(rep(0,2), quantile(pc1$scores[,2], ci2)))

  corners2 <- as.matrix(corners1) %*% MASS::ginv(pc1$loadings)

  corners <- NA*corners2
  corners[,1] <- corners2[,1]*sdx + mnx
  corners[,2] <- corners2[,2]*sdy + mny

  if(!return_coords) segments(x0=corners[c(1,3),1], y0=corners[c(1,3),2],
           x1=corners[c(2,4),1], y1=corners[c(2,4),2],
           lwd=lwd, col=col, lty=lty)
  if(return_coords) return(data.frame(x=corners[,1], y=corners[,2]))
}

# confidenceblob <- function(x, y, res=round(length(x)^.5), ci=0.95, lwd=1, col=1, lty=1) {
#
#   xmed <- median(x, na.rm=TRUE)
#   ymed <- median(y, na.rm=TRUE)
#
#   r <- sqrt(((x-xmed)^2) + ((y-ymed)^2))
#   theta <- atan2(y-ymed, x-xmed)
#
#   # anglebreaks <- seq(from=-pi, to=pi, length.out=res+1)
#   anglebreaks <- unname(quantile(theta, seq(from=0, to=1, length.out=res+1)))
#
#   thetacut <- cut(theta, breaks=anglebreaks, include.lowest=TRUE, ordered_result = TRUE)
#   # thetacut <- cut(theta, breaks=unname(quantile(theta, seq(from=0, to=1, length.out=res+1))), include.lowest=TRUE)
#   blobr <- tapply(r, thetacut, quantile, p=ci, na.rm=TRUE)
#   # blobtheta <- (anglebreaks[-1] + anglebreaks[-(res+1)])/2
#   blobtheta <- tapply(theta, thetacut, median)
#
#   blobx <- xmed + blobr*cos(blobtheta)
#   bloby <- ymed + blobr*sin(blobtheta)
#   # blobr_smooth <- predict(loess(blobr ~ seq_along(blobr)))
#   # blobx <- xmed + blobr_smooth*cos(blobtheta)
#   # bloby <- ymed + blobr_smooth*sin(blobtheta)
#   blobx <- c(blobx, blobx[1])
#   bloby <- c(bloby, bloby[1])
#
#   # segments(x0=rep(xmed, length(blobr)),
#   #          y0=rep(ymed, length(blobr)),
#   #          x1=xmed+max(blobr, na.rm=TRUE)*cos(blobtheta),
#   #          y1=ymed+max(blobr, na.rm=TRUE)*sin(blobtheta), col=adjustcolor(1,alpha.f=.2))
#
#   lines(blobx, bloby, lwd=lwd, col=col, lty=lty)
# }

confidenceblob <- function(x, y,
                           blobres=NULL, blobsmooth=NULL,
                           ci=0.95, lwd=1, col=1, lty=1, return_coords = FALSE) {

  if(is.null(blobres)) blobres <- round(length(x)/50)
  # if(is.null(blobsmooth)) blobsmooth <- round(100*blobres/length(x))
  if(is.null(blobsmooth)) blobsmooth <- round(blobres^0.25)

  xmed <- median(x, na.rm=TRUE)
  ymed <- median(y, na.rm=TRUE)
  xsd <- sd(x, na.rm=TRUE)
  ysd <- sd(y, na.rm=TRUE)

  x1 <- (x-xmed)/xsd
  y1 <- (y-ymed)/ysd

  r <- sqrt((x1^2) + (y1^2))
  theta <- atan2(y1, x1)

  anglebreaks <- unname(quantile(theta, seq(from=0, to=1, length.out=blobres+1)))

  thetacut <- cut(theta, breaks=anglebreaks, include.lowest=TRUE, ordered_result = TRUE)
  blobr <- tapply(r, thetacut, quantile, p=ci, na.rm=TRUE)
  blobtheta <- tapply(theta, thetacut, median)

  smooth <- T     ## instead of this, add tuneable smoothing parameters
  if(!smooth) {
    blobx1 <- blobr*cos(blobtheta)
    bloby1 <- blobr*sin(blobtheta)
  } else {
    # blobr_smooth <- predict(loess(blobr ~ seq_along(blobr)))
    blobr_smooth <- NA*blobr
    for(i in seq_along(blobr)) {
      mod_ind <- (i + -blobsmooth:blobsmooth) %% length(blobr)
      mod_ind[mod_ind == 0] <- length(blobr)
      blobr_smooth[i] <- mean(blobr[mod_ind])
    }
    blobx1 <- blobr_smooth*cos(blobtheta)
    bloby1 <- blobr_smooth*sin(blobtheta)
  }

  blobx1 <- c(blobx1, blobx1[1])
  bloby1 <- c(bloby1, bloby1[1])

  blobx <- blobx1*xsd + xmed
  bloby <- bloby1*ysd + ymed

  # segments(x0=rep(xmed, length(blobr)),
  #          y0=rep(ymed, length(blobr)),
  #          x1=xmed+max(blobr, na.rm=TRUE)*cos(blobtheta),
  #          y1=ymed+max(blobr, na.rm=TRUE)*sin(blobtheta), col=adjustcolor(1,alpha.f=.2))

  # lines(blobx, bloby, lwd=lwd, col=col, lty=lty)
  if(!return_coords) {
    ## include option to outline=T, style it after envelope if possible
    polygon(blobx, bloby, col=adjustcolor(col, alpha.f=.3), border=NA)
  }
  if(return_coords) return(data.frame(x=blobx, y=bloby))
}
