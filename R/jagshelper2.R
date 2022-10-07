# comparedens <- function(x1,x2,...) {
#   xdf1 <- jags_df(x1)
#   xdf2 <- jags_df(x2)
#
#   parmx1 <- cbind(pull_post(xdf1,"sig"), phi=pull_post(xdf1,"phi"))
#   parmx2 <- cbind(pull_post(xdf2,"sig"), phi=pull_post(xdf2,"phi"))
#
#   allparms <- sort(unique(c(names(parmx1),names(parmx2))))
#
#   plot(NA,xlim=c(0,length(allparms)+1), ylim=range(parmx1,parmx2,na.rm=T), ylab="",xlab="",xaxt="n",...=...)
#   axis(1,at=1:length(allparms),labels=allparms,las=2)
#   abline(v=1:length(allparms),lty=3)
#
#   for(i in 1:length(allparms)) {
#     if(allparms[i] %in% names(parmx1)) {
#       xxx <- density(parmx1[,which(names(parmx1)==allparms[i])])
#       polygon(i-xxx$y/max(xxx$y)/2, xxx$x, col=adjustcolor(2,alpha.f=.5), border=2)
#     }
#     if(allparms[i] %in% names(parmx2)) {
#       xxx <- density(parmx2[,which(names(parmx2)==allparms[i])])
#       polygon(i+xxx$y/max(xxx$y)/2, xxx$x, col=adjustcolor(4,alpha.f=.5), border=4)
#     }
#   }
# }
# # comparedens(x1,x2)
#
# comparecat <- function(x,...) {
#   xdf <- lapply(x, jags_df)
#   parmx <- list()
#   for(i in 1:length(x)) parmx[[i]] <- cbind(pull_post(xdf[[i]],"sig"), phi=pull_post(xdf[[i]],"phi"))
#
#   allparms <- sort(unique(unlist(lapply(parmx,names))))
#
#   plot(NA,xlim=c(0,length(allparms)+1), ylim=range(parmx,na.rm=T), ylab="",xlab="",xaxt="n",...=...)
#   axis(1,at=1:length(allparms),labels=allparms,las=2)
#   # abline(v=1:length(allparms),lty=3)
#
#   for(i in 1:length(allparms)) {
#     for(ii in 1:length(x)) {
#       if(allparms[i] %in% names(parmx[[ii]])) {
#         xplot <- i+(seq(-.3,.3,length.out=length(x)))[ii]
#         vec <- parmx[[ii]][,which(names(parmx[[ii]])==allparms[i])]
#         segments(x0=rep(xplot,2), y0=quantile(vec,p=c(.025,.25)), y1=quantile(vec,p=c(.975,.75)),
#                  lwd=c(1,3), lend=1, col=ii+1)
#         points(xplot,median(vec),pch=16,col=ii+1)
#       }
#     }
#   }
# }
#
#
# tracedens_jags1 <- function (x, p = NULL, parmfrow = NULL, lwd = 1, shade=T, diagnostics=F, ...)
# {
#   x_plist <- jags_plist1(x, p = p)
#   if (!is.null(parmfrow)) {
#     parmfrow1 <- par("mfrow")
#     par(mfrow = parmfrow)
#   }
#   nline <- ncol(x_plist[[1]])
#   cols <- adjustcolor(rainbow(nline), red.f = 0.9, blue.f = 0.9,
#                       green.f = 0.9, alpha.f = 0.6)
#   for (i in 1:length(x_plist)) {
#     allbw <- density(as.vector(x_plist[[i]]))$bw
#     denses <- apply(x_plist[[i]], 2, density, bw = allbw)
#     dx <- sapply(denses, function(x) x$x)
#     dy <- sapply(denses, function(x) x$y)
#     dy1 <- (dy/max(dy) * 0.3 + 1) * nrow(x_plist[[i]])
#     plot(NA, xlim = c(1, 1.3 * nrow(x_plist[[i]])), ylim = range(x_plist[[i]]),
#          main = names(x_plist)[i], xlab = "iter", ylab = "",
#          ... = ...)
#     for (j in 1:nline) lines(x_plist[[i]][, j], col = cols[j],
#                              lwd = lwd)
#     for (j in 1:nline) lines(dy1[, j], dx[, j], col = cols[j],
#                              lwd = lwd)
#     if(shade) {
#       for (j in 1:nline) polygon(dy1[, j], dx[, j], border = NA, #cols[j],
#                                  col=adjustcolor(cols[j], alpha.f=.2))
#     }
#     if(diagnostics) {    #### this doesnt work if part of a vector
#       Rhat <- unlist(x$Rhat[names(x$Rhat)==names(x_plist)[i]])
#       neff <- unlist(x$n.eff[names(x$n.eff)==names(x_plist)[i]])
#       legend("topright",legend=c(paste("Rhat:",round(Rhat,3)), paste("neff:",round(neff,0))), bty='n')
#     }
#   }
#   if (!is.null(parmfrow))
#     par(mfrow = parmfrow1)
# }
#
#
# traceworstRhat <- function(x, p) {
#   thevec <- unname(unlist(x$Rhat[names(x$Rhat)==p]))
#   whichone <- which.max(thevec)
#   thename <- paste0(p,"[",whichone,"]")
#   tracedens_jags(x,p=thename)
# }
#
#
# traceworstRhat_mat <- function(x, p) {
#   # thevec <- unname(unlist(x$Rhat[names(x$Rhat)==p]))
#   themat <- x$Rhat[names(x$Rhat)==p][[1]]
#   # whichone <- which.max(thevec)
#   whichones <- apply(themat,2,which.max)
#   thenames <- paste0(p,"[",whichones,",",1:length(whichones),"]")
#   tracedens_jags1(x,p=thenames)
# }
#
# jags_plist1 <- function (x, p = NULL) {   # need to make this addition to jags_plist!!
#   x_dflist <- lapply(x$samples, as.data.frame)
#   x2 <- lapply(1:length(x_dflist[[1]]), function(x) sapply(x_dflist,
#                                                            "[[", x))
#   names(x2) <- names(x_dflist[[1]])
#   if (!is.null(p)) {
#     these <- rep(F, length(x2))
#     for(i in 1:length(p)) {
#       these[substr(names(x2), 1, nchar(p[i])) == p[i]] <- T
#     }
#     x2 <- x2[these]
#   }
#   return(x2)
# }
#
# envelope0 <- function (df, x = NA, median = T, ci = c(0.5, 0.95), col = 4,
#                        add = F, dark = 0.3, outline = F, xlab = "", ylab = "", main = "",
#                        ...)
# {
#   ci <- sort(ci)
#   loq <- apply(df, 2, quantile, p = (1 - ci)/2, na.rm = T)
#   hiq <- apply(df, 2, quantile, p = 1 - (1 - ci)/2, na.rm = T)
#   if (length(ci) == 1) {
#     loq <- t(as.matrix(loq))
#     hiq <- t(as.matrix(hiq))
#   }
#   med <- apply(df, 2, median, na.rm = T)
#   if (all(is.na(x)))
#     x <- 1:ncol(df)
#   if (!add) {
#     plot(NA, xlim = range(x), #ylim = range(loq, hiq, na.rm = T),
#          xlab = xlab, ylab = ylab, main = main, ... = ...)
#     if (median)
#       lines(x, med, col = col)
#   }
#   else if (median)
#     lines(x, med, type = "l", col = col, ... = ...)
#   if (outline) {
#     darks <- rev(1 - ((1 - dark)^(1:length(ci))))
#     for (i in 1:length(ci)) {
#       lines(x, loq[i, ], col = adjustcolor(col, alpha.f = darks[i]),
#             ... = ...)
#       lines(x, hiq[i, ], col = adjustcolor(col, alpha.f = darks[i]),
#             ... = ...)
#     }
#   }
#   else {
#     for (i in 1:length(ci)) polygon(c(x, rev(x)), c(loq[i,
#     ], rev(hiq[i, ])), col = adjustcolor(col, alpha.f = dark),
#     border = NA)
#   }
# }
# envelope2 <- function(df1, df2, col=NA, ...) {
#   if(is.na(col)) col <- c(2,4)
#   lim1 <- apply(df1,2,quantile,p=c(.025,.975),na.rm=T)
#   lim2 <- apply(df2,2,quantile,p=c(.025,.975),na.rm=T)
#   ylim <- range(lim1,lim2,na.rm=T)
#   envelope0(df1,col=col[1],ylim=ylim,...=...)
#   envelope0(df2,col=col[2],add=T,...=...)
# }
#
# pairstrace_jags <- function (x, p = NULL, parmfrow = NULL, lwd = 1, alpha=0.2,points=F,...)
# {
#   x_plist <- jags_plist1(x, p = p)
#   if (!is.null(parmfrow)) {
#     parmfrow1 <- par("mfrow")
#     par(mfrow = parmfrow)
#   }
#   nline <- ncol(x_plist[[1]])
#   cols <- adjustcolor(rainbow(nline), red.f = 0.9, blue.f = 0.9,
#                       green.f = 0.9, alpha.f = alpha)
#   nn <- length(x_plist)
#   for (i in 1:(nn-1)) {
#     for(j in (i+1):nn) {
#       plot(NA, xlim=range(x_plist[[i]], na.rm=T), ylim=range(x_plist[[j]]),
#            main="", xlab=names(x_plist)[i], ylab=names(x_plist)[j])
#       for (k in 1:nline) {
#         if(!points) {
#           lines(x=x_plist[[i]][, k], y=x_plist[[j]][, k], col = cols[k], lwd = lwd)
#         } else {
#           points(x=x_plist[[i]][, k], y=x_plist[[j]][, k], col = cols[k])
#         }
#       }
#     }
#   }
#   if (!is.null(parmfrow))
#     par(mfrow = parmfrow1)
# }
#
# pairstrace_jags_comp <- function (x, p = NULL, parmfrow = NULL, lwd = 1, alpha=0.2,points=F,...)
# {
#   # x_plist <- jags_plist1(x, p = p)
#   x_plist <- lapply(x, FUN=jags_plist1, p = p)
#   if (!is.null(parmfrow)) {
#     parmfrow1 <- par("mfrow")
#     par(mfrow = parmfrow)
#   }
#   nline <- ncol(x_plist[[1]][[1]])
#   cols <- adjustcolor(rainbow(nline), red.f = 0.9, blue.f = 0.9,
#                       green.f = 0.9, alpha.f = alpha)
#   nn <- length(x_plist[[1]])
#   for (i in 1:(nn-1)) {
#     for(j in (i+1):nn) {
#       for(kk in 1:length(x)) {
#         plot(NA, xlim=range(x_plist[[kk]][[i]], na.rm=T), ylim=range(x_plist[[kk]][[j]]),
#              main="", xlab=names(x_plist[[1]])[i], ylab=names(x_plist[[1]])[j])
#         for (k in 1:nline) {
#           if(!points) {
#             lines(x=x_plist[[kk]][[i]][, k], y=x_plist[[kk]][[j]][, k], col = cols[k], lwd = lwd)
#           } else {
#             points(x=x_plist[[kk]][[i]][, k], y=x_plist[[kk]][[j]][, k], col = cols[k])
#           }
#         }
#       }
#     }
#   }
#   if (!is.null(parmfrow))
#     par(mfrow = parmfrow1)
# }
#
# library(MASS)
# pairscontour_jags_comp <- function (x, p = NULL, parmfrow = NULL, lwd = 1, alpha=0.2,points=F,...)
# {
#   # x_plist <- jags_plist1(x, p = p)
#   x_plist <- lapply(x, FUN=jags_plist1, p = p)
#   if (!is.null(parmfrow)) {
#     parmfrow1 <- par("mfrow")
#     par(mfrow = parmfrow)
#   }
#   nline <- ncol(x_plist[[1]][[1]])
#   cols <- adjustcolor(rainbow(nline), red.f = 0.9, blue.f = 0.9,
#                       green.f = 0.9, alpha.f = alpha)
#   nn <- length(x_plist[[1]])
#   for (i in 1:(nn-1)) {
#     for(j in (i+1):nn) {
#       for(kk in 1:length(x)) {
#         # plot(NA, xlim=range(x_plist[[kk]][[i]], na.rm=T), ylim=range(x_plist[[kk]][[j]]),
#         #      main="", xlab=names(x_plist[[1]])[i], ylab=names(x_plist[[1]])[j])
#         # for (k in 1:nline) {
#         #   if(!points) {
#         #     lines(x=x_plist[[kk]][[i]][, k], y=x_plist[[kk]][[j]][, k], col = cols[k], lwd = lwd)
#         #   } else {
#         #     points(x=x_plist[[kk]][[i]][, k], y=x_plist[[kk]][[j]][, k], col = cols[k])
#         #   }
#         # }
#         xx <- as.vector(x_plist[[kk]][[i]])
#         yy <- as.vector(x_plist[[kk]][[j]])
#         if(sum(diff(xx))!=0 & sum(diff(yy))!=0) {
#           contour(MASS::kde2d(x=xx, y=yy), main="", xlab=names(x_plist[[1]])[i], ylab=names(x_plist[[1]])[j])
#         } else {
#           plot(NA, main="", xlab=names(x_plist[[1]])[i], ylab=names(x_plist[[1]])[j], xlim=0:1, ylim=0:1)
#         }
#       }
#     }
#   }
#   if (!is.null(parmfrow))
#     par(mfrow = parmfrow1)
# }
#
# tracequantileRhat_mat <- function(x, p, q=0.5) {
#   # thevec <- unname(unlist(x$Rhat[names(x$Rhat)==p]))
#   themat <- x$Rhat[names(x$Rhat)==p][[1]]
#   # whichone <- which.max(thevec)
#   whichones <- apply(themat,2,function(y) which(rank(y)==round(q*length(y),0)))
#   thenames <- paste0(p,"[",whichones,",",1:length(whichones),"]")
#   tracedens_jags1(x,p=thenames)
# }
