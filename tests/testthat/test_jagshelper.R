### tests
library(testthat)
library(jagshelper)
# library(jagsUI)

test_that("skeleton", {
  expect_output(skeleton("a"))
})

test_that("jags_df", {
  out_df <- jags_df(asdf_jags_out)
  expect_true(inherits(out_df, "data.frame"))
  expect_equal(dim(out_df), c(1500,8))
  expect_equal(sum(out_df), 259036, tolerance = 0.1)
  expect_error(jags_df(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("pull_post", {
  out_df <- jags_df(asdf_jags_out)
  a <- pull_post(out_df, "a")
  expect_true(inherits(a, "data.frame"))
  expect_equal(dim(a), c(1500,3))
  expect_equal(dim(pull_post(out_df)), c(1500,8))
  expect_equal(sum(a), 108.1291, tolerance = 0.1)
  expect_error(pull_post(asdf_jags_out), "Input must be a data.frame")
})

test_that("jags_plist", {
  out_plist <- jags_plist(asdf_jags_out)
  expect_equal(length(out_plist),8)
  dims <- lapply(out_plist, dim)
  expect_true(all(sapply(dims,"[",1)==500))
  expect_true(all(sapply(dims,"[",2)==3))
})

test_that("trace_jags", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(trace_jags(asdf_jags_out))
  expect_silent(trace_jags(asdf_jags_out, p="a", parmfrow=c(2,2)))
  expect_error(trace_jags(asdf_jags_out, p="steve"), "No parameters with matching names")
  expect_error(trace_jags(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("chaindens_jags", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(chaindens_jags(asdf_jags_out))
  expect_silent(chaindens_jags(asdf_jags_out, p="a", parmfrow=c(2,2)))
  expect_error(chaindens_jags(asdf_jags_out, p="steve"), "No parameters with matching names")
  expect_error(chaindens_jags(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("tracedens_jags", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(tracedens_jags(asdf_jags_out))
  expect_silent(tracedens_jags(asdf_jags_out, p="a", parmfrow=c(2,2)))
  expect_error(tracedens_jags(asdf_jags_out, p="steve"), "No parameters with matching names")
  expect_error(tracedens_jags(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("nparam", {
  out_df <- jags_df(asdf_jags_out)
  expect_equal(nparam(asdf_jags_out),8)
  expect_equal(nparam(SS_out), 334)
  expect_error(nparam(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("nbyname", {
  out_df <- jags_df(asdf_jags_out)
  expect_true(inherits(nbyname(asdf_jags_out),"list"))
  expect_equal(sum(unlist(nbyname(asdf_jags_out))), 8)
  expect_equal(nbyname(asdf_jags_out)$a, 3)
  expect_equal(nbyname(SS_out)$cycle_s, c(41,2))
  expect_error(nbyname(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("check_Rhat", {
  out_df <- jags_df(asdf_jags_out)
  expect_equal(length(check_Rhat(SS_out)),12)
  expect_equal(sum(unlist(check_Rhat(asdf_jags_out))), 6)
  expect_equal(sum(unlist(check_Rhat(SS_out))), 10.5863,tolerance=0.0001)
  expect_equal(unname(check_Rhat(SS_out)[9]), 0.5121951, tolerance=0.0001)
  expect_equal(names(check_Rhat(SS_out)[9]), "cycle_s")
  expect_error(check_Rhat(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("check_neff", {
  out_df <- jags_df(asdf_jags_out)
  expect_equal(length(check_neff(SS_out)),12)
  expect_equal(sum(unlist(check_neff(asdf_jags_out))), 5)
  expect_equal(sum(unlist(check_neff(SS_out))), 1.012195,tolerance=0.0001)
  expect_equal(unname(check_neff(SS_out)[9]), 0.03658537, tolerance=0.0001)
  expect_equal(names(check_neff(SS_out)[9]), "cycle_s")
  expect_error(check_neff(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("logit", {
  expect_equal(logit(0.5),0)
})

test_that("expit", {
  expect_equal(expit(0),0.5)
})

test_that("trace_line", {
  expect_silent(trace_line(1:30, nline=3))
  out_df <- jags_df(asdf_jags_out)
  b1 <- pull_post(out_df, "b1")
  a <- pull_post(out_df, "a")
  expect_silent(trace_line(b1, nline=3))
})

test_that("chaindens_line", {
  expect_silent(chaindens_line(1:30, nline=3))
  out_df <- jags_df(asdf_jags_out)
  b1 <- pull_post(out_df, "b1")
  a <- pull_post(out_df, "a")
  expect_silent(chaindens_line(b1, nline=3))
})

test_that("trace_df", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(trace_df(out_df, nline=3))
  expect_silent(trace_df(as.matrix(out_df), nline=3))
  expect_silent(trace_df(out_df, nline=3, parmfrow=c(2,2)))
  expect_error(trace_df(asdf_jags_out), "Input must be a data.frame")
})

test_that("chaindens_df", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(chaindens_df(out_df, nline=3))
  expect_silent(chaindens_df(as.matrix(out_df), nline=3))
  expect_silent(chaindens_df(out_df, nline=3, parmfrow=c(2,2)))
  expect_error(chaindens_df(asdf_jags_out), "Input must be a data.frame")
})

test_that("envelope", {
  SS_df <- jags_df(SS_out)
  trend <- pull_post(SS_df, "trend")
  expect_silent(envelope(trend, x=SS_data$x))
  expect_silent(envelope(trend, x=SS_data$x ,ci=.5))
  expect_silent(envelope(trend, x=SS_data$x ,ci=c(.1,.5,.9)))
  expect_silent(envelope(SS_out, p="trend"))
  expect_silent(envelope(SS_out, p="cycle_s", column=1))
  expect_silent(envelope(SS_out, p="cycle_s", column=1, main="cycle_s"))
  expect_silent(envelope(SS_out, p="cycle_s", column=2, col=2, add=TRUE))

  expect_error(envelope(SS_out, p=c("trend","rate")), "Need single parameter name in p= argument")
  expect_error(envelope(SS_out), "Need single parameter name in p= argument")
  expect_error(envelope(SS_out, p="steve"), "No parameters with matching names")

  expect_silent(envelope(trend, x=SS_data$x, transform="exp"))
  expect_silent(envelope(trend, x=SS_data$x, transform="exp", log="y"))
  expect_silent(envelope(SS_out, p="trend", transform="exp", log="y"))
  expect_silent(envelope(SS_out, p="trend", transform="expit", log="y", ylab="transform it"))
  expect_error(envelope(SS_out, p="trend", transform="somethingsalwayswrong", log="y"))
})

test_that("overlayenvelope", {
  expect_silent(overlayenvelope(df=list(SS_out$sims.list$cycle_s[,,1],SS_out$sims.list$cycle_s[,,2])))
  expect_silent(overlayenvelope(df=SS_out$sims.list$cycle_s))
  expect_silent(overlayenvelope(df=SS_out, p="cycle_s"))
  expect_silent(overlayenvelope(df=SS_out, p=c("trend","rate")))
  expect_silent(overlayenvelope(df=SS_out, p=c("trend","cycle_s"), column=2))
  expect_silent(overlayenvelope(df=SS_out, p=c("trend","rate"), legendnames=c("bob","larry")))
  expect_error(overlayenvelope(df=SS_out))

  expect_silent(overlayenvelope(df=SS_out, p="cycle_s", transform="exp"))
  expect_silent(overlayenvelope(df=SS_out, p="cycle_s", transform="exp", log="y"))
  expect_error(overlayenvelope(df=SS_out, p="cycle_s", transform="bob"))
})

test_that("caterpillar", {
  SS_df <- jags_df(SS_out)
  trend <- pull_post(SS_df, "trend")
  expect_silent(caterpillar(trend, x=SS_data$x))
  expect_silent(caterpillar(trend, x=SS_data$x ,ci=.5))
  expect_silent(caterpillar(trend, x=SS_data$x ,ci=c(.1,.5,.9)))
  expect_silent(caterpillar(SS_out, p="trend"))
  expect_silent(caterpillar(SS_out, p="trend", xax=c(letters, LETTERS)[1:41]))
  expect_silent(caterpillar(SS_out, p="trend", xax=c(letters, LETTERS)[1:41], las=2))
  expect_silent(caterpillar(SS_out, p="trend", ylim=c(-10, 10)))
  expect_silent(caterpillar(SS_out, p="cycle_s", column=1))
  expect_silent(caterpillar(SS_out, p="cycle_s", column=1, main="cycle_s"))
  expect_silent(caterpillar(SS_out, p="cycle_s", column=2, col=2, add=TRUE))

  expect_error(caterpillar(SS_out, p=c("trend","rate")), "Need single parameter name in p= argument")
  expect_error(caterpillar(SS_out), "Need single parameter name in p= argument")
  expect_error(caterpillar(SS_out, p="steve"), "No parameters with matching names")

  expect_silent(caterpillar(trend, x=SS_data$x, transform="exp", log="y"))
  expect_silent(caterpillar(SS_out, p="trend", transform="exp"))
  expect_silent(caterpillar(SS_out, p="trend", transform="expit"))
  expect_error(caterpillar(SS_out, p="trend", transform="somethingelse"))
})

test_that("traceworstRhat", {
  expect_silent(traceworstRhat(SS_out, parmfrow=c(3,2)))
  expect_silent(traceworstRhat(SS_out, parmfrow=c(3,2), n.eff=TRUE))
  expect_silent(traceworstRhat(x=SS_out, p="cycle_s", margin=2, parmfrow=c(2,2)))
  expect_silent(traceworstRhat(x=SS_out, p="cycle_s", margin=2, parmfrow=c(2,2), n.eff=TRUE))

  SS_df <- jags_df(SS_out)
  expect_error(traceworstRhat(SS_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("rcolors", {
  expect_equal(length(rcolors(10)),10)
})

test_that("plotRhats", {
  expect_silent(plotRhats(SS_out))
  expect_silent(plotRhats(SS_out, n.eff=TRUE))
  expect_silent(plotRhats(SS_out))
  expect_silent(plotRhats(SS_out, splitarr=TRUE, n.eff=TRUE))
  expect_silent(plotRhats(SS_out, splitarr=TRUE, margin=2))
  expect_silent(plotRhats(SS_out, p=c("trend", "cycle"), splitarr=TRUE, plotsequence=TRUE))
  SS_df <- jags_df(SS_out)
  expect_error(plotRhats(SS_df), "Input must be an output object returned from jagsUI::jags().")
  expect_error(plotRhats(SS_out, p="steve"), "No parameters with matching names")
})

test_that("comparedens", {
  expect_silent(comparedens(x1=asdf_jags_out, x2=asdf_jags_out, p=c("a","b","sig")))
  out_df <- jags_df(asdf_jags_out)
  expect_silent(comparedens(x1=out_df, x2=asdf_jags_out, p=c("a","b","sig")))
  expect_silent(comparedens(x1=out_df, x2=asdf_jags_out, p=c("a","b","sig"), ylim=c(-10,10)))
  expect_silent(comparedens(x1=out_df, x2=asdf_jags_out, p=c("a","b","sig"), col=2:3))
  expect_silent(comparedens(x2=out_df, x1=asdf_jags_out, p=c("a","b","sig")))
  expect_error(comparedens(x2=1:10, x1=asdf_jags_out, p=c("a","b","sig")),"Inputs must be data.frames or output objects returned from jagsUI::jags().")
})

test_that("comparecat", {
  expect_silent(comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),p=c("a","b","sig")))
  expect_silent(comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),p=c("a","b","sig"), col=1:3))
  expect_silent(comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),p=c("a","b","sig"), col=1:3, transform="exp"))
  expect_silent(comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),p=c("a","b","sig"), col=1:3, transform="expit"))
  expect_error(comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),p=c("a","b","sig"), col=1:3, transform="badbadbad"))
  expect_error(comparecat(x=list(1:10, asdf_jags_out, asdf_jags_out),p=c("a","b","sig")))
  expect_error(comparecat(x=asdf_jags_out))
})

test_that("pairstrace_jags", {
  expect_silent(pairstrace_jags(SS_out, p="sig", parmfrow=c(2,3), lwd=2))
  expect_silent(pairstrace_jags(SS_out, p="sig", parmfrow=c(2,3), points=TRUE))
  expect_silent(pairstrace_jags(SS_out, p="sig", parmfrow=c(2,3), contour=TRUE))
  expect_silent(pairstrace_jags(asdf_jags_out, parmfrow=c(3,3)))
  expect_silent(pairstrace_jags(asdf_jags_out, parmfrow=c(3,3), points=TRUE))
  expect_silent(pairstrace_jags(asdf_jags_out, parmfrow=c(3,3), contour=TRUE))
  out_df <- jags_df(asdf_jags_out)
  expect_error(pairstrace_jags(out_df),"Input must be an output object returned from jagsUI::jags().")
})

test_that("cor_jags", {
  expect_equal(sum(cor_jags(asdf_jags_out)), 8.291098, tolerance=0.0001)
  expect_equal(dim(cor_jags(asdf_jags_out)), c(8,8))
  expect_equal(dim(cor_jags(asdf_jags_out, p=c("a","b"))), c(5,5))
})

test_that("plotcor_jags", {
  expect_silent(plotcor_jags(asdf_jags_out))
  expect_silent(plotcor_jags(jags_df(asdf_jags_out)))
  expect_silent(asdf_jags_out$sims.list$a)
  expect_silent(plotcor_jags(asdf_jags_out, p=c("a","b")))
  expect_silent(plotcor_jags(asdf_jags_out, legend=F, mincor=0.1, maxn=1))
})

test_that("plotdens", {
  expect_silent(plotdens(asdf_jags_out, p="b1"))
  expect_silent(plotdens(asdf_jags_out, p="b1", add=T, shade=F, lwd=F))
  expect_silent(plotdens(asdf_jags_out, p="a", minCI=.95, col=2:4))
  expect_silent(plotdens(asdf_jags_out, p=c("a[1]","a[2]","a[3]"), legend=F))
  expect_silent(plotdens(jags_df(asdf_jags_out, p="a"),legendnames=c("albert","betty","chuck")))
  expect_silent(plotdens(list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p="b1"))
  expect_error(plotdens(list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p="a"),"No parameter names are an exact match to p= argument.")
  expect_silent(plotdens(list(asdf_jags_out,asdf_jags_out,asdf_jags_out), p=c("a[1]","a[2]","a[3]")))
})

test_that("qq_postpred", {
  expect_silent(qq_postpred(ypp=SS_out, p="ypp", y=SS_data$y))
  expect_silent(qq_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y))
  expect_silent(qq_postpred(ypp=jags_df(x=SS_out, p="ypp"), y=SS_data$y))
  expect_silent(qq_postpred(ypp=SS_out, p="ypp", y=SS_data$y, add=T, pch="+", col=4))
  expect_error(qq_postpred(ypp=SS_out, p="ypp", y=SS_data$y[1]))
  expect_error(qq_postpred(ypp=SS_out, p="ypp", y=SS_data$y[1:2]), "Posterior matrix ypp has more columns than length of data matrix y")
  expect_warning(qq_postpred(ypp=SS_out$sims.list$ypp[,1:2], y=SS_data$y), "Posterior matrix ypp has fewer columns than length of data matrix y")
  expect_silent(qq_postpred(ypp=SS_out$sims.list$ypp[,1:2], y=SS_data$y[1:2]))
  expect_error(qq_postpred(ypp=SS_out, y=SS_data$y), "Parameter name must be supplied to p= argument if jagsUI object is used in argument ypp")
})

test_that("ts_postpred", {
  expect_silent(ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y))
  expect_silent(ts_postpred(ypp=SS_out, x=SS_data$x, p="ypp", y=SS_data$y))
  expect_silent(ts_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y))
  expect_silent(ts_postpred(ypp=jags_df(x=SS_out, p="ypp"), y=SS_data$y))
  expect_silent(ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y, add=T, col=3))
  expect_error(ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y[1]))
  expect_error(ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y[1:2]), "Posterior matrix ypp must have the same number of columns as length of data matrix y")
  expect_error(ts_postpred(ypp=SS_out, y=SS_data$y), "Parameter name must be supplied to p= argument if jagsUI object is used in argument ypp")

  expect_silent(ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y, transform="exp"))
  expect_silent(ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y, transform="expit"))
  expect_error(ts_postpred(ypp=SS_out, p="ypp", y=SS_data$y, transform="larryboy"))
})

test_that("comparepriors", {
  expect_silent(comparepriors(x=asdf_prior_jags_out))
  expect_silent(comparepriors(x=asdf_prior_jags_out, parmfrow=c(3,2)))
  expect_silent(comparepriors(x=asdf_prior_jags_out, parmfrow=c(3,2), col=3:2, minCI=0.7, legendpos="bottomleft"))
  expect_warning(comparepriors(x=asdf_jags_out), 'No parameter names ending in "_prior"')
})

trend1 <- trend2 <- SS_out$sims.list$trend
rate1 <- rate2 <- SS_out$sims.list$rate
trend1[, 2:3] <- rate1[, 2:3] <- NA
trend2[, 2:3] <- rate2[, 2:3] <- 42
test_that("crossplot", {
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend[,1], dfy=SS_out$sims.list$rate[,1]))
  expect_silent(crossplot(dfx=as.data.frame(SS_out$sims.list$trend),
                          dfy=as.data.frame(SS_out$sims.list$rate)))
  expect_silent(crossplot(dfx=log(20+SS_out$sims.list$trend),
                          dfy=log(20+SS_out$sims.list$rate),
                          transformx="exp", transformy="exp"))
  expect_silent(crossplot(dfx=suppressWarnings(log(SS_out$sims.list$trend)),
                          dfy=suppressWarnings(log(SS_out$sims.list$rate))))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate,
                          drawcross=FALSE))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate,
                          drawcross=FALSE, drawx=TRUE, drawblob=FALSE))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate,
                          drawcross=FALSE, drawx=FALSE, drawblob=TRUE))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate,
                          drawcross=FALSE, drawx=FALSE, drawblob=TRUE,
                          outline=TRUE, lwd=2))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate,
                          drawcross=TRUE, drawx=TRUE, drawblob=TRUE,
                          lwd=1, link=TRUE, col=2, linklwd=3, labels=TRUE))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate,
                          drawcross=TRUE, drawx=TRUE, drawblob=TRUE,
                          lwd=1, link=TRUE, col="random", linklwd=3, labels=TRUE))
  expect_silent(crossplot(dfx=SS_out$sims.list$trend, dfy=SS_out$sims.list$rate,
                          drawcross=TRUE, drawx=TRUE, drawblob=TRUE,
                          labels=SS_data$x, labelpos=1, labelcex=1.2))
  expect_silent(crossplot(dfx=trend1, dfy=rate1, drawblob=TRUE, drawx=TRUE))
  expect_silent(crossplot(dfx=trend2, dfy=rate2, drawblob=TRUE, drawx=TRUE))
  expect_silent(crossplot(dfx=SS_out, p=c("trend","rate")))
  expect_silent(crossplot(dfx=SS_out, p=c("trend","rate"), whichx=7, whichy=7))
  expect_silent(crossplot(dfx=SS_out, p=c("trend","rate"), whichx=7:10, whichy=7:10))
  expect_silent(crossplot(dfx=SS_out, p=c("trend","sig_eps"), whichx=7))
  expect_silent(crossplot(dfx=SS_out, p=rev(c("trend","sig_eps")), whichy=7))
  expect_error(crossplot(dfx=SS_out, p=c("trend","sig_eps")),
               "Dimension or length mismatch between X and Y")
  expect_error(crossplot(dfx=SS_out, p=c("trend","cycle_s")),
               "Dimension or length mismatch between X and Y")
  expect_silent(crossplot(dfx=SS_out, p=rev(c("trend","cycle_s")), columnx = 1))
  expect_silent(crossplot(dfx=SS_out, p=c("trend","cycle_s"), columny = 1))
})

test_that("plot_postpred", {
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              lines=TRUE))
  expect_silent(plot_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y, x=SS_data$x))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              plot_residuals=FALSE))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              plot_data=FALSE))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              plot_sd=FALSE))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              whichplots=1))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              whichplots=2))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              whichplots=3))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              whichplots=4))
  expect_silent(plot_postpred(ypp=SS_out, p="ypp", y=SS_data$y, x=SS_data$x,
                              whichplots=1:4,
                              pointcol=2+(SS_data$y>2),
                              pch=2+(SS_data$x>2018)))
})

x_withNA <- SS_data$x
y_withNA <- SS_data$y
ypp_withNA <- SS_out$sims.list$ypp
x_withNA[sample(1:41, 2)] <- NA
y_withNA[sample(1:41, 2)] <- NA
ypp_withNA[,sample(1:41, 2)] <- NA
test_that("NA cases in _postpred", {
  expect_silent(qq_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y))
  expect_silent(qq_postpred(ypp=ypp_withNA, y=SS_data$y))
  expect_silent(qq_postpred(ypp=SS_out$sims.list$ypp, y=y_withNA))
  expect_silent(qq_postpred(ypp=ypp_withNA, y=y_withNA))

  expect_silent(ts_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y, x=SS_data$x))
  expect_silent(ts_postpred(ypp=ypp_withNA, y=SS_data$y, x=SS_data$x))
  expect_silent(ts_postpred(ypp=SS_out$sims.list$ypp, y=y_withNA, x=SS_data$x))
  expect_silent(ts_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y, x=x_withNA))
  expect_silent(ts_postpred(ypp=ypp_withNA, y=y_withNA, x=x_withNA))

  expect_silent(envelope(df=SS_out$sims.list$ypp, x=SS_data$x))
  expect_silent(envelope(df=ypp_withNA, x=SS_data$x))
  expect_silent(envelope(df=ypp_withNA, x=x_withNA))

  expect_silent(plot_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y, x=SS_data$x))
  expect_silent(plot_postpred(ypp=ypp_withNA, y=SS_data$y, x=SS_data$x))
  expect_silent(plot_postpred(ypp=SS_out$sims.list$ypp, y=y_withNA, x=SS_data$x))
  expect_silent(plot_postpred(ypp=SS_out$sims.list$ypp, y=SS_data$y, x=x_withNA))
  expect_silent(plot_postpred(ypp=ypp_withNA, y=y_withNA, x=x_withNA))
})
