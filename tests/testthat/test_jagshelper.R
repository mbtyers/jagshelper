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
  expect_equal(dim(out_df), c(3000,8))
  expect_equal(sum(out_df), 303953, tolerance = 0.1)
  expect_error(jags_df(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("pull_post", {
  out_df <- jags_df(asdf_jags_out)
  a <- pull_post(out_df, "a")
  expect_true(inherits(a, "data.frame"))
  expect_equal(dim(a), c(3000,3))
  expect_equal(dim(pull_post(out_df)), c(3000,8))
  expect_equal(sum(a), 34.72302, tolerance = 0.1)
  expect_error(pull_post(asdf_jags_out), "Input must be a data.frame")
})

test_that("jags_plist", {
  out_plist <- jags_plist(asdf_jags_out)
  expect_equal(length(out_plist),8)
  dims <- lapply(out_plist, dim)
  expect_true(all(sapply(dims,"[",1)==1000))
  expect_true(all(sapply(dims,"[",2)==3))
})

test_that("trace_jags", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(trace_jags(asdf_jags_out))
  expect_silent(trace_jags(asdf_jags_out, p="a", lwd=c(2,2)))
  expect_error(trace_jags(asdf_jags_out, p="steve"), "No parameters with matching names")
  expect_error(trace_jags(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("chaindens_jags", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(chaindens_jags(asdf_jags_out))
  expect_silent(chaindens_jags(asdf_jags_out, p="a", lwd=c(2,2)))
  expect_error(chaindens_jags(asdf_jags_out, p="steve"), "No parameters with matching names")
  expect_error(chaindens_jags(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("tracedens_jags", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(tracedens_jags(asdf_jags_out))
  expect_silent(tracedens_jags(asdf_jags_out, p="a", lwd=c(2,2)))
  expect_error(tracedens_jags(asdf_jags_out, p="steve"), "No parameters with matching names")
  expect_error(tracedens_jags(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("nparam", {
  out_df <- jags_df(asdf_jags_out)
  expect_equal(nparam(asdf_jags_out),8)
  expect_equal(nparam(SS2d_jags), 501)
  expect_error(nparam(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("nbyname", {
  out_df <- jags_df(asdf_jags_out)
  expect_true(inherits(nbyname(asdf_jags_out),"list"))
  expect_equal(sum(unlist(nbyname(asdf_jags_out))), 8)
  expect_equal(nbyname(asdf_jags_out)$a, 3)
  expect_equal(nbyname(SS2d_jags)$rate, c(61,2))
  expect_error(nbyname(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("check_Rhat", {
  out_df <- jags_df(asdf_jags_out)
  expect_equal(length(check_Rhat(SS2d_jags)),9)
  expect_equal(sum(unlist(check_Rhat(asdf_jags_out))), 6)
  expect_equal(sum(unlist(check_Rhat(SS2d_jags))), 7.778689,tolerance=0.0001)
  expect_equal(unname(check_Rhat(SS2d_jags)[4]), 0.6967213, tolerance=0.0001)
  expect_equal(names(check_Rhat(SS2d_jags)[4]), "trend")
  expect_error(check_Rhat(out_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("check_neff", {
  out_df <- jags_df(asdf_jags_out)
  expect_equal(length(check_neff(SS2d_jags)),9)
  expect_equal(sum(unlist(check_neff(asdf_jags_out))), 6)
  expect_equal(sum(unlist(check_neff(SS2d_jags))), 2.188525,tolerance=0.0001)
  expect_equal(unname(check_neff(SS2d_jags)[4]), 0.2868852, tolerance=0.0001)
  expect_equal(names(check_neff(SS2d_jags)[4]), "trend")
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
  expect_silent(trace_df(out_df, nline=3, parmfrow=c(2,2)))
  expect_error(trace_df(asdf_jags_out), "Input must be a data.frame")
})

test_that("chaindens_df", {
  out_df <- jags_df(asdf_jags_out)
  expect_silent(chaindens_df(out_df, nline=3))
  expect_silent(chaindens_df(out_df, nline=3, parmfrow=c(2,2)))
  expect_error(chaindens_df(asdf_jags_out), "Input must be a data.frame")
})

test_that("envelope", {
  SS_df <- jags_df(SS_jags)
  trend <- pull_post(SS_df, "trend")
  expect_silent(envelope(trend, x=seq(2010,2015,by=1/12)))
  expect_silent(envelope(trend, x=seq(2010,2015,by=1/12) ,ci=.5))
  expect_silent(envelope(trend, x=seq(2010,2015,by=1/12) ,ci=c(.1,.5,.9)))
  expect_silent(envelope(SS_jags, p="trend"))
  expect_silent(envelope(SS2d_jags, p="trend", column=1))
  expect_silent(envelope(SS2d_jags, p="rate", column=1, main="rate"))
  expect_silent(envelope(SS2d_jags, p="rate", column=2, col=2, add=TRUE))

  expect_error(envelope(SS_jags, p=c("trend","rate")), "Need single parameter name in p= argument")
  expect_error(envelope(SS_jags), "Need single parameter name in p= argument")
  expect_error(envelope(SS_jags, p="steve"), "No parameters with matching names")
})

test_that("overlayenvelope", {
  expect_silent(overlayenvelope(df=list(SS2d_jags$sims.list$rate[,,1],SS2d_jags$sims.list$rate[,,2])))
  expect_silent(overlayenvelope(df=SS2d_jags$sims.list$rate))
  expect_silent(overlayenvelope(df=SS2d_jags, p="rate"))
  expect_silent(overlayenvelope(df=SS_jags, p=c("trend","rate")))
  expect_silent(overlayenvelope(df=SS2d_jags, p=c("trend","rate"), column=2))
  expect_silent(overlayenvelope(df=SS_jags, p=c("trend","rate"), legendnames=c("bob","larry")))
  expect_error(overlayenvelope(df=SS2d_jags))
})

test_that("caterpillar", {
  SS_df <- jags_df(SS_jags)
  trend <- pull_post(SS_df, "trend")
  expect_silent(caterpillar(trend, x=seq(2010,2015,by=1/12)))
  expect_silent(caterpillar(trend, x=seq(2010,2015,by=1/12) ,ci=.5))
  expect_silent(caterpillar(trend, x=seq(2010,2015,by=1/12) ,ci=c(.1,.5,.9)))
  expect_silent(caterpillar(SS_jags, p="trend"))
  expect_silent(caterpillar(SS2d_jags, p="trend", column=1))
  expect_silent(caterpillar(SS2d_jags, p="rate", column=1, main="rate"))
  expect_silent(caterpillar(SS2d_jags, p="rate", column=2, col=2, add=TRUE))

  expect_error(caterpillar(SS_jags, p=c("trend","rate")), "Need single parameter name in p= argument")
  expect_error(caterpillar(SS_jags), "Need single parameter name in p= argument")
  expect_error(caterpillar(SS_jags, p="steve"), "No parameters with matching names")
})

test_that("traceworstRhat", {
  expect_silent(traceworstRhat(SS_jags, parmfrow=c(3,2)))
  expect_silent(traceworstRhat(SS_jags, parmfrow=c(3,2), n.eff=TRUE))
  expect_silent(traceworstRhat(x=SS2d_jags, p=c("trend","cycle"), margin=2, parmfrow=c(2,2)))
  expect_silent(traceworstRhat(x=SS2d_jags, p=c("trend","cycle"), margin=2, parmfrow=c(2,2), n.eff=TRUE))

  SS_df <- jags_df(SS_jags)
  expect_error(traceworstRhat(SS_df), "Input must be an output object returned from jagsUI::jags().")
})

test_that("rcolors", {
  expect_equal(length(rcolors(10)),10)
})

test_that("plotRhats", {
  expect_silent(plotRhats(SS_jags))
  expect_silent(plotRhats(SS_jags, n.eff=TRUE))
  expect_silent(plotRhats(SS2d_jags))
  expect_silent(plotRhats(SS2d_jags, splitarr=TRUE, n.eff=TRUE))
  expect_silent(plotRhats(SS2d_jags, splitarr=TRUE, margin=2))
  expect_silent(plotRhats(SS2d_jags, p=c("trend", "cycle"), splitarr=TRUE, plotsequence=TRUE))
  SS_df <- jags_df(SS_jags)
  expect_error(plotRhats(SS_df), "Input must be an output object returned from jagsUI::jags().")
  expect_error(plotRhats(SS2d_jags, p="steve"), "No parameters with matching names")
})

test_that("comparedens", {
  expect_silent(comparedens(x1=asdf_jags_out, x2=asdf_jags_out, p=c("a","b","sig")))
  out_df <- jags_df(asdf_jags_out)
  expect_silent(comparedens(x1=out_df, x2=asdf_jags_out, p=c("a","b","sig")))
  expect_silent(comparedens(x2=out_df, x1=asdf_jags_out, p=c("a","b","sig")))
  expect_error(comparedens(x2=1:10, x1=asdf_jags_out, p=c("a","b","sig")),"Inputs must be data.frames or output objects returned from jagsUI::jags().")
})

test_that("comparecat", {
  expect_silent(comparecat(x=list(asdf_jags_out, asdf_jags_out, asdf_jags_out),p=c("a","b","sig")))
  expect_error(comparecat(x=list(1:10, asdf_jags_out, asdf_jags_out),p=c("a","b","sig")))
  expect_error(comparecat(x=asdf_jags_out))
})

test_that("pairstrace_jags", {
  expect_silent(pairstrace_jags(SS2d_jags, p="sig", parmfrow=c(2,3), lwd=2))
  expect_silent(pairstrace_jags(SS2d_jags, p="sig", parmfrow=c(2,3), points=TRUE))
  expect_silent(pairstrace_jags(SS2d_jags, p="sig", parmfrow=c(2,3), contour=TRUE))
  expect_silent(pairstrace_jags(asdf_jags_out, parmfrow=c(3,3)))
  expect_silent(pairstrace_jags(asdf_jags_out, parmfrow=c(3,3), points=TRUE))
  expect_silent(pairstrace_jags(asdf_jags_out, parmfrow=c(3,3), contour=TRUE))
  out_df <- jags_df(asdf_jags_out)
  expect_error(pairstrace_jags(out_df),"Input must be an output object returned from jagsUI::jags().")
})

