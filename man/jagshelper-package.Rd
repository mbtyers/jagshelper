\name{jagshelper-package}
\alias{jagshelper-package}
\alias{jagshelper}
\docType{package}
\title{
Functions for Extracting and Visualizing Output from 'jagsUI'
}
\description{
Functions are provided to help run Bayesian analyses in JAGS using
    the 'jagsUI' package.  Included are functions for extracting output in
    simpler format, functions for streamlining assessment of convergence, and
    functions for producing summary plots of output.  Also included is a
    function that provides a simple template for running JAGS from R.
}
\details{
\tabular{ll}{
Package: \tab jagshelper\cr
Type: \tab Package\cr
Version: \tab 0.1.10\cr
Date: \tab 2022-11-02\cr
License: \tab GPL-2\cr
}
The jagshelper package is intended to extend and streamline Bayesian analysis using the `jagsUI` package.

The \link{skeleton} function prints a template JAGS model with associated R code to the console, which can easily be copied & pasted to an R script and modified as needed.

Functions are also provided for visually assessing model convergence.  In particular, \link{tracedens_jags} gives a relatively simple syntax for trace plots of a collection or subset of parameter nodes, and overlays by-chain kernel densities for visual assessment of marginal posterior shapes as well as overlap between MCMC chains.  Another function that could be particularly useful to users is \link{plotRhats}, which gives a visual representation of the values of the Gelman-Rubin convergence diagnostic \code{Rhat} (or alternately effective sample size \code{n.eff}) for all saved parameters.  This may be particularly useful in the case where a model has many saved parameters.  Additionally, function \link{traceworstRhat} is a wrapper for \link{tracedens_jags}, but only produces trace plots for the parameter nodes with the worst (largest) values of \code{Rhat} or \code{n.eff}.

Functions are also provided for visualizing posterior densities; in particular, the case of a vector of parameter nodes (one-dimensional in the JAGS model, giving a two-dimensional matrix of MCMC iterations).  Notably, the \link{envelope} function is intended for a sequence of nodes (as in a time series), and the \link{caterpillar} function is intended for cases in which order may not matter (as in a collection of random effects).

Wrapper functions are also given for overlay of multiple such plots, as \link{overlayenvelope} and \link{comparecat}, and \link{comparedens} giving plots as vertically-oriented left- and right-facing kernel densities.
}
\author{
Matt Tyers

Maintainer: Matt Tyers <matttyersstat@gmail.com>
}
