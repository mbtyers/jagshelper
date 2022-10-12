
<!-- README.md is generated from README.Rmd. Please edit that file -->

# jagshelper

<!-- badges: start -->
<!-- badges: end -->

The goal of `jagshelper` is to streamline Bayesian analysis in `JAGS`
using the `jagsUI` package.  
Functions are provided for extracting output in a simpler form,
assessing model convergence, and plotting model output. Also included is
a function giving a template model in `JAGS` syntax with the associated
`jagsUI` code.

## Installation

You can install the development version of `jagshelper` like so:

``` r
devtools::install_github("mbtyers/jagshelper")
```

## Model template

The `skeleton()` function prints a `JAGS` model template to the screen,
along with an associated simulated dataset.

``` r
library(jagshelper)
skeleton("EXAMPLE")
```

## Assessing convergence - a simple model

In the example below, there are relatively few parameters saved, and it
is feasible to examine the trace plots associated with each parameter.

``` r
nparam(asdf_jags_out)  # how many parameters in total
nbyname(asdf_jags_out)  # how many parameters (or dimensions) per parameter name
tracedens_jags(asdf_jags_out, parmfrow=c(2,2))  #trace plots
check_Rhat(asdf_jags_out)  # proportion of Rhats below a threshold of 1.1
plotRhats(asdf_jags_out)  # plotting Rhat values
```

## Assessing convergence - a more complex model

In the example below, there are relatively many parameters saved, and it
is perhaps more illustrative to examine the trace plots associated with
the least- converged parameters, as measured by `Rhat` value.

``` r
nparam(SS_out)  # how many parameters in total
nbyname(SS_out)  # how many parameters (or dimensions) per parameter name
traceworstRhat(SS_out, parmfrow=c(2,2))  #trace plots
check_Rhat(SS_out)  # proportion of Rhats below a threshold of 1.1
plotRhats(SS_out)  # plotting Rhat values
```

## Extracting output as data.frame

The `jags_df()` function extracts the full posterior from an output
object returned by `jagsUI::jags()` as a `data.frame`, which may be
preferable for some users.

``` r
out_df <- jags_df(asdf_jags_out)
str(out_df)
```

## Visualizing posteriors of vectors of parameter nodes

The `caterpillar()` and `envelope()` functions plot output for vectors
of parameter nodes, and `comparecat()`, `comparedens()` and
`overlayenvelope()` functions allow comparison between multiple models
or parameter vectors.

``` r
caterpillar(asdf_jags_out,"a")
envelope(SS_out, "rate", x=SS_data$x)
```
