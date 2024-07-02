#' Number of parameters
#' @description Total number of individual parameter nodes saved in `jagsUI` output.
#' @param x Output object from `jagsUI::jags()`
#' @return A single numeric value giving the number of parameter nodes.
#' @seealso \link{nbyname}
#' @author Matt Tyers
#' @examples
#' head(jags_df(asdf_jags_out))
#'
#' nparam(asdf_jags_out)
#' @export
nparam <- function(x) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(ncol(jags_df(x)))
}


#' Number of parameter nodes by parameter name
#' @description Returns a list of the numbers of parameter nodes saved in `jagsUI` output, by parameter name.
#' As a default, what is returned for each list element is a vector of the array dimensions within the JAGS model
#'  (that is, excluding the dimension associated with the number of MCMC samples for each parameter node),
#' or alternately, just the total number of parameter nodes.
#' @param x Output object from `jagsUI::jags()`
#' @param justtotal Whether to just report the total number of parameters, as opposed to dimensions.
#' @return A `list` with an element associated with each parameter.  Each element
#' can be interpreted as the vector length or array dimension associated with the
#' given parameter.
#' @seealso \link{nparam}
#' @author Matt Tyers
#' @examples
#' head(jags_df(asdf_jags_out))
#'
#' nbyname(asdf_jags_out)
#'
#' nparam(SS_out)
#' nbyname(SS_out)
#' @export
nbyname <- function(x, justtotal=FALSE) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  # sapply(x$sims.list, function(x) ncol(as.matrix(x)))
  out1 <- sapply(x$sims.list, function(x) dim(x)[-1])
  out1[sapply(out1, is.null)] <- 1
  if(!justtotal) return(out1) #comment
  if(justtotal) return(sapply(out1, prod))
}


#' Quick summary of Rhat values by parameter name
#' @description Returns the mean number of `Rhat` values for each parameter (by each parameter)
#' that are less than a specified threshold criterion.
#'
#' `Rhat` (Gelman-Rubin Convergence Diagnostic, or Potential Scale Reduction Factor)
#' is calculated within 'JAGS', and is
#' commonly used as a measure of convergence for a given parameter node.  Values close
#' to 1 are seen as evidence of adequate convergence.
#' @param x Output object from `jagsUI::jags()`
#' @param thresh Threshold value (defaults to 1.1)
#' @return Numeric (named) giving the proportion of Rhat values below the given threshold.
#' @seealso \link{check_neff}, \link{traceworstRhat}, \link{plotRhats}, \link{qq_postpred}, \link{ts_postpred}
#' @author Matt Tyers
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from Iterative Simulation
#' Using Multiple Sequences. *Statistical Science, 7*(4), 457–472. http://www.jstor.org/stable/2246093
#' @examples
#' check_Rhat(SS_out)
#' @export
check_Rhat <- function(x, thresh=1.1) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(sapply(x$Rhat, function(x) mean(x<thresh, na.rm=TRUE)))
}


#' Quick summary of n.eff values by parameter name
#' @description Returns the mean number of `n.eff` values (by each parameter) that are greater than a specified threshold criterion.
#'
#' `n.eff` is calculated within 'JAGS', and may be interpreted as a crude measure of
#' effective sample size for a given parameter node.
#' @param x Output object from `jagsUI::jags()`
#' @param thresh Threshold value (defaults to 500)
#' @return Numeric (named) giving the proportion of `n.eff` values above the given threshold.
#' @seealso \link{check_Rhat}, \link{traceworstRhat}, \link{plotRhats}, \link{qq_postpred}, \link{ts_postpred}
#' @author Matt Tyers
#' @examples
#' check_neff(SS_out)
#' @export
check_neff <- function(x, thresh=500) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  return(sapply(x$n.eff, function(x) mean(x>thresh, na.rm=TRUE)))
}


#' Correlation matrix from a JAGS object
#' @description Computes a correlation matrix of all MCMC samples from an object
#' returned by 'jagsUI', or an optional subset of parameter nodes.
#' @param x Output object returned from `jagsUI`
#' @param p Optional string to begin posterior names.  If `NULL` is used, all parameters will be used
#' @param exact Whether name must be an exact match (`TRUE`) or with initial sub-string matching only supplied characters (`FALSE`).
#' Defaults to `FALSE.`
#' @return A 2-dimensional correlation matrix (n X n, where n is the number of parameter nodes)
#' @seealso \link{plotcor_jags}
#' @author Matt Tyers
#' @examples
#' cor_jags(asdf_jags_out)
#' @export
cor_jags <- function(x, p=NULL, exact=FALSE) {
  if(!inherits(x,"jagsUI")) stop("Input must be an output object returned from jagsUI::jags().")
  df <- jags_df(x, p=p, exact=exact)   ### make this jags_df
  suppressWarnings(dfcor <- cor(df))
  return(dfcor)
}
