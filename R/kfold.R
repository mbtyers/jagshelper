library(jagshelper)
library(jagsUI)

# specify model, which is written to a temporary file
asdf_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    ypp[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0 + b1*x[i] + a[grp[i]]
  }

  for(j in 1:ngrp) {
    a[j] ~ dnorm(0, tau_a)
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  tau_a <- pow(sig_a, -2)
  sig_a ~ dunif(0, 10)
}', file=asdf_jags)


# simulate data to go with the example model
n <- 60
x <- rnorm(n, sd=3)
grp <- sample(1:3, n, replace=T)
y <- rnorm(n, mean=grp-x)

# bundle data to pass into JAGS
asdf_data <- list(x=x,
                  y=y,
                  n=length(x),
                  grp=as.numeric(as.factor(grp)),
                  ngrp=length(unique(grp)))

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  asdf_jags_out <- jagsUI::jags(model.file=asdf_jags, data=asdf_data,
                                parameters.to.save=c("b0","b1","sig","a","sig_a"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(asdf_jags_out)
plotRhats(asdf_jags_out)
traceworstRhat(asdf_jags_out, parmfrow = c(3, 3))




#### test case where y is a matrix

asdf_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    for(j in 1:ngrp) {
      y[i,j] ~ dnorm(mu[i,j], tau)
      mu[i,j] <- b0 + b1*x[i,j] + a[j]
    }

  }

  for(j in 1:ngrp) {
    a[j] ~ dnorm(0, tau_a)
  }

  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  tau_a <- pow(sig_a, -2)
  sig_a ~ dunif(0, 10)
}', file=asdf_jags)


# simulate data to go with the example model
n <- 60
x <- matrix(rnorm(n, sd=3),
            nrow=20, ncol=3)
# grp <- sample(1:3, n, replace=T)
y <- matrix(rnorm(n, mean=rep(1:3, each=20)-x),
            nrow=20, ncol=3)

# bundle data to pass into JAGS
asdf_data <- list(x=x,
                  y=y,
                  n=nrow(x),
                  # grp=as.numeric(as.factor(grp)),
                  ngrp=ncol(x))

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  asdf_jags_out <- jagsUI::jags(model.file=asdf_jags, data=asdf_data,
                                parameters.to.save=c("b0","b1","sig","a","sig_a","mu"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(asdf_jags_out)
plotRhats(asdf_jags_out)
traceworstRhat(asdf_jags_out, parmfrow = c(3, 3))






### actually starting on the function!
# define which data objects to withold data (char)
# define which model output params will correspond (char - will be ypp i think, or maybe y)
# how many folds (default to 5)
# model file

rmse <- function(x1, x2) sqrt(mean((x1-x2)^2, na.rm=TRUE))
mae <- function(x1, x2) mean(abs(x1-x2), na.rm=TRUE)

allocate <- function(n, k) {
  if(k < n) {
    grp <- rep(NA, n)
    i <- 1
    resample <- function(x, ...) x[sample.int(length(x), ...)]  # stolen from ?sample
    while(!all(!is.na(grp))) {
      grp[resample(which(is.na(grp)), 1)] <- i
      i <- i+1
      if(i>k) i <- 1
    }
  } else {
    grp <- seq(n)
  }
  return(grp)
}

# return what:  ---> i think all of these!!
# - just RMSE, just MAE
# - prediction vector (from medians)
# - MCMC matrix? maybe this is more fair
# - also, MCMC of specified parameters, for the purpose of lppd

# think if there will be a use case for multiple data vectors  ---> I think just one

# what if there is a data matrix (not vector)?

# document etc

# any difference between y and ypp? assuming there won't be


## to do:
# + structure output object with $rms_pred, $ma_pred, $pred
#   - also $postpred if post=TRUE ?
#   - posteriors for addl parameters if !is.null(addl_p)??? - dims will be weird
# + figure out how to handle when data_y is a matrix
#   - add optional fold_byrow and fold_bycolumn
#   - fold & allocate should work ok as-is when allocation is random
#   ====> I feel like i should be able to generalize to an array

#   - need a test case where data is matrix


# # trying a generalization thingy
# data_y <- array(dim=c(2,3,4))
# fold_dims <- 2:1
# rpt_dims <- (1:length(dim(data_y)))[-fold_dims]
# # rpt_bydim <- dim(data_y)
# # rpt_bydim[fold_dims] <- 1
# nfold <- prod(dim(data_y)[fold_dims])
# fold <- aperm(a = array(1:nfold,
#                         dim=c(dim(data_y)[fold_dims], dim(data_y)[rpt_dims])),
#               perm = order(c(fold_dims, rpt_dims)))
# fold


#' Automated K-fold or Leave One Out Cross Validation
#' @description Runs k-fold or Leave One Out Cross Validation for a specified
#' component of a JAGS data object, for a specified JAGS model.
#'
#' JAGS is run internally `k` times (or alternately, the size of the dataset),
#' withholding each of `k` "folds" of the input data and drawing posterior predictive
#' samples corresponding to the withheld data, which can then be compared to the
#' input data to assess model predictive power.
#'
#' Global measures of predictive power are provided in output: Root Mean Square
#' (Prediction) Error and Mean Absolute (Prediction) Error.  However, it is likely
#' that these measures will not be meaningful by themselves; rather, as a metric
#' for scoring a set of candidate models.
#' @param df Output object returned from `jagsUI::jags()`; or alternately,
#' two-dimensional `data.frame` or matrix in which parameter node element is
#' given by column and MCMC iteration is given by row.
#' @param model.file Path to file containing the model written in BUGS code,
#' passed directly to \link[jagsUI]{jags}.
#' @param data The named list of data objects,
#' passed directly to \link[jagsUI]{jags}.
#' @param p The name of the data object to use for K-fold or LOO CV.
#' @param addl_p Names of additional parameters to save from JAGS output,
#' if a metric such as Log Pointwise Predictive Density is to be calculated from
#' cross-validation results.  Defaults to `NULL`, indicating no additional parameters.
#' @param save_postpred Whether to save all posterior predictive samples,
#' in addition to posterior medians.  Defaults to `FALSE`.
#' @param k How many folds to use for cross-validation.  Defaults to `10`.
#' If this is set to a number equal to (or greater than) the sample size, LOOCV
#' behavior will result.
#' @param loocv Whether to perform Leave One Out (rather than k-fold) Cross
#' Validation.  Setting this to `TRUE` will override the input to `k=`.  Defaults
#' to `FALSE`.
#' @param fold_dims A vector of margins to use for selecting folds, if the data
#' object used for cross validation is a matrix or array.  For example, if the
#' data consists of a two-dimensional matrix, setting `fold_dims=1` will result
#' in whole rows being selected in each fold, or setting `fold_dims=2` will result
#' in whole columns.  However, this is generalized to accept vectors of
#' multiple `fold_dims` and higher-dimensional arrays of data.
#' @param ... additional arguments to \link[jagsUI]{jags}.  These may (or must)
#' include `n.chains`, `n.iter`, `n.burnin`, `n.thin`, `parallel`, etc.
#' @return A named list, which may consist of the following:
#' * `$pred_y`: Point estimates of predicted values corresponding to each data
#' element, calculated as the posterior predictive median value
#' * `$data_y`: Original data used for cross validation
#' * `$postpred_y`: All posterior predictive samples corresponding to each data
#' element, if `save_postpred=TRUE`
#' * `$rmse_pred`: Root Mean Square (Prediction) Error
#' * `$mae_pred`: Mean Absolute (Prediction) Error
#' * `$addl_p`: A list with length equal to `k` (or the number of folds), with
#' each list element containing all posterior samples for additional parameters,
#' if these are supplied in argument `addl_p=`.
#' * `$fold`: A vector, matrix, or array corresponding to the original data,
#' giving the numerical values of the corresponding fold used
#' @seealso \link{overlayenvelope}, \link{caterpillar}
#' @author Matt Tyers
#' @examples
#' ## I need examples
#' @export
kfold <- function(model.file, data,
                  p, addl_p=NULL, save_postpred=FALSE,
                  k=10, loocv=FALSE,
                  fold_dims=NULL,
                  # fold_byrow=FALSE, fold_bycolumn=FALSE,    ### or fold_by = c(NA,"row","column")
                  ...) {
  if(length(p) > 1) stop("Only one data object or parameter may be used at once.")

  data_y <- data[[p]]

  if(is.na(k) | loocv) {
    k <- length(data_y)
  }
  fold_dims <- fold_dims[fold_dims <= length(dim(data_y))]

  ## generalize fold_dim?
  # if((!fold_byrow & !fold_bycolumn) | is.null(dim(data_y)) | min(dim(data_y)==1)) {
  #   fold <- allocate(n=length(data_y), k=k)
  #   if(length(dim(data_y))==2) {
  #     fold <- matrix(fold, nrow=nrow(data_y), ncol=ncol(data_y))
  #   }
  # } else {    ## see if i can generalize to array in this section
  #   if(fold_byrow) {
  #     fold <- matrix(allocate(n=nrow(data_y), k=k),
  #                    nrow=nrow(data_y), ncol=ncol(data_y))
  #   }
  #   if(fold_bycolumn) {   # should make sure they can't both be true
  #     fold <- matrix(allocate(n=ncol(data_y), k=k),
  #                    nrow=nrow(data_y), ncol=ncol(data_y),
  #                    byrow=TRUE)
  #   }
  # }
  if(is.null(dim(data_y)) | min(dim(data_y)==1)) {
    # if data_y is a vector

    fold <- allocate(n=length(data_y), k=k)
  } else {

    # if data_y is a matrix or array
    if(is.null(fold_dims)) {
      fold <- array(allocate(n=length(data_y), k=k), dim=dim(data_y))
    } else {
      rpt_dims <- (1:length(dim(data_y)))[-fold_dims]   # dims where fold is repeated
      nfold <- prod(dim(data_y)[fold_dims])
      fold <- aperm(a = array(allocate(n=nfold, k=k),
                              dim=c(dim(data_y)[fold_dims], dim(data_y)[rpt_dims])),
                    perm = order(c(fold_dims, rpt_dims)))
    }
  }

  pred_y <- NA*data_y

  if(!is.null(addl_p)) {
    addl_p_post <- list()
  }

  if(interactive()) pb <- txtProgressBar(style=3)

  for(i_fold in seq(k)) {
    data_fold <- data
    data_fold[[p]][fold==i_fold] <- NA

    out_fold <- jagsUI::jags(model.file=model.file,
                             data=data_fold,
                             parameters.to.save=c(p, addl_p),
                             verbose=FALSE,
                                                       # n.chains=ncores, parallel=T, n.iter=niter,
                                                       # n.burnin=niter/2, n.thin=niter/2000,
                             codaOnly = FALSE, bugs.format = FALSE,
                             ...=...
                             )
    pred_fold <- out_fold$q50[[p]]
    pred_y[fold==i_fold] <- pred_fold[fold==i_fold]

    if(save_postpred) {
      if(i_fold==1) {   # initialize with the same dims as post pred
        postpred_y <- NA*out_fold$sims.list[[p]]
      }
      for(irep in 1:dim(postpred_y)[1]) {  # this is a hack
        # this is hilarious!!  can go as many commas as I want, but would be better to generalize
        if(length(dim(postpred_y))==2) {
          postpred_y[irep,][fold==i_fold] <- out_fold$sims.list[[p]][irep,][fold==i_fold]
        }
        if(length(dim(postpred_y))==3) {
          postpred_y[irep,,][fold==i_fold] <- out_fold$sims.list[[p]][irep,,][fold==i_fold]
        }
        if(length(dim(postpred_y))==4) {
          postpred_y[irep,,,][fold==i_fold] <- out_fold$sims.list[[p]][irep,,,][fold==i_fold]
        }
        if(length(dim(postpred_y))==5) {
          postpred_y[irep,,,,][fold==i_fold] <- out_fold$sims.list[[p]][irep,,,,][fold==i_fold]
        }
        if(length(dim(postpred_y))==6) {
          postpred_y[irep,,,,,][fold==i_fold] <- out_fold$sims.list[[p]][irep,,,,,][fold==i_fold]
        }
      }
    }

    if(!is.null(addl_p)) {
      addl_p_post[[i_fold]] <- out_fold$sims.list[addl_p]
    }

    if(interactive()) setTxtProgressBar(pb=pb, value=i_fold/k)
  }
  out <- list(pred_y=pred_y, data_y=data_y)
  if(save_postpred) {
    out$postpred_y <- postpred_y
  }
  out$rmse_pred <- rmse(x1=data_y, x2=pred_y)
  out$mae_pred <- mae(x1=data_y, x2=pred_y)
  if(!is.null(addl_p)) {
    out$addl_p <- addl_p_post
  }
  out$fold <- fold
  return(out)
}
aa <- kfold(p="y",
            save_postpred = TRUE,
            addl_p = c("a","sig_a"),
  model.file=asdf_jags, data=asdf_data,
      n.chains=ncores, parallel=T, n.iter=niter,
      n.burnin=niter/2, n.thin=niter/2000)
str(aa)





# kk <- c(3,5,10,30,60)
# nrep <- 50
# rmsemat <- matrix(nrow=nrep, ncol=length(kk))
# for(j in seq_along(kk)) {
#   for(i in 1:nrep) {
#     rmsemat[i,j] <- kfold(k=kk[j],
#                           pdata="y", pmodel="y", #try y here too
#                           model.file=asdf_jags, data=asdf_data,
#                           n.chains=ncores, parallel=T, n.iter=niter,
#                           n.burnin=niter/2, n.thin=niter/2000)
#     print(c(j,i))
#   }
# }
# par(mfrow=c(1,1))
# boxplot(rmsemat)
# nrep <- 25
# for(j in seq_along(kk)) {
#   for(i in 1:nrep) {
#     rmsemat[i,j] <- kfold(k=kk[j],
#                           pdata="y", pmodel="ypp", #try y here too
#                           model.file=asdf_jags, data=asdf_data,
#                           n.chains=ncores, parallel=T, n.iter=niter,
#                           n.burnin=niter/2, n.thin=niter/2000)
#     print(c(j,i))
#   }
# }
# boxplot(rmsemat)
