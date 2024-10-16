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
  grp <- rep(NA, n)
  i <- 1
  resample <- function(x, ...) x[sample.int(length(x), ...)]  # stolen from ?sample
  while(!all(!is.na(grp))) {
    grp[resample(which(is.na(grp)), 1)] <- i
    i <- i+1
    if(i>k) i <- 1
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

kfold <- function(model.file, data,
                  p, addl_p=NULL, save_postpred=FALSE,
                  k=5,
                  fold_byrow=FALSE, fold_bycolumn=FALSE,    ### or fold_by = c(NA,"row","column")
                  ...) {
  data_y <- data[[p]]     ### figure out if there will be a case for multiple vectors or matrices or something
  # fold <- sample(x=seq(k), size=length(data_y), replace=TRUE)  ## figure out a more robust way to do this!!

  if((!fold_byrow & !fold_bycolumn) | is.null(dim(data_y)) | min(dim(data_y)==1)) {
    fold <- allocate(n=length(data_y), k=k)
  } else {
    if(fold_byrow) {
      fold <- matrix(allocate(n=nrow(data_y), k=k),
                     nrow=nrow(data_y), ncol=ncol(data_y))
    }
    if(fold_bycolumn) {   # should make sure they can't both be true
      fold <- matrix(allocate(n=ncol(data_y), k=k),
                     nrow=nrow(data_y), ncol=ncol(data_y),
                     byrow=TRUE)
    }
  }

  pred_y <- NA*data_y #rep(NA, length(data_y))

  if(interactive()) pb <- txtProgressBar(style=3)

  for(i_fold in seq(k)) {
    data_fold <- data
    data_fold[[p]][fold==i_fold] <- NA

    out_fold <- jagsUI::jags(model.file=model.file,
                             data=data_fold,
                             parameters.to.save=c(p,addl_p),
                             verbose=FALSE,
                                                       n.chains=ncores, parallel=T, n.iter=niter,
                                                       n.burnin=niter/2, n.thin=niter/2000,
                             # ...=...
                             )
    pred_fold <- out_fold$q50[[p]]
    pred_y[fold==i_fold] <- pred_fold[fold==i_fold]

    if(save_postpred) {
      if(i_fold==1) {   # initialize with the same dims as post pred
        postpred_y <- NA*out_fold$sims.list[[p]]
      }
      for(irep in 1:dim(postpred_y)[1]) {  # this is a hack
        # actually this only works for 2d y
        postpred_y[irep,,][fold==i_fold] <- out_fold$sims.list[[p]][irep,,][fold==i_fold]
      }
    }

    if(interactive()) setTxtProgressBar(pb=pb, value=i_fold/k)
  }
  out <- list(pred_y=pred_y, data_y=data_y)
  if(save_postpred) {
    # out$postpred <- y_postpred
  }
  out$rmse_pred <- rmse(x1=data_y, x2=pred_y)
  out$mae_pred <- mae(x1=data_y, x2=pred_y)
  return(out)
}
aa <- kfold(p="y",
  model.file=asdf_jags, data=asdf_data,
      n.chains=ncores, parallel=T, n.iter=niter,
      n.burnin=niter/2, n.thin=niter/2000)   # might be able to get this stuff from mcmc.info
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
