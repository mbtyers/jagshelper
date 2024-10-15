library(jagshelper)
skeleton("asdf")
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



### actually starting on the function!
# define which data objects to withold data (char)
# define which model output params will correspond (char - will be ypp i think, or maybe y)
# how many folds (default to 5)
# model file

rmse <- function(x1, x2) sqrt(mean((x1-x2)^2, na.rm=TRUE))

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

kfold <- function(model.file, data,
                  pdata, pmodel,
                  k=5,
                  ...) {
  data_y <- data[[pdata]]     ### figure out if there will be a case for multiple vectors or matrices or something
  # fold <- sample(x=seq(k), size=length(data_y), replace=TRUE)  ## figure out a more robust way to do this!!
  fold <- allocate(n=length(data_y), k=k)
  pred_y <- rep(NA, length(data_y))

  if(interactive()) pb <- txtProgressBar(style=3)

  for(i_fold in seq(k)) {
    data_fold <- data
    data_fold[[pdata]][fold==i_fold] <- NA

    out_fold <- jagsUI::jags(model.file=model.file,
                             data=data_fold,
                             parameters.to.save=pmodel,
                             verbose=FALSE,     # maybe do a progress bar here
                             ...=...)
    pred_fold <- out_fold$q50[[pmodel]]

    pred_y[fold==i_fold] <- pred_fold[fold==i_fold]
    if(interactive()) setTxtProgressBar(pb=pb, value=i_fold/k)
  }
  return(rmse(x1=data_y, x2=pred_y))
}
kfold(pdata="y", pmodel="ypp", #try y here too
  model.file=asdf_jags, data=asdf_data,
      n.chains=ncores, parallel=T, n.iter=niter,
      n.burnin=niter/2, n.thin=niter/2000)   # might be able to get this stuff from mcmc.info

kk <- c(3,5,10,30,60)
nrep <- 50
rmsemat <- matrix(nrow=nrep, ncol=length(kk))
for(j in seq_along(kk)) {
  for(i in 1:nrep) {
    rmsemat[i,j] <- kfold(k=kk[j],
                          pdata="y", pmodel="y", #try y here too
                          model.file=asdf_jags, data=asdf_data,
                          n.chains=ncores, parallel=T, n.iter=niter,
                          n.burnin=niter/2, n.thin=niter/2000)
    print(c(j,i))
  }
}
boxplot(rmsemat)
for(j in seq_along(kk)) {
  for(i in 1:nrep) {
    rmsemat[i,j] <- kfold(k=kk[j],
                          pdata="y", pmodel="ypp", #try y here too
                          model.file=asdf_jags, data=asdf_data,
                          n.chains=ncores, parallel=T, n.iter=niter,
                          n.burnin=niter/2, n.thin=niter/2000)
    print(c(j,i))
  }
}
boxplot(rmsemat)
