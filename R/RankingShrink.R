RankingShrink <- function(y, X, lambda, external.scores, internal.obj="gaussian",
                          discrepancy="spearman", nu=NULL, maxiter=500) {

  ## first standardize X
  Xnorm <- scale(X)
  na.cols <- apply(Xnorm, 2, function(x) sum(is.na(x)))
  intercept <- FALSE
  if(sum(na.cols > 0) == 1) {
      Xnorm[,na.cols > 0] <- rep(1, nrow(Xnorm))
      intercept <- TRUE
  } else if(sum(na.cols > 0) > 1) {
      stop("More than 1 column of X has no variation")
  }

  initial.pars <- GetInitialPars(y=y, X=Xnorm, internal.obj=internal.obj)
  if(is.null(nu)) {
    ## Create a method to calculate a default value of nu
    if(intercept) {
        beta.avg <- mean(initial.pars[-1])
    } else {
        beta.avg <- mean(initial.pars)
    }
    nu <- 0.4*beta.avg
  }

  p <- ncol(Xnorm)
  external.ranking <- rank(external.scores)
  Amat <- ConstructAmat(Xnorm, nu)
  Vmat <- ConstructVmat(y, external.ranking, lambda, discrepancy, internal.obj)
  ## Vmat is really a vector.
  Lf <- FindLipschitz(Xnorm, Amat, Vmat, internal.obj)
  stplngth <- 1/Lf
  ## it's probably safe, in practice, to use 8 times this step length
  stplngth <- 8*stplngth

  objfnvals <- rep(NA, maxiter + 1)
  beta.old <- initial.pars
  objfnvals[1] <- RankPenalizedObj(par=beta.old, y=y, X=Xnorm, A=Amat,
                                   V=Vmat, nu=nu, internal.obj=internal.obj,
                                   stplngth=stplngth)[1]
  penaltyvals <- rep(NA, maxiter + 1)
  penaltyvals[1] <- RankPenalizedObj(par=beta.old, y=y, X=Xnorm, A=Amat,
                                   V=Vmat, nu=nu, internal.obj=internal.obj,
                                   stplngth=stplngth)[2]
  for(k in 1:maxiter) {
      beta.new <- GradDescUpdate(par=beta.old, y=y, X=Xnorm, A=Amat,
                                 V=Vmat, nu=nu, internal.obj=internal.obj,
                                 stplngth=stplngth)
      tmp <- RankPenalizedObj(par=beta.new, y=y, X=Xnorm, A=Amat,
                                         V=Vmat, nu=nu, internal.obj=internal.obj,
                                         stplngth=stplngth)
      objfnvals[k+1] <- tmp[1]
      penaltyvals[k+1] <- tmp[2]
      beta.old <- beta.new
  }
  fitted.vals <- NULL
  if(internal.obj=="gaussian") {
      fitted.vals <- as.vector(Xnorm%*%beta.old)
  }
  return(list(coef=beta.old, objfnvals=objfnvals, fitted.values=fitted.vals,
              penaltyvals=penaltyvals))
}


