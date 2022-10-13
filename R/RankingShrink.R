RankingShrink <- function(y, X, lambda, internal.obj="gaussian",
                          discrepancy="spearman", nu=NULL) {

  if(is.null(nu)) {
      ## Create a method to calculate a default value of nu
      nu <- 1
  }
  p <- ncol(X)
  Amat <- ConstructA(X, nu)
  Vmat <- ConstructVmat(y, external.ranking, lambda, discrepancy, internal.obj)
  Lf <- FindLipschitz(X, Amat, Vmat, internal.obj)
  stplngth <- 1/Lf
  ## it's probably safe, in practice, to use 8 times this step length
  stplngth <- 8*stplngth

  ## Finding initial parameters
  ## Write a function called GetInitialPars to get initial parameters
  lm.mod <- lm(y ~ X - 1)
  par.init <- lm.mod$coefficients/sqrt(sum(lm.mod$coefficients*lm.mod$coefficients))

  ans <- fpiter(par=par.init, fixptfn=GradDescUpdate, objfn=RankObjFn,
                X=X, A=Amat, nu=nu, W=Wmat, stplngth=stplngth,
                control=list(maxiter=5))
  return(ans)
}
