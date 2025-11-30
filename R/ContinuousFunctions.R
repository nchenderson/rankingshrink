ContinuousObjFn <- function(par, y, X, Amat, Vmat, XtX, Xty) {
  n <- nrow(X)
  p1 <- (1/2.0)*crossprod(y - X%*%par, y - X%*%par)
  p2 <- sum(Vmat*plogis(Amat%*%par))
  ans <- p1 + p2
  return(ans)
}

ContinuousObjFnNew <- function(par, y, X, Amat, Vmat, XtX, Xty, lambda) {
  n <- nrow(X)
  p1 <- (1/2.0)*crossprod(y - X%*%par, y - X%*%par)
  p2 <- -lambda*log(sum(Vmat*plogis(Amat%*%par)))
  ans <- p1 + p2
  return(ans)
}

ContinuousGradFn <- function(objfn, par, y, X, Amat, Vmat, XtX, Xty, Interval) {
  n <- nrow(X)
  mu <- plogis(Amat%*%par)
  ans <- XtX%*%par - Xty + crossprod(Amat, Vmat*mu*(1 - mu))
  return(as.numeric(ans))
}

ContinuousGradFnShort <- function(par, y, X, Amat, Vmat, XtX, Xty) {
  n <- nrow(X)
  mu <- plogis(Amat%*%par)
  ans <- (1/n)*XtX%*%par - (1/n)*Xty + crossprod(Amat, Vmat*mu*(1 - mu))
  return(as.numeric(ans))
}

ContinuousHessianShort <- function(par, y, X, Amat, Vmat, XtX, Xty) {
  n <- nrow(X)
  mu <- plogis(Amat%*%par)
  ans <- (1/n)*XtX #+ crossprod(Amat, Vmat*mu*(1 - mu))
  return(as.numeric(ans))
}
