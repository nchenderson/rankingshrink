

FindLipschitz <- function(X, A, V, internal.obj) {

  n <- nrow(X)
  if(internal.obj=="gaussian") {
    XtX <- crossprod(X, X)
    Lh <- norm(XtX, "2")/n
  } else if(internal.obj=="binary") {
    XtX <- crossprod(X, X)
    Lh <- norm(XtX, "2")/(4*n)
  } else if(internal.obj=="auc") {
    hval <- 0
  }
  AtA <- crossprod(A, A)
  lamA <- norm(AtA, "2")
  vmax <- max(abs(V))
  Lf <- Lh + 0.97*vmax*lamA
  return(Lf)
}
