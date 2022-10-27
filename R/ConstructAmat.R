
ConstructAmat <- function(X, nu) {
  ## A should be n^2 x p
  ee <- matrix(1, nrow=nrow(X), ncol=1)
  yy <- kronecker(X, ee) - kronecker(ee, X)
  A <- yy/nu
  return(A)
}





