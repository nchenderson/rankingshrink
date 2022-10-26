
ConstructAmat <- function(X, nu) {
  ## Probably want to have a separate type of computation
  ##  if n is not much bigger than p

  ## Also note that if n is quite large (more than 10,000) or so
  ## we probably don't want to allocate space for the "full" matrix
  p <- ncol(X)
  A <- NULL
  count <- 1
  for(k in 1:p) {
     A <- cbind(A, c(t(outer(X[,k], X[,k], FUN="-"))/nu))
  }
  return(A)
}





