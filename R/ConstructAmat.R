
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

#ConstructA2 <- function(X, nu) {
#  n <- nrow(X)
#  A <- matrix(0, n*n, ncol(X))
#  count <- 1
#  for(i in 1:n) {
#    for(j in 1:n) {
#      A[count,] <- (X[i,] - X[j,])/nu
#      count <- count + 1
#    }
#  }
##  return(A)
#}
#n <- 20000
#p <- 15
#X <- matrix(rnorm(n*p), nrow=n, ncol=p)
#system.time(A1 <- ConstructA(X, nu=2))
#system.time(A2 <- ConstructA2(X, nu=2))
#all.equal(A1, A2)





