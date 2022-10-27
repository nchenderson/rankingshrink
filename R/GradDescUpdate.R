GradDescUpdate <- function(par, y, X, A, V, nu,
                           internal.obj, stplngth) {
  n <- length(y)
  if(internal.obj=="auc") {
     theta.hat <- as.numeric(X%*%par)/nu
     A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
     mu.tmp <- plogis(A.beta)
     vmu.tmp <- V*as.vector(mu.tmp*(1 - mu.tmp))
     gradf <- crossprod(A, wmu.tmp)

     beta.tmp <- par - stplngth*gradf
     beta.new <- beta.tmp/sqrt(sum(beta.tmp*beta.tmp))
  } else {
     if(internal.obj=="gaussian") {
        gradh <- crossprod(X, X%*%par - y)/n
     } else if(internal.obj=="logistic") {
        X.beta <- as.numeric(X%*%par)
        gradh <- crossprod(X, y - plogis(X.beta))/n
     }
     theta.hat <- as.numeric(X%*%par)/nu
     A.beta <- A%*%par
     mu.tmp <- plogis(A.beta)
     vmu.tmp <- V*as.vector(mu.tmp*(1 - mu.tmp))
     gradf <- gradh + crossprod(A, vmu.tmp)

     beta.new <- par - stplngth*gradf
  }
  return(as.numeric(beta.new))
}
