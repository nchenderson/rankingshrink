GradDescUpdate <- function(par, X, A, nu, W, internal.obj, stplngth) {
  # In this case, the gradient is
  # At W expit(Abeta)(1 - expit(Abeta))

  if(internal.obj=="auc") {
     theta.hat <- as.numeric(X%*%par)/nu
     A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
     mu.tmp <- plogis(A.beta)
     wmu.tmp <- ww*as.vector(mu.tmp*(1 - mu.tmp))
     gradf <- crossprod(A, wmu.tmp)

     beta.tmp <- par - stplngth*gradf
     beta.new <- beta.tmp/sqrt(sum(beta.tmp*beta.tmp))
  } else {
     if(internal.obj=="gaussian") {
        gradh <- crossprod(X, X%*%par - y)
     } else if(internal.obj=="logistic") {
        X.beta <- as.numeric(X%*%par)
        gradh <- crossprod(X, y - plogis(X.beta))
     }
     theta.hat <- as.numeric(X%*%par)/nu
     A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
     mu.tmp <- plogis(A.beta)
     wmu.tmp <- ww*as.vector(mu.tmp*(1 - mu.tmp))
     gradf <- gradh + crossprod(A, wmu.tmp)

     beta.new <- par - stplngth*gradf
  }
  return(as.numeric(beta.new))
}
