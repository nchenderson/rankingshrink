GradDescUpdate <- function(par, X, A, nu, W, response, stplngth) {
  # In this case, the gradient is
  # At W expit(Abeta)(1 - expit(Abeta))

  if(response="gaussian") {
    gradh <- crossprod(X, X%*%par - y)
  } else if(response="binary") {
    X.beta <- as.numeric(X%*%par)
    gradh <- crossprod(X, y - plogis(X.beta))
  } else if(response="nonparametric") {
    gradh <- 0
  }

  theta.hat <- as.numeric(X%*%par)/nu
  A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
  mu.tmp <- plogis(A.beta)
  wmu.tmp <- ww*as.vector(mu.tmp*(1 - mu.tmp))
  gradf <- gradh + crossprod(A, wmu.tmp)

  beta.tmp <- par - stplngth*gradf
  beta.new <- beta.tmp/sqrt(sum(beta.tmp*beta.tmp))
  return(as.numeric(beta.new))
}
