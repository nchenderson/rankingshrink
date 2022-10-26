RankPenalizedObj <- function(par, y, X, A, V, nu,
                               internal.obj, stplngth) {
  n <- length(y)
  if(internal.obj=="auc") {
     stop("Nonparametric option not available yet")
  } else if(internal.obj=="gaussian") {
     X.beta <- X%*%par
     resids <- y - X.beta
     theta.hat <- as.numeric(X.beta)/nu
     A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
     mu.tmp <- plogis(A.beta)

     hI <- mean(X.beta*y) - mean(X.beta*X.beta)/2
     ans <- rep(0, 2)
     ans[1] <- hI
     ans[2] <- sum(V*mu.tmp)
     objfn.val <- hI + sum(V*mu.tmp)
  } else if(internal.obj=="logistic") {
     X.beta <- X%*%par
     theta.hat <- as.numeric(X.beta)/nu
     A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
     mu.tmp <- plogis(A.beta)

     hI <- sum(y*X.beta) - sum(log(1 + exp(X.beta)))
     objfn.val <- hI/n + sum(V*mu.tmp)
  }
  return(ans)
}
