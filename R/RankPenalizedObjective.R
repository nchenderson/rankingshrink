RankPenalizedObj <- function(par, y, X, A, V, nu,
                               internal.obj, stplngth) {
  n <- length(y)
  if(internal.obj=="auc") {
     stop("Nonparametric option not available yet")
  } else if(internal.obj=="gaussian") {
     X.beta <- as.numeric(X%*%par)
     A.beta <- as.numeric(A%*%par)
     mu.tmp <- plogis(A.beta)

     hI <- mean(X.beta*X.beta)/2 - mean(X.beta*y)
     objfn.val <- hI + sum(V*mu.tmp)
  } else if(internal.obj=="logistic") {
     X.beta <- X%*%par
     theta.hat <- as.numeric(X.beta)/nu
     A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
     mu.tmp <- plogis(A.beta)

     hI <- sum(y*X.beta) - sum(log(1 + exp(X.beta)))
     objfn.val <- hI/n + sum(V*mu.tmp)
  }
  return(objfn.val)
}
