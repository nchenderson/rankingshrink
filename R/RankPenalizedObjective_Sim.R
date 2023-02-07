RankPenalizedObj_Sim <- function(par, y, Xlist, V, nu,
                                 internal.obj, stplngth) {
  n <- length(y)
  S <- length(Xlist)
  if(internal.obj=="auc") {
    stop("Nonparametric option not available yet")
  } else if(internal.obj=="gaussian") {
    Gmat <- matrix(NA, nrow=n*n, ncol=S)
    ee <- matrix(1, nrow=nrow(X), ncol=1)

    for(s in 1:S) {
      ## form A
      AA_diff <- kronecker(Xlist[[s]], ee) - kronecker(ee, Xlist[[s]])
      A.beta <- AA_diff%*%beta
      Gmat[,s] <- plogis(A.beta/nu)
    }
    X.beta <- as.numeric(X%*%par)
    mu.tmp <- rowMeans(Gmat)

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
