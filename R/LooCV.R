LooCV <- function(y, X, lambda.seq, external.rankings, nu=0.2) {
    nint <- length(y)
    nlam <- length(lambda.seq)
    loo_score <- rep(NA, nlam)

    for(k in 1:nlam) {
       inner_loo <- rep(0, nint)
       for(i in 1:nint) {
          ymini <- y[-i]
          Xmini <- X[-i,]
          yi <- y[i]
          Xi <- X[i,]

          beta.hat <- rpenlm(y=ymini, X=Xmini, lambda=lambda.seq[k], external.rankings=external.rankings[-i])$par
          rpenfit <- sum(Xi*beta.hat)
          inner_loo[i] <- (yi - rpenfit)^2
       }
       loo_score[k] <- mean(inner_loo)
    }
    return(loo_score)
}

