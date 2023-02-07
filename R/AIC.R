AIC <- function(par, y, X, A, V, internal.obj) {

    aic_val <- NULL
    aic_correct_val <- NULL

    if(internal.obj=="gaussian") {
      X.beta <- as.numeric(X%*%par)
      A.beta <- as.numeric(A%*%par)
      mu.tmp <- plogis(A.beta)

      XtX <- crossprod(X, X)
      ww <- V*mu.tmp*(1 - mu.tmp)*(1 - 2*mu.tmp)
      BB <- XtX + crossprod(A, ww*A)
      Phat <- solve(BB, XtX)
      df <- sum(diag(Phat))

      rss <- sum((y - X.beta)*(y - X.beta))
      aic_val <- log(rss) + (2*df)/nrow(X)
      aic_correct_val <- log(rss) + (2*df + 4)/(nrow(X) - df - 3)
    }
    return(list(aic=aic_val, aic_correct=aic_correct_val))
}

