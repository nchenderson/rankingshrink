rpenlm <- function(y, X, lambda, external.rankings, nu=0.2) {
     ## Drop intercept and standardize columns of X
     Amat <- ConstructAmat(X, 0.2)
     Vmat <- ConstructVmat(y, external.rankings, lambda=lambda, discrepancy="spearman",
                        internal.obj="gaussian")

     XtX <- crossprod(X, X)
     Xty <- crossprod(X, y)
     beta.init <- solve(XtX, Xty)
     beta.ls <- beta.init
     #rpen <- optimg(par=beta.init, fn=ContinuousObjFn, gr=ContinuousGradFn,
     #               y=y, X=X, Amat=Amat, Vmat=Vmat, XtX=XtX, Xty=Xty, Interval=1e-6,
     #               method="ADAM")
     #ans <- list()
     #ans$par <- rpen$par
     #ans$grad <- ContinuousGradFnShort(par=rpen$par, y=y, X=X, Amat=Amat, Vmat=Vmat,
     #                                  XtX=XtX, Xty=Xty)

     beta.old <- rep(0.1, length(beta.init))
     #beta.old <- beta.init
     numiters <- 100
     objfnvals <- rep(NA, numiters + 1)
     objfnvals[1] <- ContinuousObjFnNew(beta.old, y, X, Amat, Vmat, XtX, Xty, lambda)
     weights.pgamma <- rep(NA, nrow(Amat))
     for(k in 1:numiters) {
         Abeta <- as.numeric(Amat%*%beta.old)
         phat <- plogis(Abeta) # This is expit(phi)
         small.phi <- abs(Abeta) < 1e-4
         weights.pgamma[small.phi] <- 1/4 - (Abeta[small.phi]^2)/48 + (Abeta[small.phi]^4)/480
         weights.pgamma[!small.phi] <- ((phat[!small.phi] - 0.5))/Abeta[!small.phi]

         tmpQ <- Vmat*phat
         qweights <- tmpQ/sum(tmpQ)
         #Q <- diag(qweights)
         #QW <- diag(weights.pgamma*qweights)
         QWA <- weights.pgamma*qweights*Amat
         AtA <- crossprod(Amat, QWA)
         #AV <- rowSums(crossprod(Amat, Q))
         AV <- rowSums(t(Amat*qweights))
         beta.new <- solve(XtX + lambda*AtA, Xty + (lambda/2)*AV)

         objfnvals[k+1] <- ContinuousObjFnNew(beta.new, y, X, Amat, Vmat, XtX, Xty, lambda)
         beta.old <- as.numeric(beta.new)
     }
     ans <- list()
     ans$par <- beta.old
     ans$objfnvals <- objfnvals

     ## Now, compute degrees of freedom:
     Abeta <- as.numeric(Amat%*%beta.ls)
     weights.pgamma <- rep(NA, length(Abeta))
     phat <- plogis(Abeta) # This is expit(phi)
     small.phi <- abs(Abeta) < 1e-4
     weights.pgamma[small.phi] <- 1/4 - (Abeta[small.phi]^2)/48 + (Abeta[small.phi]^4)/480
     weights.pgamma[!small.phi] <- ((phat[!small.phi] - 0.5))/Abeta[!small.phi]

     tmpQ <- Vmat*phat
     qweights <- tmpQ/sum(tmpQ)
     #Q <- diag(qweights)
     #QW <- diag(weights.pgamma*qweights)
     QWA <- weights.pgamma*qweights*Amat
     AtA <- crossprod(Amat, QWA)
     Dmat <- solve(XtX + lambda*AtA, XtX)
     ans$df <- sum(diag(Dmat))
     return(ans)
}
