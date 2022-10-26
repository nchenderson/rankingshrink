GetInitialPars <- function(y, X, internal.obj) {
   if(internal.obj=="gaussian") {
       initial.par <- lm(y ~ X - 1)$coefficients
   } else if(internal.obj=="logistic") {
       initial.par <- glm(y ~ X - 1, family="binomial")$coefficients
   } else if(internal.obj=="auc") {
       stop("nonparametric objective not available yet")
   }
   return(initial.par)
}
