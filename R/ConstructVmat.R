ConstructVmat <- function(y, external.ranking, lambda,
                        discrepancy, internal.obj) {
   ## The returned matrix should be n^2 x n^2
   ## We only need to return the diagonals of the matrix V

   ## Use Reduce for the marginalized matrix multiplication.
    # https://stackoverflow.com/questions/53285930/multiplication-of-many-matrices-in-r
   n <- length(y)
   if(discrepancy=="spearman" & internal.obj!="auc") {
      Vmat <- (4*n*n)*rep(external.ranking, each=n)
      #Vmat <- lambda/(4*n)*rep(n - external.ranking, each=n)
   } else if(discrepancy=="spearman" & internal.obj == "auc") {
      CompareY <- outer(y, y) > 0
      Vmat <- -lambda*matrix(rep(external.ranking, each=n), nrow=n, ncol=n)/(4*n) + CompareY
   } else if(discrepancy=="kendall" & internal.obj != "auc") {
      CompareRanks <- outer(external.ranking, external.ranking, FUN="-") > 0
      Vmat <- (2*lambda/(n-1))*(0.5 - CompareRanks)
   } else if(discrepancy=="kendall" & internal.obj == "auc") {
      CompareRanks <- outer(external.ranking, external.ranking, FUN="-") > 0
      CompareY <- outer(y, y) > 0
      Vmat <- (2*lambda/(n-1))*(0.5 - CompareRanks) + CompareY
   }
   return(c(Vmat))
}
