\name{RankingShrink}
\alias{RankingShrink}
\title{
Ranking-informed shrinkage estimation
}
\description{
Performs penalized regression using a penalty based rankings from an external model
}
\usage{
RankingShrink(y, X, lambda, external.scores, internal.obj="gaussian",
              discrepancy="spearman", nu=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the vector of responses
}
\item{X}{
the design matrix
}
\item{lambda}{
the nonnegative penalty term
}
\item{external.scores}{
The risk scores generated from an external risk model.
}
\item{internal.obj}{
The type of internal objective function to be optimized when fitting the regression model. Options are \code{gaussian}, \code{logistic}, or \code{auc}.
}
\item{discrepancy}{
Type of discrepancy measures used for the ranking-based penalty function. Options are \code{spearman}, or \code{kendall}.
}
\item{maxiter}{
The maximum number of gradient descent iterations to perform.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list with the following components
\item{coef}{The estimated regression coefficients}
\item{objfnval}{The value of the objective function at each iteration of the fitting procedure}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(diabetes)
diabetes.old <- diabetes[1:100,]
diabetes.new <- diabetes[101:442,]

prior.model <- lm(prog ~ age + sex + map, data=diabetes.old)

y <- diabetes.new$prog
X <- cbind(rep(1, nrow(diabetes.new)), as.matrix(diabetes.new[,1:10]))
Z <- model.matrix(prog ~ age + sex + map, data=diabetes.new)
external.scores <- as.vector(Z \%*\% prior.model$coefficients)

tst <- RankingShrink(y=y, X=X, lambda=0,
                     external.scores=external.scores,
                     optimization="optimg", maxiter=50)
rank_fitted <- tst$fitted.values
ols <- lm(y ~ X - 1)$fitted.values
plot(ols, rank_fitted)
cor(ols, rank_fitted)
cor(external.scores, rank_fitted)

\dontrun{
lam.seq <- exp(seq(0, 10, length.out=20))
cols.seq <- cext.seq <- rep(NA, 20)
for(k in 1:20) {
  tst <- RankingShrink(y=y, X=X, lambda=lam.seq[k],
                       external.scores=external.scores,
                       maxiter=50, optimization="optimg")
  rank_fitted <- tst$fitted.values

  cols.seq[k] <- cor(ols, rank_fitted, method="spearman")
  cext.seq[k] <- cor(external.scores, rank_fitted, method="spearman")
}

plot(lam.seq, cols.seq, type="n", ylim=c(0.5, 1),
     xlab="lambda", ylab="Rank Correlation",
     main="Correlation with Other Methods")
lines(lam.seq, cols.seq)
lines(lam.seq, cext.seq, col="red")
legend("bottomleft", legend=c("OLS", "External"),
       col=c("black","red"), bty='n', lwd=2)
}

}
\keyword{ regression }
\keyword{ robust }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
