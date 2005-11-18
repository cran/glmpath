\name{cv.glmpath}
\alias{cv.glmpath}
\title{
  Computes cross-validated prediction errors for glmpath
}
\description{  
  This function computes cross-validation errors for \code{glmpath.}
  0-1 loss is used for \code{family=binomial}; mean squared error loss
  is used for \code{family=gaussian} and \code{family=poisson.}
}
\usage{
  cv.glmpath(x, y, family = binomial, weight = rep(1, length(y)),
             nfold = 10, fraction = seq(from=0, to=1, length=100),
             mode = c("norm","lambda"), plot.it = TRUE, se = TRUE, ...)
}
\arguments{
\item{x}{
  matrix of predictors
}
\item{y}{
  response
}
\item{family}{
  name of a family function that represents the distribution of y to be
  used in the model. It must be \code{binomial}, \code{gaussian}, or
  \code{poisson}. For each one, the canonical link function is used;
  \code{logit} for binomial, \code{identity} for gaussian, and
  \code{log} for poisson distribution. Default is \code{binomial.}
}
\item{weight}{
  an optional vector of weights for observations
}
\item{nfold}{
  number of folds to be used in the cross-validation. Default is
  \code{nfold=10.}
}
\item{fraction}{
  fraction of the L1 norm or log(\eqn{\lambda}) with respect to their
  maximum values at which CV values are computed. Default is
  \code{seq(0,1,length=100).}
}
\item{mode}{
  If \code{mode=norm,} cross-validation is done at certain grids of L1
  norm. If \code{mode=lambda,} cross-validation is done at certain grids of
  log(\eqn{\lambda}). Default is \code{norm.}
}
\item{plot.it}{
  If \code{TRUE,} the CV curve is plotted.
}
\item{se}{
  If \code{TRUE,} the standard errors are computed.
}
\item{...}{
  other options for glmpath
}
}
\references{
  Mee Young Park and Trevor Hastie (2005) L1 Regularization Path
Algorithm for Generalized Linear Models - This will be soon available at the
authors' websites, \url{http://www.stanford.edu/~mypark} or
\url{http://stat.stanford.edu/~hastie/pub.htm}.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
glmpath, plot.glmpath, predict.glmpath
}
\examples{
data(heart.data)
attach(heart.data)
cv <- cv.glmpath(x, y, family=binomial)
detach(heart.data)
}
\keyword{models}
\keyword{regression}