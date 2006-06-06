\name{cv.coxpath}
\alias{cv.coxpath}
\title{
  Computes cross-validated (minus) log-partial-likelihoods for coxpath
}
\description{  
  This function computes cross-validated (minus) log-partial-likelihoods
  for \code{coxpath.}
}
\usage{
  cv.coxpath(data, method = c("breslow", "efron"), nfold = 5,
             fraction = seq(from=0, to=1, length=100),
             mode = c("norm","lambda"), plot.it = TRUE, se = TRUE, ...)
}
\arguments{
  \item{data}{
    a list consisting of \code{x:} a matrix of features, \code{time:}
    the survival time, and \code{status:} censor status with 1 if died
    and 0 if censored. 
  }
  \item{method}{
    approximation method for tied survival times. Approximations derived
    by Breslow (1974) and Efron (1977) are available. Default is
    \code{breslow.}
  }
  \item{nfold}{
    number of folds to be used in cross-validation. Default is
    \code{nfold=5.}
  }
  \item{fraction}{
    fraction of the L1 norm or log(\eqn{\lambda}) with respect to their
    maximum values at which CV values are computed. Default is
    \code{seq(0,1,length=100).}
  }
  \item{mode}{
    If \code{mode=norm,} cross-validation is done at certain values of
    L1 norm. If \code{mode=lambda,} cross-validation is done at certain
    values of log(\eqn{\lambda}). Default is \code{norm.}
  }
  \item{plot.it}{
    If \code{TRUE,} the CV curve is plotted.
  }
  \item{se}{
    If \code{TRUE,} the standard errors are plotted.
  }
  \item{...}{
    other options for coxpath
  }
}
\references{
  Mee Young Park and Trevor Hastie (2006) L1 Regularization Path
  Algorithm for Generalized Linear Models - available at the authors'
  websites, \url{http://www.stanford.edu/~mypark} or
  \url{http://stat.stanford.edu/~hastie/pub.htm}.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
coxpath, plot.coxpath, predict.coxpath
}
\examples{
data(lung.data)
attach(lung.data)
cv <- cv.coxpath(lung.data)
detach(lung.data)
}
\keyword{models}
\keyword{regression}