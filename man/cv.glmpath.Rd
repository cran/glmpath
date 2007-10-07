\name{cv.glmpath}
\alias{cv.glmpath}
\title{
  Computes cross-validated (minus) log-likelihoods or prediction errors
  for glmpath
}
\description{  
  This function computes cross-validated (minus) log-likelihoods or
  prediction errors for \code{glmpath.}
}
\usage{
  cv.glmpath(x, y, data, family = binomial, weight = rep(1, n),
             offset = rep(0, n), nfold = 10,
             fraction = seq(from=0, to=1, length=100),
             type = c("loglik","response"), mode = c("norm","lambda"),
             plot.it = TRUE, se = TRUE, ...)
}
\arguments{
  \item{x}{
    matrix of features
  }
  \item{y}{
    response
  }
  \item{data}{
    a list consisting of \code{x:} a matrix of features and \code{y:}
    response. \code{data} is not needed if above \code{x} and \code{y}
    are input separately.
  }
  \item{family}{
    name of a family function that represents the distribution of y to
    be used in the model. It must be \code{binomial}, \code{gaussian},
    or \code{poisson}. For each one, the canonical link function is
    used; \code{logit} for binomial, \code{identity} for gaussian, and
    \code{log} for poisson distribution. Default is \code{binomial.}
  }
  \item{weight}{
    an optional vector of weights for observations
  }
  \item{offset}{
    an optional vector of offset. If a column of \code{x} is used as
    offset, the corresponding column must be excluded from \code{x.}
  }
  \item{nfold}{
    number of folds to be used in cross-validation. Default is
    \code{nfold=10.}
  }
  \item{fraction}{
    the fraction of L1 norm or log(\eqn{\lambda}) with respect to their
    maximum values at which the CV errors are computed. Default is
    \code{seq(0,1,length=100).}
  }
  \item{type}{
    If \code{type=loglik,} cross-validated minus log-likelihoods are
    computed. If \code{type=response,} cross-validated prediction errors
    are computed. Default is \code{loglik.}
  }
  \item{mode}{
    If \code{mode=norm,} cross-validation is run at certain values of
    L1 norm. If \code{mode=lambda,} cross-validation is run at certain
    values of log(\eqn{\lambda}). Default is \code{norm.}
  }
  \item{plot.it}{
    If \code{TRUE,} CV curve is plotted.
  }
  \item{se}{
    If \code{TRUE,} standard errors are plotted.
  }
  \item{...}{
    other options for glmpath
  }
}
\references{
  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
glmpath, plot.glmpath, predict.glmpath
}
\examples{
data(heart.data)
attach(heart.data)
cv.a <- cv.glmpath(x, y, family=binomial)
cv.b <- cv.glmpath(x, y, family=binomial, type="response")
detach(heart.data)
}
\keyword{models}
\keyword{regression}
