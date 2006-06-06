\name{predict.glmpath}
\alias{predict.glmpath}
\title{
  Makes predictions at particular points along the fitted glmpath
}
\description{
  This function makes predictions at particular points along the fitted
  \code{glmpath.} The linear predictor, estimated response,
  log-likelihood, or the coefficients can be computed.
}
\usage{
  predict.glmpath(object, newx, newy, s, type = c("link", "response",
                  "loglik", "coefficients"), mode = c("step",
                  "norm.fraction", "norm", "lambda.fraction", "lambda"),
                  weight = NULL, eps = .Machine$double.eps, ...)
}
\arguments{
  \item{object}{
    a \code{glmpath} object
  }
  \item{newx}{
    a matrix of features at which the predictions are made. If
    \code{type=link,} \code{type=response,} or \code{type=loglik,}
    \code{newx} is required.
  }
  \item{newy}{
    a vector of responses corresponding to \code{newx.} If
    \code{type=loglik,} \code{newy} is required.
  }
  \item{s}{
    the values of \code{mode} at which the predictions are made
  }
  \item{type}{
    If \code{type=link,} the linear predictors are returned; if
    \code{type=response,} the estimated responses are returned; if
    \code{type=loglik,} the log-likelihoods are returned, and if
    \code{type=coefficients,} the coefficients are returned. The
    coefficients for the initial input variables are returned (rather
    than the standardized coefficients). Default is \code{link.}
  }
  \item{mode}{
    what \code{mode=s} refers to. If \code{mode=step,} \code{s} is the
    number of steps taken; if \code{mode=norm.fraction,} \code{s} is the
    fraction of the L1 norm of the standardized coefficients (with
    respect to the largest norm); if \code{mode=norm,} \code{s} is the
    L1 norm of the standardized coefficients; if
    \code{mode=lambda.fraction,} \code{s} is the fraction of
    log(\eqn{\lambda}); and if \code{mode=lambda,} \code{s} is
    \eqn{\lambda}. Default is \code{step.}
  }
  \item{weight}{
    an optional vector of weights for observations. \code{weight} is
    effective only if \code{type=loglik.}
  }
  \item{eps}{
    an effective zero
  }
  \item{...}{
    other options for the prediction
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
cv.glmpath, glmpath, plot.glmpath
}
\examples{
data(heart.data)
attach(heart.data)
fit <- glmpath(x, y, family=binomial)
pred <- predict(fit, x, s = seq(0, 1, length=10), mode="norm.fraction")
detach(heart.data)
}
\keyword{models}
\keyword{regression}
