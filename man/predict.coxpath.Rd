\name{predict.coxpath}
\alias{predict.coxpath}
\title{
  Makes predictions at particular points along the fitted coxpath
}
\description{
  This function makes predictions at particular points along the fitted
  \code{coxpath.} The coefficients, log-partial-likelihood, linear
  predictor or the risk can be computed. A \code{coxph} object can be
  returned at one particular value of \eqn{\lambda.}
}
\usage{
  \method{predict}{coxpath}(object, data, s, type = c("coefficients", "loglik",
        "lp", "risk", "coxph"), mode = c("step",
        "norm.fraction", "norm", "lambda.fraction", "lambda"),
        eps = .Machine$double.eps, ...)
}
\arguments{
  \item{object}{
    a \code{coxpath} object
  }
  \item{data}{
    a list containing \code{x,} \code{time,} and \code{status,} with
    which the predictions are made. If \code{type=lp} or
    \code{type=risk,} then \code{x} is required. If \code{type=loglik}
    or \code{type=coxph,} then \code{x,} \code{time,} and \code{status}
    are required.
  }
  \item{s}{
    the values of \code{mode} at which the predictions are made. If
    \code{type=coxph,} only the first element of \code{s} is used. If
    \code{s} is missing, then the steps at which the active set changed
    are used, and thus, \code{mode} is automatically switched to
    \code{step.}
  }
  \item{type}{
    If \code{type=coefficients,} the coefficients are returned; if
    \code{type=loglik,} log-partial-likelihoods are returned; if
    \code{type=lp,} linear predictors (\eqn{x'\beta}) are returned; if
    \code{type=risk,} risks (\eqn{e^{x'\beta}}) are returned; and if
    \code{type=coxph,} a \code{coxph} object (as in \pkg{survival}
    package) at the first element of \code{s} is returned. (i.e. the
    components of a \code{coxph} object such as coefficients, variance,
    and the test statistics are adjusted to the shrinkage corresponding
    to \code{s.} A \code{coxph} object can be further used as an
    argument to the functions in \pkg{survival} package.) Default is
    \code{coefficients.} The coefficients for the initial input
    variables are returned (rather than the standardized coefficients).
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
  \item{eps}{
    an effective zero
  }
  \item{...}{
    other options for the prediction
  }
}
\references{
  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
cv.coxpath, coxpath, plot.coxpath
}
\examples{
data(lung.data)
attach(lung.data)
fit <- coxpath(lung.data)
pred.a <- predict(fit, x, s = seq(0, 1, length=10),
                  mode = "norm.fraction")
library(survival)
pred.b <- predict(fit, lung.data, s = 0.5, type = "coxph",
                  mode = "lambda.fraction")
pred.s <- survfit(pred.b)
plot(pred.s)
detach(lung.data)
}
\keyword{models}
\keyword{regression}
