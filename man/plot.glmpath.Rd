\name{plot.glmpath}
\alias{plot.glmpath}
\title{
  Plots the entire regularization path computed from glmpath
}
\description{
  This function takes a \code{glmpath} object and visualizes the
  regularization path. The horizontal axis can be \code{norm,}
  \code{lambda} or \code{step.} The vertical axis can be
  \code{coefficients,} \code{aic} or \code{bic.}
}
\usage{
  plot.glmpath(x, xvar = c("norm", "lambda", "step"),
               type = c("coefficients", "aic", "bic"), xlimit = NULL,
               predictor = FALSE, omit.zero = TRUE, breaks = TRUE,
               mar = NULL, eps = .Machine$double.eps, main = NULL, ...)
}
\arguments{
  \item{x}{
    a \code{glmpath} object
  }
  \item{xvar}{
    horizontal axis. \code{xvar=norm} plots against the L1 norm of the
    coefficients (to which L1 norm penalty was applied);
    \code{xvar=lambda} plots against \eqn{\lambda}; and \code{xvar=step}
    plots against the number of steps taken. Default is \code{norm.}
  }
  \item{type}{
    type of the plot, or the vertical axis. Default is
    \code{coefficients.}
  }
  \item{xlimit}{
    When the user wants to visualize a (beginning) sub-part of the plot,
    \code{xlimit} sets an upper limit to the L1 norm or the number of
    steps, or a lower limit to \eqn{\lambda}.
  }
  \item{predictor}{
    If \code{TRUE} and \code{type=coefficients,} the predictor step
    estimates are connected with dotted lines. If \code{FALSE,} only the
    corrector step estimates are connected with solid lines.
  }
  \item{omit.zero}{
    If \code{TRUE} and \code{type=coefficients,} the predictors that
    were never in the active set are omitted.
  }
  \item{breaks}{
    If \code{TRUE,} vertical lines are drawn at the points where the
    active set changes and numbered with the degrees of freedom.
  }
  \item{mar}{
    margin relative to the current font size
  }
  \item{eps}{
    an effective zero
  }
  \item{main}{
    title of the plot
  }
  \item{...}{
    other options for the plot
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
cv.glmpath, glmpath, predict.glmpath
}
\examples{
data(heart.data)
attach(heart.data)
fit <- glmpath(x, y, family=binomial)
par(mfrow=c(3,2))
plot(fit)
plot(fit,xvar="lambda")
plot(fit,xvar="step")
plot(fit,xvar="step",xlimit=8)
plot(fit,type="aic")
plot(fit,type="bic")
detach(heart.data)
}
\keyword{models}
\keyword{regression}
