\name{plot.coxpath}
\alias{plot.coxpath}
\title{
  Plots the regularization path computed from coxpath
}
\description{
  This function takes a \code{coxpath} object and visualizes the
  regularization path. The horizontal axis can be \code{norm,}
  \code{lambda} or \code{step.} The vertical axis can be
  \code{coefficients,} \code{aic} or \code{bic.}
}
\usage{
  \method{plot}{coxpath}(x, xvar = c("norm", "lambda", "step"),
     type = c("coefficients", "aic", "bic"),
     plot.all.steps = FALSE, xlimit = NULL, predictor = FALSE,
     omit.zero = TRUE, breaks = TRUE, mar = NULL, main = NULL,
     eps = .Machine$double.eps, ...)
}
\arguments{
  \item{x}{
    a \code{coxpath} object
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
  \item{plot.all.steps}{
    If \code{TRUE,} all the steps taken along the path are marked on the
    plot. If \code{FALSE,} which is the default, only the steps at which
    the active set changed are shown on the plot.
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
    If \code{TRUE,} the predictors that were never in the active set are
    omitted.
  }
  \item{breaks}{
    If \code{TRUE,} vertical lines are drawn at the points where the
    active set changes and numbered with the degrees of freedom.
  }
  \item{mar}{
    margin relative to the current font size
  }
  \item{main}{
    title of the plot
  }
  \item{eps}{
    an effective zero
  }
  \item{...}{
    other options for the plot
  }
}
\references{
  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
cv.coxpath, coxpath, predict.coxpath
}
\examples{
data(lung.data)
attach(lung.data)
fit <- coxpath(lung.data)
par(mfrow=c(3, 2))
plot(fit)
plot(fit, xvar="lambda")
plot(fit, xvar="step")
plot(fit, xvar="step", xlimit=8)
plot(fit, type="aic")
plot(fit, type="bic")
detach(lung.data)
}
\keyword{models}
\keyword{regression}
