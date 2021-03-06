\name{plot.bootpath}
\alias{plot.bootpath}
\title{
  Generates the histograms or the pairwise scatter plots of the
  bootstrap coefficients computed from bootstrap.path
}
\description{
  This function takes a \code{bootpath} object from
  \code{bootstrap.path} and generates the histograms or the pairwise
  scatter plots of the bootstrap coefficients.
}
\usage{
  \method{plot}{bootpath}(x, type = c("histogram", "pairplot"),
     mfrow = NULL, mar = NULL, ...)
}
\arguments{
  \item{x}{
    a \code{bootpath} object from \code{bootstrap.path.}
  }
  \item{type}{
    If \code{type=histogram,} the histograms of bootstrap coefficients
    for individual features are generated. The red vertical bar
    indicates the coefficient computed using the whole data. The thick
    bar at zero indicates the frequency of the zero coefficients. If
    \code{type=pairplot,} the pairwise scatter plots of the bootstrap
    coefficients are generated. The red solid dot indicates the pair of
    coefficients computed using the whole data. Default is
    \code{histogram.}
  }
  \item{mfrow}{
    determines the numbers of rows and columns of the histograms on a
    page. 2 rows are generated as a default.
  }
  \item{mar}{
    margin relative to the current font size
  }
  \item{...}{
    other options for the plot
  }
}
\details{
  Fitting \code{glmpath} or \code{coxpath} gives a series of solution
  sets with a varying size of the active set. Once we select an
  appropriate value of the regularization parameter, and, thus a set of
  coefficients, we may then validate the chosen coefficients through a
  bootstrap analysis. \code{plot.bootstrap} summarizes the bootstrap
  results by generating the histograms or the pairwise scatter plots of
  the bootstrap coefficients.
}
\references{
  Bradley Efron and Robert Tibshirani (1993) \emph{An Introduction to
    the Bootstrap} CHAPMAN & HALL/CRC, Boca Raton.

  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
bootstrap.path, coxpath, glmpath
}
\examples{
data(heart.data)
attach(heart.data)
bootstrap.a <- bootstrap.path(x, y, B=5)
plot(bootstrap.a)
plot(bootstrap.a, type="pairplot")
detach(heart.data)
}
\keyword{models}
\keyword{regression}
