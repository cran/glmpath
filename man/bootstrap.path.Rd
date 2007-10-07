\name{bootstrap.path}
\alias{bootstrap.path}
\title{
  Generates a set of bootstrap coefficients for either glmpath or
  coxpath
}
\description{
  This function generates a set of bootstrap coefficients for either
  \code{glmpath} or \code{coxpath.} For each bootstrap run, the
  regularization parameter may be determined based on either \code{aic}
  or \code{bic.}
  }
\usage{
  bootstrap.path(x, y, data, B, index=NULL, path=c("glmpath","coxpath"),
                 method=c("aic","bic"), trace=FALSE, ...)
}
\arguments{
  \item{x}{
    matrix of features
  }
  \item{y}{
    response 
  }
  \item{data}{
    a list of data components. If \code{path=glmpath,} \code{data}
    consists of \code{x:} a matrix of features and \code{y:}
    response. \code{data} is not needed if \code{x} and \code{y} are
    input separately. If \code{path=coxpath,} \code{data} must be
    provided, including \code{x:} a matrix of features, \code{time:} the
    survival time, and \code{status:} censor status with 1 if died and 0
    if censored.
  }
  \item{B}{
    number of bootstrap runs
  }
  \item{index}{
    matrix (\code{B} rows, \code{ncol(x)} columns) of bootstrap sample
    indices. Each row is a vector of indices for a bootstrap run. If
    \code{index=NULL,} the indices are randomly chosen.
  }
  \item{path}{
    Bootstrap coefficients for either \code{glmpath} or \code{coxpath}
    are computed. Default is \code{path=glmpath.} 
  }
  \item{method}{
    For each bootstrap run, the regularization parameter is determined
    based on either \code{aic} or \code{bic.} Default is \code{aic.}
  }
  \item{trace}{
    If \code{TRUE,} the number of bootstrap runs is printed out.
  }
  \item{...}{
    other options for \code{glmpath} or \code{coxpath}
  }
}
\value{
  \code{bootstrap.path} returns a \code{bootpath} object, which is a
  matrix (\code{B} by \code{ncol(x)}) of bootstrap
  coefficients. Coefficients computed from the whole data are stored as
  an attribute \code{coefficients.}
}
\details{
  Fitting \code{glmpath} or \code{coxpath} gives a series of solution
  sets with a varying size of the active set. Once we select an
  appropriate value of the regularization parameter, and thus a set of
  coefficients, we may then validate the chosen coefficients through
  bootstrap analysis. \code{plot.bootstrap} summarizes the bootstrap
  results by generating the histograms or the pair scatter plots of the
  bootstrap coefficients.
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
coxpath, glmpath, plot.bootstrap
}
\examples{
data(heart.data)
attach(heart.data)
bootstrap.a <- bootstrap.path(x, y, B=100, method="bic")
detach(heart.data)
data(lung.data)
attach(lung.data)
bootstrap.b <- bootstrap.path(data=lung.data, B=100, path="coxpath")
detach(lung.data)
}
\keyword{models}
\keyword{regression}
