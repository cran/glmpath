\name{summary.coxpath}
\alias{summary.coxpath}
\title{
  Produces an anova-type summary for a coxpath object
}
\description{
  This function produces an anova-type summary for a coxpath object.
}
\usage{
  \method{summary}{coxpath}(object, ...)
}
\arguments{
  \item{object}{
    a \code{coxpath} object
  }
  \item{...}{
    additional arguments
  }
}
\value{
  A data.frame is returned, with the following components at transition points:
  \item{Df}{
    degrees of freedom at each step
  }
  \item{Log.p.lik}{
    log-partial-likelihood at each step
  }
  \item{AIC}{
    AIC value at each step
  }
  \item{BIC}{
    BIC value at each step
  }
}
\details{
  An anova type of summary is returned, including Df,
  Log-partial-likelihood, AIC, and BIC values for the steps where the
  active set changed.
}
\references{
  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
coxpath, plot.coxpath, print.coxpath
}
\examples{
data(lung.data)
attach(lung.data)
fit <- coxpath(lung.data)
summary(fit)
detach(lung.data)
}
\keyword{models}
\keyword{regression}
