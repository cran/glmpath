\name{summary.glmpath}
\alias{summary.glmpath}
\title{
  Produces an anova-type summary for a glmpath object
}
\description{
  This function produces an anova-type summary for a glmpath object.
}
\usage{
  \method{summary}{glmpath}(object, ...)
}
\arguments{
  \item{object}{
    a \code{glmpath} object
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
  \item{Deviance}{
    deviance computed at each step
  }
  \item{AIC}{
    AIC value at each step
  }
  \item{BIC}{
    BIC value at each step
  }
}
\details{
  An anova type of summary is returned, including Df, Deviance, AIC, and
  BIC values for the steps where the active set changed.
}
\references{
  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
glmpath, plot.glmpath, print.glmpath
}
\examples{
data(heart.data)
attach(heart.data)
fit <- glmpath(x, y)
summary(fit)
detach(heart.data)
}
\keyword{models}
\keyword{regression}
