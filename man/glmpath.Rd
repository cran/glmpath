\name{glmpath}
\alias{glmpath}
\title{
  Fits the entire L1 regularization path for generalized linear models
}
\description{
  This algorithm uses predictor-corrector method to compute the
  entire regularization path for generalized linear models with L1
  penalty.
}
\usage{
  glmpath(x, y, family = binomial, weight = rep(1,length(y)),
          max.steps = NULL, max.norm = 100*ncol(x),
          min.lambda = 0, max.arclength = Inf,
          standardize = TRUE, add.newvars = 1,
          bshoot.threshold = 0.5, relax.lambda = 1e-9,
          eps = .Machine$double.eps, trace = FALSE)
}
\arguments{
\item{x}{
  matrix of predictors
}
\item{y}{
  response 
}
\item{family}{
  name of a family function that represents the distribution of y to be
  used in the model. It must be \code{binomial}, \code{gaussian}, or
  \code{poisson}. For each one, the canonical link function is used;
  \code{logit} for binomial, \code{identity} for gaussian, and
  \code{log} for poisson distribution. Default is \code{binomial.}
}
\item{weight}{
  an optional vector of weights for observations
}
\item{max.steps}{
  an optional bound for the number of steps to be taken
}
\item{max.norm}{
  an optional bound for the L1 norm of the coefficients. Default is
  \code{100 * ncol(x).}
}
\item{min.lambda}{
  an optional (lower) bound for the size of \eqn{\lambda}. When
  \code{ncol(x)} is relatively large, the coefficient estimates are
  prone to numerical precision errors at extremely small
  \eqn{\lambda}. In such cases, early stopping is recommended (using
  \code{max.steps,} \code{max.norm,} or \code{min.lambda}). Default is
  \code{0} for \code{ncol(x) < nrow(x)} cases and \code{1e-6} otherwise.
}
\item{max.arclength}{
  an optional bound for arc length (L1 norm) of a step. If
  \code{max.arclength} is extremely small, an exact nonlinear path is
  produced. Default is \code{Inf.}
}
\item{standardize}{
  If \code{TRUE,} each column of x is standardized to have a unit
  variance before fitting.
}
\item{add.newvars}{
  \code{add.newvars} candidate variables (that are currently not in
  the active set) are used in the corrector step as potential active
  variables. Default is \code{1.}
}
\item{bshoot.threshold}{
  If the absolute value of a coefficient is larger than
  \code{bshoot.threshold} at the first corrector step it becomes nonzero
  (therefore when \eqn{\lambda} is considered to have been decreased too
  far), \eqn{\lambda} is increased again. i.e. A backward distance in
  \eqn{\lambda} that makes the coefficient zero is computed. Default is
  \code{0.5.}
}
\item{relax.lambda}{
  A variable joins the active set if \eqn{|l'(\beta)| >
    \lambda}(1-relax.lambda). Default is \code{1e-9.} If no variable
  joins the active set even after many (>20) steps, the user should
  increase \code{relax.lambda} to \code{1e-8} or \code{1e-7,} but not
  more than that. This adjustment is sometimes needed because of the
  numerical precision/error propagation problems. Overall, the paths are
  less accurate with relaxed lambda.     
}
\item{eps}{
  an effective zero
}
\item{trace}{
  If \code{TRUE,} the algorithm prints out its progress.
}
}
\value{
  A \code{glmpath} object is returned.
  \item{lambda}{
    vector of \eqn{\lambda} values for which exact coefficients are computed
  }
  \item{step.length}{
    vector of step lengths in \eqn{\lambda}
  }
  \item{new.df}{
    vector of degrees of freedom (to be used in the plot function)
  }
  \item{corr}{
    matrix of \eqn{l'(\beta)} values (derivatives of the log-likelihood)
  }
  \item{aic}{
    vector of AIC values 
  }
  \item{bic}{
    vector of BIC values 
  }
  \item{b.predictor}{
    matrix of coefficient estimates from the predictor steps
  }
  \item{b.corrector}{
    matrix of coefficient estimates from the corrector steps
  }
  \item{actions}{
    actions taken at each step
  }
  \item{meanx}{
    means of the columns of x
  }
  \item{sdx}{
    standard deviations of the columns of x
  }
  \item{xnames}{
    column names of x
  }
  \item{family}{
    family used
  }
  \item{weight}{
    weights used
  }
}
\details{
  This algorithm implements the predictor-corrector method to determine
  the entire path of the coefficient estimates as the amount of
  regularization varies; it computes a series of solution sets,
  each time estimating the coefficients with less regularization, based
  on the previous estimate. The coefficients are estimated with
  no error at the knots, and the values are connected, thereby making
  the paths piecewise linear.
}
\references{
  Mee Young Park and Trevor Hastie (2005) L1 Regularization Path
Algorithm for Generalized Linear Models - This will be soon available at the
authors' websites, \url{http://www.stanford.edu/~mypark} or
\url{http://stat.stanford.edu/~hastie/pub.htm}.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
cv.glmpath, plot.glmpath, predict.glmpath
}
\examples{
data(heart.data)
attach(heart.data)
fit.a <- glmpath(x, y, family=binomial)
fit.b <- glmpath(x, y, family=gaussian)
detach(heart.data)
}
\keyword{models}
\keyword{regression}