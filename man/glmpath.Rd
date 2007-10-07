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
  glmpath(x, y, data, nopenalty.subset = NULL, family = binomial,
          weight = rep(1, n), offset = rep(0, n), lambda2 = 1e-5,
          max.steps = 10*min(n, m), max.norm = 100*m,
          min.lambda = (if (m >= n) 1e-6 else 0), max.vars = Inf,
          max.arclength = Inf, frac.arclength = 1, add.newvars = 1,
          bshoot.threshold = 0.1, relax.lambda = 1e-8,
          standardize = TRUE, function.precision = 3e-13,
          eps = .Machine$double.eps, trace = FALSE)
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
    response. \code{data} is not needed if \code{x} and \code{y} are
    input separately.
  }
  \item{nopenalty.subset}{
    a set of indices for the predictors that are not subject to the L1
    penalty
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
    offset, the corresponding column must be removed from \code{x.}
  }
  \item{lambda2}{
    regularization parameter for the L2 norm of the
    coefficients. Default is \code{1e-5.}
  }
  \item{max.steps}{
    an optional bound for the number of steps to be taken. Default is
    \code{10 * min{nrow(x), ncol(x)}.}
  }
  \item{max.norm}{
    an optional bound for the L1 norm of the coefficients. Default is
    \code{100 * ncol(x).}
  }
  \item{min.lambda}{
    an optional (lower) bound for the size of \eqn{\lambda}. Default is
    \code{0} for \code{ncol(x) < nrow(x)} cases and \code{1e-6}
    otherwise.    
  }
  \item{max.vars}{
    an optional bound for the number of active variables. Default is
    \code{Inf.}
  }
  \item{max.arclength}{
    an optional bound for arc length (L1 norm) of a step. If
    \code{max.arclength} is extremely small, an exact nonlinear path is
    produced. Default is \code{Inf.}
  }
  \item{frac.arclength}{
    Under the default setting, the next step size is computed so that
    the active set changes right at the next value of lambda. When
    \code{frac.arclength} is assigned some fraction between 0 and 1, the
    step size is decreased by the factor of \code{frac.arclength} in arc
    length. If \code{frac.arclength=0.2,} the step length is adjusted so
    that the active set would change after five smaller steps. Either
    \code{max.arclength} or \code{frac.arclength} can be used to force
    the path to be more accurate. Default is \code{1.}
  }
  \item{add.newvars}{
    \code{add.newvars} candidate variables (that are currently not in
    the active set) are used in the corrector step as potential active
    variables. Default is \code{1.}
  }
  \item{bshoot.threshold}{
    If the absolute value of a coefficient is larger than
    \code{bshoot.threshold} at the first corrector step it becomes
    nonzero (therefore when \eqn{\lambda} is considered to have been
    decreased too far), \eqn{\lambda} is increased again. i.e. A
    backward distance in \eqn{\lambda} that makes the coefficient zero
    is computed. Default is \code{0.1.}
  }
  \item{relax.lambda}{
    A variable joins the active set if \eqn{|l'(\beta)| >
      \lambda}*(1-\code{relax.lambda}). Default is \code{1e-8.} If no
    variable joins the active set even after many (>20) steps, the user
    should increase \code{relax.lambda} to \code{1e-7} or \code{1e-6,}
    but not more than that. This adjustment is sometimes needed because
    of the numerical precision/error propagation problems. In general,
    the paths are less accurate with relaxed lambda.     
  }
  \item{standardize}{
    If \code{TRUE,} predictors are standardized to have a unit variance.
  }
  \item{function.precision}{
    \code{function.precision} parameter used in the internal
    solver. Default is \code{3e-13.} The algorithm is faster, but less
    accurate with relaxed, larger \code{function precision.}
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
    vector of \eqn{\lambda} values for which the exact coefficients are
    computed
  }
  \item{lambda2}{
    \eqn{\lambda_2} used
  }
  \item{step.length}{
    vector of step lengths in \eqn{\lambda}
  }
  \item{corr}{
    matrix of \eqn{l'(\beta)} values (derivatives of the log-likelihood)
  }
  \item{new.df}{
    vector of degrees of freedom (to be used in the plot function)
  }
  \item{df}{
    vector of degrees of freedom at each step
  }
  \item{deviance}{
    vector of deviance computed at each step
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
  \item{new.A}{
    vector of boolean values indicating the steps at which the active
    set changed (to be used in the plot/predict functions)
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
  \item{offset}{
    offset used
  }
  \item{nopenalty.subset}{
    nopenalty.subset used
  }
  \item{standardize}{
    \code{TRUE} if the predictors were standardized before fitting
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

  We thank Michael Saunders of SOL, Stanford University for providing
  the solver used for the convex optimization in corrector steps of
  glmpath.
}
\references{
  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
cv.glmpath, plot.glmpath, predict.glmpath, summary.glmpath
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
