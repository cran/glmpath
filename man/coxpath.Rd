\name{coxpath}
\alias{coxpath}
\title{
  Fits the entire L1 regularization path for Cox proportional hazards
  model
}
\description{
  This algorithm uses predictor-corrector method to compute the
  entire regularization path for Cox proportional hazards model with L1
  penalty.
}
\usage{
  coxpath(data, nopenalty.subset = NULL, method = c("breslow", "efron"),
          lambda2 = 1e-5, max.steps = 10 * min(n, m), max.norm = 100 * m,
          min.lambda = (if (m >= n) 1e-3 else 0), max.vars = Inf,
          max.arclength = Inf, frac.arclength = 1, add.newvars = 1,
          bshoot.threshold = 0.1, relax.lambda = 1e-7,
          approx.Gram = FALSE, standardize = TRUE,
          eps = .Machine$double.eps, trace = FALSE)
}
\arguments{
  \item{data}{
    a list consisting of \code{x:} a matrix of features, \code{time:}
    the survival time, and \code{status:} censor status with 1 if died
    and 0 if censored. 
  }
  \item{nopenalty.subset}{
    a set of indices for the predictors that are not subject to the L1
    penalty
  }
  \item{method}{
    approximation method for tied survival times. Approximations derived
    by Breslow (1974) and Efron (1977) are available. Default is
    \code{breslow.}
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
    an optional (lower) bound for the size of \eqn{\lambda}. When
    \code{ncol(x)} is relatively large, the coefficient estimates are
    prone to numerical precision errors at extremely small
    \eqn{\lambda}. In such cases, early stopping is recommended. Default
    is \code{0} for \code{ncol(x) < nrow(x)} cases and \code{1e-3}
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
      \lambda}*(1-\code{relax.lambda}). Default is \code{1e-7.} If no
    variable joins the active set even after many (>20) steps, the user
    should increase \code{relax.lambda} to \code{1e-6} or \code{1e-5,}
    but not more than that. This adjustment is sometimes needed because
    of the numerical precision/error propagation problems. In general,
    the paths are less accurate with relaxed lambda.     
  }
  \item{approx.Gram}{
    If \code{TRUE,} an approximated Gram matrix is used in predictor
    steps; each step takes less number of computations, but the total
    number of steps usually increases. This might be useful when the
    number of features is large.    
  }
  \item{standardize}{
    If \code{TRUE,} predictors are standardized to have a unit variance.
  }
  \item{eps}{
    an effective zero
  }
  \item{trace}{
    If \code{TRUE,} the algorithm prints out its progress.
  }
}
\value{
  A \code{coxpath} object is returned.
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
    matrix of \eqn{l'(\beta)} values (derivatives of the
    log-partial-likelihood)
  }
  \item{new.df}{
    vector of degrees of freedom (to be used in the plot function)
  }
  \item{df}{
    vector of degrees of freedom at each step
  }
  \item{loglik}{
    vector of log-partial-likelihood computed at each step
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
  \item{method}{
    method used
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
}
\references{
  Mee Young Park and Trevor Hastie (2007) L1 regularization path
  algorithm for generalized linear models. \emph{J. R. Statist. Soc.} B,
  69, 659-677.
}
\author{Mee Young Park and Trevor Hastie}
\seealso{
cv.coxpath, plot.coxpath, predict.coxpath, summary.coxpath
}
\examples{
data(lung.data)
attach(lung.data)
fit.a <- coxpath(lung.data)
fit.b <- coxpath(lung.data, method="efron")
detach(lung.data)
}
\keyword{models}
\keyword{regression}
