\name{aftgee}
\alias{aftgee}
\alias{summary.aftgee}
\title{
  Accelerate Failure Time with Generalize Estimating Equation
}
\description{
  A package that uses Generalized Estimating Equations (GEE) to estimate
  Multivariate Accelerated Failure Time Model (AFT).
  This package allows dependence working correlation structure for GEE.
  Moreover, it gives two options of initial estimator,
  the simple linear regression and the smooth-hegan weight estimator.
}
\usage{
aftgee(formula, data, subset, id,
       contrasts = NULL,
       weights = NULL,
       margin  = NULL,
       corstr = "independence",
       binit = "lm",
       B = 100,
       control = aftgee.control())
}
\arguments{
  \item{formula}{ A formula expression as for \code{glm}, of the form
    \code{response ~ predictors}. Response can be in class of survival data.
    See the documentation of lm, coxph and formula for details.
  }
  \item{id}{ A vector which identifies the clusters.
	The length of 'id' should be the same as the number of observations.
  }
  \item{data}{ An optional data frame in which to interpret the variables
    occurring in the \code{formula}, along with the \code{id}.
  }
  \item{subset}{ An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{corstr}{ a character string specifying the correlation structure.
    The following are permitted:
    \code{'independence'},
    \code{'exchangeable'},
    \code{'ar1'},
    \code{'unstructured'},
    \code{'userdefined'}, and
    \code{'fixed'}
  }
  \item{B}{ is the sample size for resampling estimation of the variance
  }
  \item{contrasts}{ an optional list.
  }
  \item{binit}{ a character string specifying the initial estimator.
    The following are permitted:
    \code{"lm"} for simple linear regression,
    \code{"srrgehan"} for smoothed gehan weight estimator.
  }
  \item{weights}{ is a weight vector for covariates; default at 1.
  }
  \item{margin}{ a sformula vector; default at 1.
  }
  \item{control}{ Controls maxiter and tolerance.
  }
}

\details{
  \tabular{ll}{
    Package:  aftgee
    Type:  Package
    Version:  0.2-19
    Date:  2012-03-12
    License: GPL (>=3)
    LazyLoad:  yes
  }
}

\value{
  An object of class "aftgee" representing the fit.
  An object of class "aftgee" is a list containing at elast the
  following components:
  \item{coefficients}{a vector of initial value and a vector of point estimates}
  \item{coef.res}{a vector of point estimates}
  \item{var.res}{estimated covariance matrix}
  \item{coef.init}{a vector of initial value}
  \item{var.init.mat}{estimated initial covariance matrix}
  \item{initial}{a character string specifying the initial estimator.}
  \item{conv}{An integer code indicating type of convergence after GEE
  iteration. 0 indicates successful convergence; 1 indicates that the
  iteration limit \code{maxit} has been reached }
  \item{ini.conv}{An integer code indicating type of convergence for
  initial value. 0 indicates successful convergence; 1 indicates that the
  iteration limit \code{maxit} has been reached }
  \item{conv.step}{An integer code indicating the step until convergence}
}

\author{
  Sy Han Chiou, Sangwook Kang, Jun Yan
}
\references{
  Chiou, S. H. and Kim, H. and Yan, J. (2012) Semiparametric
  Multivariate Accelerated Failure Time Model with Generalized
  Estimating Equations. University of Connecticut. %% report number?

  Jin, Z. and Lin, D. Y. and Ying, Z. (2006) On Least-squares Regression
  with Censored Data. Biometrika, 90, 341--353.

  Johnson, L. M. and Strawderman, R. L. (2009) Induced Smoothing for the
  Semiparametric Accelerated Failure Time Model: Asymptotics and
  Extensions to Clustered Data. Biometrika, 96, 577-590.
}

\keyword{ aftgee }

\examples{
library(survival)
data(kidney)
fit <- aftgee(Surv(time, status) ~ age + sex, id = id, data = kidney, binit = "srrgehan")
fit
summary(fit)
}
