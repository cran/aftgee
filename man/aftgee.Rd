\name{aftgee}
\alias{aftgee}
\alias{summary.aftgee}
\title{
  Accelerated Failure Time with Generalize Estimating Equation
}
\description{
  Fits a semiparametric accelerated failure time (AFT) model with
  least-squares approach.
  Generalized estimating equation is generalized to multivariate AFT
  modeling to account for multivariate dependence through working
  correlation structures to improve efficiency. 
}
\usage{
aftgee(formula, data, subset, id = NULL,
       contrasts = NULL, weights = NULL,
       margin = NULL, corstr = "independence",
       binit = "srrgehan", B = 100,
       control = aftgee.control())
}
\arguments{
  \item{formula}{ a formula expression, of the form \code{response ~ predictors}.
    The \code{response} is a \code{Surv} object object with right censoring.
    In the case of no censoring, \code{aftgee} will return an ordinary
    least estimate when \code{corstr = "independence"}. 
    See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
  }
  \item{data}{ an optional data frame in which to interpret the variables
    occurring in the \code{formula}.
  }
  \item{id}{ an optional vector used to identify the clusters.
    If missing, then each individual row of \code{data} is presumed to
    represent a distinct subject.
    The length of \code{id} should be the same as the number of observations.
  }
  \item{subset}{ an optional vector specifying a subset of observations
    to be used in the fitting process.
  }
  \item{corstr}{ a character string specifying the correlation structure.
    The following are permitted:
    \code{independence},
    \code{exchangeable},
    \code{ar1},
    \code{unstructured},
    \code{userdefined}, and
    \code{fixed}
  }
  \item{B}{ a numeric value specifies the resampling number.
    When B = 0, only the beta estimate will be displayed.
  }
  \item{contrasts}{ an optional list.
  }
  \item{binit}{ can be either a vector or a character string specifying
    the initial slope estimator.
    When \code{binit} is a vector, its length should be the same as the
       number of covariates. 
    When \code{binit} is a character string, the following are permitted:
    \code{lm} for simple linear regression,
    \code{srrgehan} for smoothed gehan weight estimator.
  }
  \item{weights}{ an optional vector of observation weights.
  }
  \item{margin}{ a sformula vector; default at 1.
  }
  \item{control}{ controls maxiter and tolerance.
  }
}

\details{
  \tabular{ll}{
    Package:  aftgee
    Type:  Package
    Version:  0.4-3
    Date:  2014-04-08
    License: GPL (>=3)
    LazyLoad:  yes
  }
}

\value{
  An object of class "aftgee" representing the fit.
  An object of class "aftgee" is a list containing at least the
  following components:
  \item{coefficients}{a vector of initial value and a vector of point estimates}
  \item{coef.res}{a vector of point estimates}
  \item{var.res}{estimated covariance matrix}
  \item{coef.init}{a vector of initial value}
  \item{var.init.mat}{estimated initial covariance matrix}
  \item{binit}{a character string specifying the initial estimator.}
  \item{conv}{An integer code indicating type of convergence after GEE
  iteration. 0 indicates successful convergence; 1 indicates that the
  iteration limit \code{maxit} has been reached }
  \item{ini.conv}{An integer code indicating type of convergence for
  initial value. 0 indicates successful convergence; 1 indicates that the
  iteration limit \code{maxit} has been reached }
  \item{conv.step}{An integer code indicating the step until convergence}
}

\author{
  Sy Han Chiou, Sangwook Kang, Jun Yan.
}
\references{
  Chiou, S., Kim, J. and Yan, J. (2014) Marginal Semiparametric
  Multivariate Accelerated Failure Time Model with Generalized
  Estimating Equation. \emph{Life Time Data}, \bold{20}(4): 599--618.
  
  Jin, Z. and Lin, D. Y. and Ying, Z. (2006) On Least-squares Regression
  with Censored Data. \emph{Biometrika}, \bold{90}, 341--353.
}

\keyword{ aftgee }

\examples{
library(survival)
library(copula)
datgen <- function(n = 100, tau = 0.3, cen = 75.4, dim = 2) {
    kt <- iTau(claytonCopula(1), tau)
    copula <- claytonCopula(kt, dim = dim)
    id <- rep(1:n, rep(dim, n))
    x1 <- rbinom(dim * n, 1, 0.5)
    x2 <- rnorm(dim * n)
    ed <- mvdc(copula, rep("weibull", dim), rep(list(list(shape = 1)), dim))
    e <- c(t(rMvdc(n, ed)))
    T <- exp(2 + x1 + x2 + e)
    cstime <- runif(n, 0, cen)
    delta <- (T < cstime) * 1
    Y <- pmin(T, cstime)
    out <- data.frame(T = T, Y = Y, delta = delta, x1 = x1, x2 = x2, id = rep(1:dim, n))
    out
}
set.seed(1)
mydata <- datgen(n = 50, dim = 2)
summary(aftgee(Surv(Y, delta) ~ x1 + x2, data = mydata, 
               id = id, corstr = "ind", B = 10))
summary(aftgee(Surv(Y, delta) ~ x1 + x2, data = mydata, 
               id = id, corstr = "ex", B = 10))
}
