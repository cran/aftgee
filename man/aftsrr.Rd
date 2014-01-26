\name{aftsrr}
\alias{aftsrr}
\alias{summary.aftsrr}
\title{Accelerated Failure Time with Smooth Rank Regression}
\description{
  Fit semiparametric accelerated failure time (AFT) model with
  rank-based approach.
  Fast sandwich variance estimation, general weight functions and
  additional sampling weights are also incorporated.
}
\usage{
aftsrr(formula, data, id, subset, contrasts = NULL,
       strata = NULL, weights = NULL, rankWeights = "gehan",
       method = "sm", variance = "ISMB", B = 100, SigmaInit = NULL, 
       control = aftgee.control())
       }

\arguments{
  \item{formula}{ A formula expression, of the form \code{response ~ predictors}.
    Response can be in class of survival data.
    See the documentation of lm, coxph and formula for details.
  }
  \item{data}{ An optional data frame in which to interpret the variables
    occurring in the \code{formula}, along with the \code{id}.
  }
  \item{id}{ A vector which identifies the clusters.
    The length of 'id' should be the same as the number of observation.
  }
  \item{subset}{ An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{contrasts}{ an optional list.
  }
  \item{strata}{ A vector which identifies the strata. This can also
    be used to distinct case-cohort sampling and stratified sampling.
  }
  \item{weights}{ A weight vector for covariates; default at 1.
  }
  \item{rankWeights}{ A character string specifying the type of general
    weights.
    The following are permitted:
    \code{'logrank'}: logrank weight,
    \code{'gehan'}: Gehan's weight,
    \code{'PW'}: Prentice-Wilcoxon weight,
    \code{'GP'}: GP class weight.
  }
  \item{method}{ A character string specifying the methods to estimate
    the regression parameter.
    The following are permitted:
    \code{'nonsm'}: Regression parameters are estimated by directly solving the nonsmooth
estimating equations.
    \code{'sm'}: Regression parameters are estimated by directly solving the smoothed
estimating equations.
    \code{'monosm'}:Regression parameters are estimated by iterating the
  monotonic smoothed estimating equations. This is typical when
  \code{rankWeights = "PW"} and \code{rankWeights = "GP"}.
  }  
  \item{variance}{ A character string specifying the covariance estimating method.
    The following are permitted:
    \code{'MB'}: multiplier resampling,
    \code{'ZLCF'}: Zeng and Lin's approach with closed form Si,
    \code{'ZLMB'}: Zeng and Lin's approach with empirical Si,
    \code{'sHCF'}: Huang's approach with closed form Si,
    \code{'sHMB'}: Huang's approach with empirical Si,
    \code{'ISCF'}: Johnson and Strwderman's sandwich variance estimates with closed form Si,
    \code{'ISMB'}: Johnson and Strwderman's sandwich variance estimates with empirical Si,
    \code{'js'}: Johnson and Strwderman's iterating approach.
  }
  \item{B}{Resampling size. When M = 0, only the beta estimate will be displayed.
  }
  \item{SigmaInit}{ The initial covariance matrix; default at identity matrix.
  }
  \item{control}{ Controls maxiter and tolerance.
  }
}

\value{
  \code{aftsrr} returns an object of class "aftsrr" representing
  the fit.\\
  An object of class "aftsrr" is a list containing at least the
  following components:
  \item{beta}{A vector of beta estimates}
  \item{covmat}{A list of covariance estimates}
  \item{convergence}{An integer code indicating type of convergence. 0
    indicates successful convergence. Error codes are 1 indicates that
    the iteration limit \code{maxit} has been reached; 2 is failure
    due to stagnation; 3 indicates error in function evaluation; 4 is
    failure due to exceeding 100 step length reductions in line-search;
    and 5 indicates lack of improvement in objective function.}
  \item{bhist}{When \code{variance = "MB"}, \code{bhist} gives the
  bootstrap samples. 
  }
}

\references{
  Chiou, S. H. and Kang, S. and Yan, J. (2012)
  Fast Accelerated Failure Time Modeling for Case-Cohort Data.

  Zeng, D. and Lin, D. Y. (2008)
  Efficient Resampling Methods for Nonsmooth Estimating Functions.
  \emph{Biostatistics}, \bold{9}, 355--363

  Huang, Y. (2002)
  Calibration Regression of Censored Lifetime Medical Cost.
  \emph{Journal of American Statistical Association}, \bold{97}, 318--327

  Johnson, L. M. and Strawderman, R. L. (2009)
  Induced Smoothing for the Semiparametric Accelerated Failure Time Model: Asymptotic and Extensions to Clustered Data.
  \emph{Biometrika}, \bold{96}, 577 -- 590
}

\author{Steven Chiou \email{steven.chiou@uconn.edu}}
\examples{

#### kidney data
library(survival)
data(kidney)
foo <- aftsrr(Surv(time, status) ~ age + sex - 1, id = id,
                data = kidney, variance = c("MB", "ISMB"), B = 10)
foo

#### nwtco data
library(survival)
data(nwtco)
subinx <- sample(1:nrow(nwtco), 668, replace = FALSE)
nwtco$subcohort <- 0
nwtco$subcohort[subinx] <- 1
pn <- table(nwtco$subcohort)[[2]] / sum(table(nwtco$subcohort))
nwtco$hi <- nwtco$rel + ( 1 - nwtco$rel) * nwtco$subcohort / pn
nwtco$age12 <- nwtco$age / 12
nwtco$edrel = nwtco$edrel / 12
nwtco$study = nwtco$study - 3
nwtco$stage2 = ifelse(nwtco$stage == 2, 1, 0)
nwtco$stage3 = ifelse(nwtco$stage == 3, 1, 0)
nwtco$stage4 = ifelse(nwtco$stage == 4, 1, 0)
nwtco$histol = nwtco$histol - 1
sub <- nwtco[subinx,]
fit <- aftsrr(Surv(edrel, rel) ~ histol + age12 + study - 1, id = seqno,
       weights = hi, data = sub, B = 7, variance = c("MB", "ISMB"),
       subset = stage4 == 1)
summary(fit)
}

\keyword{aftsrr}

