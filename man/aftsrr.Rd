\name{aftsrr}
\alias{aftsrr}
\alias{summary.aftsrr}
\title{Accelerated Failure Time with Smooth Rank Regression}
\description{
  Fits a semiparametric accelerated failure time (AFT) model with
  rank-based approach.
  General weights, additional sampling weights and fast sandwich
  variance estimations are also incorporated. 
  Estimating equations are solved with Barzilar-Borwein spectral method
  implemented as \code{BBsolve} in package \pkg{BB}. 
}
\usage{
aftsrr(formula, data, subset, id = NULL, contrasts = NULL,
       strata = NULL, weights = NULL, rankWeights = "gehan",
       method = "sm", variance = "ISMB", B = 100, SigmaInit = NULL, 
       control = aftgee.control())
       }

\arguments{
  \item{formula}{ a formula expression, of the form \code{response ~ predictors}.
    The \code{response} is a \code{Surv} object object with right censoring.
    See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
  }
  \item{data}{ an optional data frame in which to interpret the variables
    occurring in the \code{formula}.
  }
  \item{id}{ an optional vector used to identify the clusters.
    If missing, then each individual row of \code{data} is presumed to
    represent a distinct subject.
    The length of \code{id} should be the same as the number of observation.
  }
  \item{subset}{ an optional vector specifying a subset of observations
    to be used in the fitting process.
  }
  \item{contrasts}{ an optional list.
  }
  \item{strata}{ a vector which identifies strata.
    This can also be used to distinct case-cohort sampling and stratified sampling.
  }
  \item{weights}{ an optional vector of observation weights.
  }
  \item{rankWeights}{ a character string specifying the type of general
    weights.
    The following are permitted:
    \code{logrank}: logrank weight,
    \code{gehan}: Gehan's weight,
    \code{PW}: Prentice-Wilcoxon weight,
    \code{GP}: GP class weight.
  }
  \item{method}{ a character string specifying the methods to estimate
    the regression parameter.
    The following are permitted:
    \code{nonsm}: Regression parameters are estimated by directly solving the nonsmooth
estimating equations.
    \code{sm}: Regression parameters are estimated by directly solving the smoothed
estimating equations.
    \code{monosm}:Regression parameters are estimated by iterating the
  monotonic smoothed estimating equations. This is typical when
  \code{rankWeights = "PW"} and \code{rankWeights = "GP"}.
  }  
  \item{variance}{ a character string specifying the covariance estimating method.
    The following are permitted:
    \code{MB}: multiplier resampling,
    \code{ZLCF}: Zeng and Lin's approach with closed form Si,
    \code{ZLMB}: Zeng and Lin's approach with empirical Si,
    \code{sHCF}: Huang's approach with closed form Si,
    \code{sHMB}: Huang's approach with empirical Si,
    \code{ISCF}: Johnson and Strwderman's sandwich variance estimates with closed form Si,
    \code{ISMB}: Johnson and Strwderman's sandwich variance estimates with empirical Si,
    \code{js}: Johnson and Strwderman's iterating approach.
  }
  \item{B}{ a numeric value specifies the resampling number.
    When B = 0, only the beta estimate will be displayed.
  }
  \item{SigmaInit}{ the initial covariance matrix; default at identity matrix.
  }
  \item{control}{ controls maxiter and tolerance.
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
  Chiou, S., Kang, S. and Yan, J. (2014) Fast Accelerated Failure Time
  Modeling for Case-Cohort Data. \emph{Statistics and Computing},
  \bold{24}(4): 559--568

  Chiou, S., Kang, S. and Yan, J. (2014) Fitting Accelerated Failure
  Time Model in Routine Survival Analysis with {R} Package
  \pkg{Aftgee}. 
  \emph{Journal of Statistical Software}, \bold{61}(11): 1--23
 
  Huang, Y. (2002)
  Calibration Regression of Censored Lifetime Medical Cost.
  \emph{Journal of American Statistical Association}, \bold{97}, 318--327

  Johnson, L. M. and Strawderman, R. L. (2009)
  Induced Smoothing for the Semiparametric Accelerated Failure Time
  Model: Asymptotic and Extensions to Clustered Data.
  \emph{Biometrika}, \bold{96}, 577 -- 590

  Varadhan, R. and Gilbert, P. (2009)
  BB: An R Package for Solving a Large System of Nonlinear Equations and
  for Optimizing a High-Dimensional Nonlinear Objective Function.
  \emph{Journal of Statistical Software}, \bold{32}(4): 1--26
  
  Zeng, D. and Lin, D. Y. (2008)
  Efficient Resampling Methods for Nonsmooth Estimating Functions.
  \emph{Biostatistics}, \bold{9}, 355--363
}

\author{
  Sy Han Chiou, Sangwook Kang, Jun Yan.
}
\examples{

#### kidney data
library(survival)
data(kidney)
foo <- aftsrr(Surv(time, status) ~ age + sex - 1, id = id,
                data = kidney, variance = c("ISMB", "ZLMB"), B = 10)
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
nwtco$study <- nwtco$study - 3
nwtco$histol = nwtco$histol - 1
sub <- nwtco[subinx,]
fit <- aftsrr(Surv(edrel, rel) ~ histol + age12 + study, id = seqno,
       weights = hi, data = sub, B = 10, variance = c("ISMB", "ZLMB"),
       subset = stage == 4)
summary(fit)
}

\keyword{aftsrr}

