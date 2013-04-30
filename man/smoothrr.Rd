\name{smoothrr}
\alias{smoothrr}
\alias{summary.smoothrr}
\title{Smooth Rank Regression}
\description{
  The smoothed gehan weight is used for the initial value. As proved in Johnson and Strawderman (2009), the smooth gehan estimate is consistent.
}
\usage{
smoothrr(formula, data, subset, contrasts = NULL, id,
         weights = NULL,
	 rankweights = "gehan", binit = "lm", sigmainit = NULL,
         variance = "ISMB", B = 100, gpweight = 1, strataid = NULL,
         iter = TRUE, control = aftgee.control())
}
\arguments{
  \item{formula}{ A formula expression, of the form \code{response ~ predictors}.
    Response can be in class of survival data.
    See the documentation of lm, coxph and formula for details.
  }
  \item{id}{ A vector which identifies the clusters.
    The length of 'id' should be the same as the number of observation.
  }
  \item{data}{ An optional data frame in which to interpret the variables
    occurring in the \code{formula}, along with the \code{id}.
  }
  \item{subset}{ An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{contrasts}{ an optional list.
  }
  \item{strataid}{ A vector which identifies the strata. This can also
    be used to distinct case-cohort sampling and stratified sampling.
  }
  \item{binit}{ A vector of initial estimator; default at simple linear regression.
  }
  \item{sigmainit}{ The initial covariance matrix; default at identity matrix.
  }
  \item{weights}{ A weight vector for covariates; default at 1.
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
  \item{rankweights}{ A character string specifying between "gehan" and "logrank" weight.}
  \item{B}{Resampling size. When M = 0, only the beta estimate will be displayed.
  }
  \item{gpweight}{ Power for general weight function when rankweights = GP
  }
  \item{iter}{ A logical value specifying whether iterative procedure is used for general weights  function
  }
  \item{control}{Controls maxiter, abstol and restol.
  }
}

\value{
  \code{smoothrr} returns an object of class "smoothrr" representing
  the fit.\\
  An object of class "smoothrr" is a list containing at least the
  following components:
  \item{beta}{a vector of beta estimates}
  \item{covmat}{a list of covariance estimates}
  \item{tbeta}{Timing for beta estimates}
  \item{beta.conv}{An integer code indicating type of convergence. 0
    indicates successful convergence. Error codes are 1 indicates that
    the iteration limit \code{maxit} has been reached; 2 is failure
    due to stagnation; 3 indicates error in function evaluation; 4 is
    failure due to exceeding 100 step length reductions in line-search;
    and 5 indicates lack of improvement in objective function.}
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
foo <- smoothrr(Surv(time, status) ~ age + sex - 1, id = id, 
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
fit <- smoothrr(Surv(edrel, rel) ~ histol + age12 + study - 1, id = seqno, 
       weights = hi, data = sub, B = 7, variance = c("MB", "ISMB"), 
       subset = stage4 == 1)
summary(fit)
}

\keyword{smoothrank}

