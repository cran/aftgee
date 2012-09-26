\name{smoothrr}
\alias{smoothrr}
\alias{summary.smoothrr}
\title{Smooth Rank Regression}
\description{
  The smoothed gehan weight is used for the initial value. As proved in Johnson and Strawderman (2009), the smoooth gehan estimate is consistent.
}
\usage{
smoothrr(formula, data, subset, id,
         Sigma = NULL,
         weight = NULL,
         binit = NULL,
         nres = 100, p0 = 1, p1 = 1,
         contrasts = NULL,
         variance = "ISMB",
	 rankweight = "gehan",
         control = aftgee.control())
}
\arguments{
  \item{formula}{ Aa formula expression, of the form \code{respone ~ predictors}. 
    Response can be in class of survival data. 
    See thye documentation of lm, coxph and formula for details.
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
  \item{binit}{ A vector of initial estimator; default at linear model.
  } 
  \item{Sigma}{ The initial covariance matrix; default at identity matrix.
  }
  \item{weight}{ A weight vector for covariates; default at 1.
  }
  \item{variance}{ A character string specifying the covariance estimating method.
    The following are permitted:
    \code{'MB'}: multipiler resmapling,
    \code{'ZLCF'}: Zeng and Lin's approach with closed form Si,
    \code{'ZLMB'}: Zeng and Lin's approach with empirical Si,
    \code{'sHCF'}: Huang's appraoch with closed form Si,
    \code{'sHMB'}: Huang's approach with empirical Si,
    \code{'ISCF'}: Johnson and Strwderman's sandwish variance estimates with closed form Si,
    \code{'ISMB'}: Johnson and Strwderman's sandwish variance estimates with empirical Si,
    \code{'js'}: Johnson and Strwderman's iterating approach.
  }
  \item{rankweight}{ A character string specifying between "gehan" and "logrank" weight.}
  \item{nres}{Resampling size. When M = 0, only the beta estimate will be displayed.
  }
  \item{p0}{Case-cohort size ratio.
  }
  \item{p1}{Case-cohort size ratio.
  }
  \item{control}{Controls maxiter, abstol and restol.
  }
}

\value{
  \code{smoothrr} returns an object of class "smoothrr" representing
  the fit.\\
  An object of class "smoothrr" is a list containing at elast the
  following components:
  \item{beta}{a vector of beta estimates}
  \item{covmat}{a list of covariance estimates}
  \item{tbeta}{Timing for beta estimates}
  \item{beta.conv}{An integer code indicating type of convergence. 0
    indicates successful convergence. Error codes are 1 indicates that
    the iteration limit \code{maxit} has been reached; 2 is failure
    due to stagnation; 3 indicates error in function evaluation; 4 is
    failure due to exceeding 100 steplength reductions in line-search;
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
  Induced Smoothing for the Semiparametric Accelerated Failure Time Model: Asymptotics and Extensions to Clustered Data.
  \emph{Biometrika}, \bold{96}, 577 -- 590
}

\author{Steven Chiou \email{steven.chiou@uconn.edu}}
\examples{

#### kidney data
library(survival)
data(kidney)
foo <- smoothrr(Surv(log(time), status) ~ age + sex - 1, id = id, data = kidney, variance = c("MB", "ISMB"), nres = 10)
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
fit <- smoothrr(Surv(edrel, rel) ~ histol + age12 + study - 1, id = seqno, weight = hi, data = sub, nres = 7, variance = c("MB", "ISMB"), subset = stage4 == 1)
summary(fit)
} 

\keyword{smoothrank}

