\name{aftgee-package}
\alias{aftgee-package}
%\alias{aftgee}
\docType{package}
\title{
  Accelerate Failure Time with Generalize Estimating Equations
}
\description{
  A package that uses Generalized Estimating Equations (GEE) to estimate
  Multivariate Accelerated Failure Time Model (AFT).
  This package implements recently developed inference procedures for
  AFT models with both the rank-based approach and the least squares
  approach.
  For the rank-based approach, the package allows various weight choices
  and uses an induced smoothing procedure that leads to much more
  efficient computation than the linear programming method.
  With the rank-based estimator as an initial value, the generalized
  estimating equation approach is used as an extension of the least
  squares approach to the multivariate case.
  Additional sampling weights are incorporated to handle missing data
  needed as in case-cohort studies or general sampling schemes.
}
\details{
  \tabular{ll}{
    Package: aftgee
    Type: Package
    Version: 1.0-0
    Date: 2014-11-12
    License: GPL(>=3)
    LazyLoad: yes
  }
}
\author{
  Sy Han Chiou \email{schiou@d.umn.edu},
  Jun Yan \email{jun.yan@uconn.edu},
  Sangwook Kang \email{kanggi1@yonsei.ac.kr}
}
\references{
  Chiou, S., Kim, J. and Yan, J. (2014) Marginal Semiparametric
  Multivariate Accelerated Failure Time Model with Generalized
  Estimating Equation. \emph{Life Time Data}, \bold{20}(4): 599--618.
  
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

  Jin, Z. and Lin, D. Y. and Ying, Z. (2006) On Least-squares Regression
  with Censored Data. \emph{Biometrika}, \bold{90}, 341--353.

  Johnson, L. M. and Strawderman, R. L. (2009)
  Induced Smoothing for the Semiparametric Accelerated Failure Time
  Model: Asymptotic and Extensions to Clustered Data.
  \emph{Biometrika}, \bold{96}, 577 -- 590

  Zeng, D. and Lin, D. Y. (2008)
  Efficient Resampling Methods for Nonsmooth Estimating Functions.
  \emph{Biostatistics}, \bold{9}, 355--363
}
\keyword{ aftgee }
