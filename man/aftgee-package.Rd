\name{aftgee-package}
\alias{aftgee-package}
%\alias{aftgee}
\docType{package}
\title{
  Accelerate Failure Time with Generalize Estimating Equations
}
\description{
  This package features both rank-based approach and least squares
  approach to fit accelerated failure time model.
  For the rank-based approach, the package allows various weight choices
  and uses an induced smoothing procedure to improve computational
  efficiency.
  With the rank-based estimator as an initial value, the generalized
  estimating equation approach is used as an extension of the least
  squares approach to the multivariate case.
  Additional sampling weights are incorporated to handle missing data
  needed as in case-cohort studies.
  }
\details{
  \tabular{ll}{
    Package: aftgee
    Type: Package
    Version: 0.4-1
    Date: 2014-01-26
    License: GPL(>=3)
    LazyLoad: yes
  }
}
\author{
  Sy Han Chiou \email{steven.chiou@uconn.edu}, Jun Yan
  \email{jun.yan@uconn.edu}, Sangwook Kang \email{sangwook.kang@uconn.edu}
}
\references{
  Chiou, S. H. and Kang, S. and Yan, J. (2013) Fast Accelerated Failure
  Time Modeling for Case-cohort Data. Statistics and Computing,
  pp. 1-10.
  doi:10.1007/s11222-013-9388-2.

  Jin, Z. and Lin, D. Y. and Ying, Z. (2006) On Least-squares Regression
  with Censored Data. Biometrika, 90, 341--353.

  Johnson, L. M. and Strawderman, R. L. (2009) Induced Smoothing for the
  Semiparametric Accelerated Failure Time Model: Asymptotic and
  Extensions to Clustered Data. Biometrika, 96, 577--590.
}
\keyword{ aftgee }
