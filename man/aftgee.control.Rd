\name{aftgee.control}
\alias{aftgee.control}
\title{Auxiliary for Controlling AFTGEE Fitting}
\description{
  Auxiliary function as user interface for \code{aftgee} and \code{aftsrr}
  fitting. 
}
\usage{
aftgee.control(maxiter = 50,
               reltol = 0.001,
               trace = FALSE)
}
\arguments{
  \item{maxiter}{max number of iteration.}
  \item{reltol}{relative error tolerance.}
  \item{trace}{binary variable, determine whether to display output for
    each iteration.}
}
\details{
  When \code{trace} is true, output for each iteration is printed to the
  screen.
}
\value{
  A list with the arguments as components.
}
\author{
    Sy Han Chiou, Sangwook Kang, Jun Yan.
  }
\seealso{\code{aftgee}, \code{aftgee.fit} and \code{aftgee.est}.}

