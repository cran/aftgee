\name{aftgee.control}
\alias{aftgee.control}
\title{Auxiliary for Controlling AFTGEE Fitting}
\description{
  Auxiliary function as user interface for 'aftgee' and 'aftsrr'
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
  When 'trace' is true, output for each iteration is printed to the
  screen.
}
\value{
  A list with the arguments as components.
}
\author{Chiou, S. H. \email{steven.chiou@uconn.edu}, Kang,
  S. \email{sangwook.kang@uconn.edu}, Yan, J. \email{jun.yan@uconn.edu}}
\seealso{'aftgee', 'aftgee.fit' and 'aftgee.est'.}

