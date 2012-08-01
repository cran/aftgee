\name{prepX}
\alias{prepX}
\title{
  Prepare Block Matrix for AFTGEE with multiple margins
}
\description{
  Prepare block matrix to be fitted with aftgee()
}
\usage{
prepX(xmat, margin, name = NULL)
}
\arguments{
  \item{xmat}{ a matrix to be prepared into a block diagonal matrix.
  }
  \item{margin}{ a vector specifying the error margin within cluster.
    The length of margin should be the same as the number of observations.
  }
  \item{name}{a optional character string specifying the column
    names of \code{xmat}.
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
\author{
  Sy Han Chiou, Sankwook Kang, Jun Yan
}
\examples{
data(kidney)
x = cbind(kidney$sex, kidney$age)
## head(x)
marg = rep(1:2, 38)
xmat = prepX(x, marg, c("sex", "age"))
## head(xmat)
}
