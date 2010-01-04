\name{qrBlockApply}
\alias{qrBlockApply}
\alias{qrBlock}
\title{
  Blockwise Cholesky Factorization
}
\description{
  Cholesky factorization of a crossproduct matrix via blockwise 
  orthogonal transformation.
}
\usage{
qrBlockApply(x, y = NULL, blockSize = NULL)
}
\arguments{
  \item{x}{
    A numeric matrix or \code{ff}.
    The number of columns should not be larger than
    the number of rows. 
  }
  \item{y}{
    An optional numeric vector or 1 column \code{ff} of responses.
  }
  \item{blockSize}{
     This is used for setting the block size for sequential access. 
     The default is determined by methods from class \code{ff}.
  }
}
\value{
  A list with the following elements:
  \item{R}{Upper triangular Cholesky factor of \code{crossprod(x)}.}
  \item{Qty}{The corresponding transformation of \code{y}.}
}
\references{
  C. Fraley and T. Hesterberg (2007),
  "Least-Angle Regression for Large Datasets",
   Technical Report, Insightful Corporation.
}
\seealso{
  \code{\link{biglars.fit}}
}
\examples{

data(diabetes)

qrx <- qr(diabetes$x)

list(R = qr.R(qrx), Qty = crossprod(qr.Q(qrx), diabetes$y))

qrBlockApply( diabetes$x, diabetes$y, blockSize = 34)
}
\keyword{regression}

