\name{biglars.fit}
\alias{biglars.fit}
\alias{biglars.fit.lasso}
\alias{biglars.fit.lar}
\alias{biglars.fit.stepwise}
\title{
  Least-Angle, Lasso and Stepwise Regression
}
\description{
  Least-angle, lasso and stepwise methods for linear regression.
}
\usage{
biglars.fit(x, y, type = "lasso", removeColumns = TRUE,
    eps = sqrt(.Machine$double.eps), blockSize = NULL, maxStages = NULL)
}
\arguments{
  \item{x}{
    A matrix or \code{ff} of numeric predictors.
    The number of columns should not be larger than
    the number of rows. 
  }
  \item{y}{
    A numeric response vector or 1-column \code{ff}.
  }
  \item{type}{
    The type of regression to be performed. The usual choices are
    \describe{
      \item{lasso}{lasso method (L1 penalty)}
      \item{lar}{least-angle regression}
      \item{stepwise}{forward stepwise regression}
    }
    but see the details section below for other possibilities.
    Default is \code{"lasso"}.
  }
  \item{removeColumns}{
    A logical scalar indicating whether columns with small variance
    should be removed from consideration as predictors; default `TRUE'.
  }
  \item{eps}{
    Numerical tolerance used for assessment of sign, equality, rank
    determination, column removal, etc. The default is the
    square root of \code{.Machine\$double.eps},
    the relative machine precision.
  }
  \item{blockSize}{
    If \code{NULL}, the block size is determined using methods of class
     \code{ff}.  Otherwise the passed value is used.
  }
  \item{maxStages}{
     The maximum number of stages allowed in the algorithm.
     This argument applies only to the \code{lasso} option.
     The defaults is \code{2p} for \code{lasso}, where \code{p}
     if the number of predictors.
  }
}
\value{
  A list with the following elements:
  \item{coef}{An array of regression coefficients for each stage.}
  \item{moves}{Any array describing variables added or removed at each stage.}
  \item{RSS}{Residual sum of squares.}
}
\details{
   An intercept is always included in the regression.
   This function calls other routines to do the core calculations, one of
  \code{biglars.fit.lasso},
  \code{biglars.fit.lar}, or
  \code{biglars.fit.stepwise}.
  These functions are associated with Fraley et~al. (2007) and will not
  be undergoing further development except for things like bug fixes.
  For ongoing development of least-angle regression, see the \code{glars}
  library.
}
\references{
  B. Efron, T. Hastie, I. Johnstone and R. Tibshirani (2004),
  "Least Angle Regression" (with discussion),
  \emph{Annals of Statistics}
  \bold{32}, 407-499.

  C. Fraley and T. Hesterberg (2007),
  " Least-Angle Regression for Large Datasets",
   Technical Report, Insightful Corporation.
}
\seealso{
  \code{\link{qrBlockApply}}
}
\examples{

data(diabetes)

larFit <- biglars.fit(diabetes$x, diabetes$y, type = "lar")
larFitBlocked <- biglars.fit(diabetes$x, diabetes$y, type = "lar",
                           blockSize = 50)

lassoFit <- biglars.fit(diabetes$x, diabetes$y, type = "lasso")
lassoFitBlocked <- biglars.fit(diabetes$x, diabetes$y, type = "lasso",
                              blockSize = 34)

}
\keyword{regression}
