\name{diabetes}
\docType{data}
\alias{diabetes}
\alias{diabetes1}
\title{Blood and other measurements in diabetics}
\description{
  The \code{diabetes1} data frame has 442 rows and 1 columns.
  These are the data used in the Efron et al "Least Angle Regression" paper.

  The \code{diabetes} data frame has 442 rows and 3 matrices, containing
  predictors, response, and interactions.
}
\format{
  \code{diabetes1} is a data frame with 442 observations on the following 11 variables.
  \describe{
    \item{\code{age}}{a numeric vector}
    \item{\code{sex}}{a numeric vector}
    \item{\code{bmi}}{a numeric vector}
    \item{\code{map}}{a numeric vector}
    \item{\code{tc}}{a numeric vector}
    \item{\code{ldl}}{a numeric vector}
    \item{\code{hdl}}{a numeric vector}
    \item{\code{tch}}{a numeric vector}
    \item{\code{ltg}}{a numeric vector}
    \item{\code{glu}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
  In the \code{sex} variable, 1 indicates female and 2 male.

  \code{diabetes} is a data frame containing the following objects:
  \describe{
    \item{x}{a matrix with 10 columns--the first 10 columns from
      \code{diabetes1}, standardized}
    \item{y}{a numeric vector}
    \item{x2}{a matrix with 64 columns--main effects and second-order interactions}
  }
}
\details{
  The \code{x} matrix is standardized to have unit L2 norm in each column
  and zero mean. The matrix \code{x2} consists of \code{x} plus
  second-order powers and interactions, also standardized.
}
\source{
  \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.ps}
}
\references{
  B. Efron and T. Hastie (2003),
  "LARS software for R and Splus",
  \code{http://www-stat.stanford.edu/~hastie/Papers/LARS}

  B. Efron, T. Hastie, I. Johnstone and R. Tibshirani (2004),
  "Least Angle Regression" (with discussion),
  \emph{Annals of Statistics}
  \bold{32}, 407-499.
}
\examples{

data(diabetes)

stepFit <- biglars.fit(diabetes$x, diabetes$y, type = "stepwise")
stepFitBlocked <- biglars.fit(diabetes$x, diabetes$y, type = "stepwise",
                           blockSize = 50)

lassoFit <- biglars.fit(diabetes$x, diabetes$y)
lassoFitBlocked <- biglars.fit(diabetes$x, diabetes$y, blockSize = 34)
}
\keyword{datasets}

