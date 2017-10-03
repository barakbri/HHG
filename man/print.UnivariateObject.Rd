\name{print.UnivariateObject}
\alias{print.UnivariateObject}

\title{Print function for Univariate Null Table Object}

\description{Print description of univariate object. }

\usage{
\method{print}{UnivariateObject}(x, ...)
}
\arguments{
  \item{x}{A univariate null table object, of type 'UnivariateObject'}
  \item{...}{Additional arguments can be sent to function. Currently not supported.}
}

\details{
Function prints description of a null table object, including sample size (for the independence hhg.univariate statsitics)
or group sizes (for hhg.univarate statistics). Also prints statistic type in terms of \code{variant}, \code{aggregation.type} and \code{score.type}.
}

\value{
Does not return value. Only prints description of null table object.
}


\author{
  Barak Brill
}


\examples{
#univariate objects are found inside null tables:
nt = hhg.univariate.ks.nulltable(group.sizes = c(20,20))

print(nt$univariate.object)
}
