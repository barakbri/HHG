\name{print.HHG.Test.Result}
\alias{print.HHG.Test.Result}

\title{Print function for result of HHG tests}

\description{Print description of for result object of HHG tests }

\usage{
\method{print}{HHG.Test.Result}(x, ...)
}
\arguments{
  \item{x}{result of \code{hhg.test}, \code{hhg.test.2.sample} or 
  
  \code{hhg.test.k.sample}}
  \item{...}{Additional arguments can be sent to function. Currently not supported.}
}

\details{
Function prints description of results for the \code{hhg.test}, \code{hhg.test.2.sample} and
\code{hhg.test.k.sample} functions. Displays: test statistics, pvalues (if permutations were performed) and description
of sample size (also displays group sizes and equality of distribution tests).
}

\value{
Does not return value. Only prints description of test statistic and results.
}


\author{
  Barak Brill
}


\examples{
#output for independence test
n = 50
X = hhg.example.datagen(n, '4indclouds') 

Dx = as.matrix(dist((X[1,]), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((X[2,]), diag = TRUE, upper = TRUE))

hhg = hhg.test(Dx, Dy, nr.perm = 200)

#output for k-sample test
n = 50
D = hhg.example.datagen(n, 'FourClassUniv')
Dx = as.matrix(dist(D$x, diag = TRUE, upper = TRUE))

hhg = hhg.test.k.sample(Dx, D$y, nr.perm = 200)
}
