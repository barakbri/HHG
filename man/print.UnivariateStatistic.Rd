\name{print.UnivariateStatistic}
\alias{print.UnivariateStatistic}

\title{Print function for Univariate Statistic Test Object}

\description{Print description of univariate statistic result & test object. }

\usage{
\method{print}{UnivariateStatistic}(x, ...)
}
\arguments{
  \item{x}{result of \code{hhg.univariate.ind.stat}, \code{hhg.univariate.ks.stat}, 
  
  \code{hhg.univariate.ind.combined.test}  or \code{hhg.univariate.ks.combined.test}}
  \item{...}{Additional arguments can be sent to function. Currently not supported.}
}

\details{
Function prints description of univariate test statistic and test results. Displays: test statistics, pvalues (for the combined test function), partition sizes and type of test statistic used in terms of \code{variant}, \code{aggregation.type} and \code{score.type}.
}

\value{
Does not return value. Only prints description of test statistic and results.
}


\author{
  Barak Brill
}


\examples{

#generate statistics and test results, and print them
ind.stat = hhg.univariate.ind.stat(1:20,1:20,variant = 'ADP',
  aggregation.type = 'sum',score.type = 'both',mmax = 5)
print(ind.stat)

ks.stat = hhg.univariate.ks.stat(1:50,sample(c(rep(0,25),rep(1,25))),
  aggregation.type = 'both',score.type = 'both',mmax = 10)
print(ks.stat)

\dontrun{

ind.combined = hhg.univariate.ind.combined.test(1:20,1:20,
  combining.type = 'Both',mmax = 5,nr.perm = 100)
print(ind.combined)

}

ks.combined = hhg.univariate.ks.combined.test(1:50,
  sample(c(rep(0,25),rep(1,25))),combining.type = 'Both')
print(ks.combined)

}
