\name{hhg.univariate.ks.stat}
\alias{hhg.univariate.ks.stat}


\title{The K-sample test statistics for all partition sizes}

\description{These statistics are used in the omnibus distribution-free test of equality of distributions among K groups, as described in Heller et al. (2014). }

\usage{
hhg.univariate.ks.stat(x, y,variant = 'KSample-Variant',aggregation.type='sum',
score.type='LikelihoodRatio', mmax = max(4,round(min(table(y))/3)),mmin=2,
nr.atoms= nr_bins_equipartition(length(x)))
}

\arguments{
  \item{x}{a numeric vector of data values. Tied observations are broken at random.}
  \item{y}{for \code{k} groups, a vector of integers with values \code{0:(k-1)} which specify the group each observation belongs to.}
  \item{variant}{Default value is \code{'KSample-Variant'}. Setting the variant to \code{'KSample-Equipartition'} performs the K-sample tests over partitions of the data where splits between cells are at least \eqn{n/nr.atoms} apart.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default), \code{"max"}, or \code{"both".}} 
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is 1/3 the number of observations in the smallest group.}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
  \item{nr.atoms}{For \code{variant=='KSample-Equipartition'} type tests, sets the number of possible split points in the data. The default value is the minimum between \eqn{n} and \eqn{60+0.5*\sqrt{n}}.}
}

\details{
  For each partition size \eqn{m= mmin,\ldots,mmax}, the function computes the scores in each of the partitions (according to score type), and aggregates all scores according to the aggregation type (see details in Heller et al. , 2014). If the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}, then the computed statistic will be in \code{statistic}, otherwise the computed statistics will be in the appropriate subset of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}.
  
  For the 'sum' aggregation type (default), The test statistic is the sum of log likelihood (or Pearson Chi-square) scores, of all partitions of size \eqn{m} of the data, normalized by the number of partitions and the data size (thus, being an estimator of the Mutual Information). For the 'max' aggregation type, the test statistic is the maximum log likelihood (or Pearson Chi-square) score acheived by a partition of data of size \code{m}.
  
  Variant type \code{"KSample-Equipartition"} is the computationally efficient version of the K-sample test. calculation time is reducing by aggregating over a subset of partitions, where a split between cells may be performed only every \eqn{n/nr.atoms} observations. This allows for a complexity of O(nr.atoms^2) (instead of O(n^2)). Computationly efficient versions are available for \code{aggregation.type=='sum'} and \code{aggregation.type=='max'} variants.
}

\value{
  Returns a \code{UnivariateStatistic} class object, with the following entries:
  
  \item{statistic}{The value of the computed statistic if the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}. One of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}. }
  
  \item{sum.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all Pearson chi-squared statistics from all the \eqn{K X m} contingency tables considered, divided by the total number of observations.} 
  
  \item{sum.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all LikelihoodRatio statistics from all the \eqn{K X m} contingency tables considered, divided by the total number of observations.} 
  
 \item{max.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{K X m} contingency tables considered.} 
  
  \item{max.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{K X m} contingency tables considered.} 
  
  \item{type}{"KSample".}
  
  \item{stat.type}{"KSample".}
  
  \item{size}{A vector of size K of the ordered group sample sizes.}
  
  \item{score.type}{The input \code{score.type}.}
  
  \item{aggregation.type}{The input \code{aggregation.type}.}
  
  \item{mmin}{The input \code{mmin}.}
  
  \item{mmax}{The input \code{mmax}.}
  
  \item{nr.atoms}{The input \code{nr.atoms}.}

}

\references{

Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, JMLR 17(29):1-54

Brill B. (2016) Scalable Non-Parametric Tests of Independence (master's thesis)

http://primage.tau.ac.il/libraries/theses/exeng/free/2899741.pdf

}

\author{
  Barak Brill and Shachar Kaufman.
}



\examples{
#Example of computing the test statisics for data from a two-sample problem:

#Two groups, each from a different normal mixture:
X = c(c(rnorm(25,-2,0.7),rnorm(25,2,0.7)),c(rnorm(25,-1.5,0.5),rnorm(25,1.5,0.5)))
Y = (c(rep(0,50),rep(1,50)))
plot(Y,X)


#I) Computing test statistics , with default parameters:
hhg.univariate.Sm.Likelihood.result = hhg.univariate.ks.stat(X,Y)

hhg.univariate.Sm.Likelihood.result

#II) Computing test statistics , with max aggregation type:
hhg.univariate.Mm.likelihood.result = hhg.univariate.ks.stat(X,Y,aggregation.type = 'max')

hhg.univariate.Mm.likelihood.result


#III) Computing statistics, which are computationaly efficient for large data:

#Two groups, each from a different normal mixture, total sample size is 10^4:
X_Large = c(c(rnorm(2500,-2,0.7),rnorm(2500,2,0.7)),
c(rnorm(2500,-1.5,0.5),rnorm(2500,1.5,0.5)))
Y_Large = (c(rep(0,5000),rep(1,5000)))
plot(Y_Large,X_Large)

# for these variants, make sure to change mmax so that mmax<= nr.atoms

hhg.univariate.Sm.EQP.Likelihood.result = hhg.univariate.ks.stat(X_Large,Y_Large,
variant = 'KSample-Equipartition',mmax=30)

hhg.univariate.Sm.EQP.Likelihood.result

hhg.univariate.Mm.EQP.likelihood.result = hhg.univariate.ks.stat(X_Large,Y_Large,
aggregation.type = 'max',variant = 'KSample-Equipartition',mmax=30)

hhg.univariate.Mm.EQP.likelihood.result

}
