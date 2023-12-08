\name{hhg.univariate.ind.stat}
\alias{hhg.univariate.ind.stat}


\title{The independence test statistics for all partition sizes}

\description{These statistics are used in the omnibus distribution-free test of independence between two univariate random variables, as described in Heller et al. (2016).}

\usage{
hhg.univariate.ind.stat(x, y, variant = 'ADP',aggregation.type='sum',
score.type='LikelihoodRatio', mmax = max(floor(sqrt(length(x))/2),2),
mmin =2, w.sum = 0, w.max = 2,nr.atoms = nr_bins_equipartition(length(x)))
}

\arguments{
  \item{x}{a numeric vector with observed \code{X} values (tied observations are broken at random).}
  \item{y}{a numeric vector with observed \code{Y} values (tied observations are broken at random).}
  \item{variant}{a character string specifying the partition type, must be one of \code{"ADP"} (default) or \code{"DDP"}, \code{"ADP-ML"}, \code{"ADP-EQP"}, \code{"ADP-EQP-ML"}.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default), \code{"max"}, or \code{"both".}}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}. }
  \item{mmax}{The partition size of the ranked observations. The default size is half the square root of the number of observations}
  \item{mmin}{The partition size of the ranked observations. The default size is half the square root of the number of observations}
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}, default value 0. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}, default value 2.}
  \item{nr.atoms}{For \code{"ADP-EQP"} and \code{"ADP-EQP-ML"} type tests, sets the number of possible split points in the data. The default value is the minimum between \eqn{n} and \eqn{60+0.5*\sqrt{n}}.}
}

\details{
  For each partition size \eqn{m= mmin,\ldots,mmax}, the function computes the scores in each of the partitions (according to score type), and aggregates all scores according to the aggregation type (see details in Heller et al. , 2014). If the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}, then the computed statistic will be in \code{statistic}, otherwise the computed statistics will be in the appropriate subset of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}. Note that if the variant is \code{"ADP"}, all partition sizes are computed together in O(N^4), so the score computational complexity is O(N^4). For \code{"DDP"} and mmax>4,the score computational complexity is O(N^4)*(mmax-mmin+1).
  
  
  For the 'sum' aggregation type (default), The test statistic is the sum of log likelihood (or Pearson Chi-square) scores, of all partitions of size \eqn{m X m} of the data, normalized by the number of partitions and the data size (thus, being an estimator of the Mutual Information). For the 'max' aggregation type, the test statistic is the maximum log likelihood (or Pearson Chi-square) score acheived by a partition of data of size \code{m}, normalized by the data size. For variant type \code{"ADP-ML"}, the statistics calculated include not only the sum over \eqn{mXm} tables (symmetric tables, same number of cells on each axis), but also assymetric tables (i.e. \eqn{mXl} tables).
  
  Variant types \code{"ADP-EQP"} and \code{"ADP-EQP-ML"}, are the computationally efficient versions of the \code{"ADP"} and \code{"ADP-ML"}. EQP type variants reduce calculation time by summing over a subset of partitions, where a split between cells may be performed only every \eqn{n/nr.atoms} observations. This allows for a complexity of O(nr.atoms^4). These variants are only available for \code{aggregation.type=='sum'} type aggregation. 
  
  For large data (n>100), it is recommended to used \code{\link{Fast.independence.test}}, which is an optimized version of the \code{\link{hhg.univariate.ind.stat}} and \code{\link{hhg.univariate.ind.combined.test}} tests.
  
}

\value{
  Returns a \code{UnivariateStatistic} class object, with the following entries:
  
  \item{statistic}{The value of the computed statistic if the score type is one of \code{"LikelihoodRatio"} or \code{"Pearson"}, and the aggregation type is one of \code{"sum"} or \code{"max"}. One of \code{sum.chisq}, \code{sum.lr},  \code{max.chisq}, and \code{max.lr}.}
  
  \item{sum.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all Pearson chi-squared statistics from all the \eqn{m X m} contingency tables considered, divided by the total number of observations.}
  
  \item{sum.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the average over all LikelihoodRatio statistics from all the \eqn{m X m} contingency tables considered, divided by the total number of observations.}
  
  \item{max.chisq}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{m X m} contingency tables considered.}
  
  \item{max.lr}{A vector of size \eqn{mmax-mmin+1}, where the \eqn{m-mmin+1} entry is the maximum over all Pearson chi-squared statistics from all the \eqn{m X m} contingency tables considered.}
  
  \item{type}{"Independence"}
  
  \item{stat.type}{"Independence-Stat"}
  
  \item{size}{The sample size}
  
  \item{score.type}{The input \code{score.type}.}
  
  \item{aggregation.type}{The input \code{aggregation.type}.}
  
  \item{mmin}{The input \code{mmin}.}
  
  \item{mmax}{The input \code{mmax}.}
  
  \item{additional}{A vector with the input \code{w.sum} and \code{w.max}.}
  
  \item{nr.atoms}{The input \code{nr.atoms}.}
 

}

\references{
  
Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, JMLR 17(29):1-54
\url{https://www.jmlr.org/papers/volume17/14-441/14-441.pdf}


Brill B. (2016) Scalable Non-Parametric Tests of Independence (master's thesis)
\url{https://tau.userservices.exlibrisgroup.com/discovery/delivery/972TAU_INST:TAU/12397000130004146?lang=he&viewerServiceCode=AlmaViewer}


}

\author{
  Barak Brill and Shachar Kaufman.
}


\examples{

\dontrun{
N = 35
data = hhg.example.datagen(N, 'Parabola')
X = data[1,]
Y = data[2,]
plot(X,Y)


#I) Computing test statistics , with default parameters(ADP statistic):

hhg.univariate.ADP.Likelihood.result = hhg.univariate.ind.stat(X,Y)

hhg.univariate.ADP.Likelihood.result

#II) Computing test statistics , with summation over Data Derived Partitions (DDP),
#using Pearson scores, and partition sizes up to 5:

hhg.univariate.DDP.Pearson.result = hhg.univariate.ind.stat(X,Y,variant = 'DDP',
  score.type = 'Pearson', mmax = 5)
hhg.univariate.DDP.Pearson.result

#III) Computing test statistics, for all M X L tables:
hhg.univariate.ADP.ML.Likelihood.result = hhg.univariate.ind.stat(X,Y,
variant='ADP-ML', mmax = 5)

hhg.univariate.ADP.ML.Likelihood.result

#IV) Computing test statistics, using efficient variants (for large data sets):
#Note : for independence testing with n>100, Fast.ADP.test is suggested
#rather than hhg.univariate.ind.stat.

N_Large = 1000
data_Large = hhg.example.datagen(N_Large, 'W')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)

hhg.univariate.ADP.EQP.Likelihood.result = hhg.univariate.ind.stat(X_Large
,Y_Large,variant='ADP-EQP', mmax = 20)

hhg.univariate.ADP.EQP.Likelihood.result

#note how only nr.atoms=76 are used - only 75 possible cell split locations are
#taken into consideration when computing the sum over all possible log likelihood scores.
#this can be changed using the nr.atoms argument:

hhg.univariate.ADP.EQP.Likelihood.result = hhg.univariate.ind.stat(X_Large,Y_Large,
variant='ADP-EQP',mmax = 20, nr.atoms =100)

hhg.univariate.ADP.EQP.Likelihood.result

#V) Computing the efficient sum over all MXL tables:

hhg.univariate.ADP.EQP.ML.Likelihood.result = hhg.univariate.ind.stat(X_Large,Y_Large,
variant='ADP-EQP-ML',mmax = 5)

hhg.univariate.ADP.EQP.ML.Likelihood.result
}

}
