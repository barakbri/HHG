\name{hhg.univariate.ind.nulltable}
\alias{hhg.univariate.ind.nulltable}


\title{Null tables for the  distribution-free test of independence}

\description{Functions for creating null table objects for the omnibus distribution-free test of independence between two univariate random variables.}

\usage{
hhg.univariate.ind.nulltable(size,mmin=2,mmax = max(floor(sqrt(size)/2),2),
variant = 'ADP',aggregation.type = 'sum',score.type='LikelihoodRatio',
w.sum = 0, w.max = 2,nr.replicates=1000,keep.simulation.data=F,
nr.atoms = nr_bins_equipartition(size),
compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001)
}

\arguments{
  \item{size}{The sample size}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is half the square root of the number of observations.}
  \item{variant}{a character string specifying the partition type, must be one of \code{"ADP"} (default) or \code{"DDP"}, \code{"ADP-ML"}, \code{"ADP-EQP"}, \code{"ADP-EQP-ML"}.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default) or \code{"max"}.}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default) or \code{"Pearson"}. }
   \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}, default value 0. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}, default value 2.}
  \item{nr.replicates}{The number of permutations for the null distribution.}
  \item{keep.simulation.data}{TRUE/FALSE. If TRUE, then in addition to the sorted statistics per column, the original matrix of size nr.replicates by mmax-mmin+1 is also stored.}
  \item{nr.atoms}{For \code{"ADP-EQP"} and \code{"ADP-EQP-ML"} type tests, sets the number of possible split points in the data. The default value is the minimum between \eqn{n} and \eqn{60+0.5*\sqrt{n}}.}
  \item{compress}{TRUE or FALSE. If enabled, null tables are compressed: The lower \code{compress.p} part of the null statistics is kept at a \code{compress.p0} resolution, while the upper part is kept at a \code{compress.p1} resolution (which is finer).}
  \item{compress.p0}{Parameter for compression. This is the resolution for the lower \code{compress.p} part of the null distribution.}
  \item{compress.p}{Parameter for compression. Part of the null distribution to compress.}
  \item{compress.p1}{Parameter for compression. This is the resolution for the upper value of the null distribution.}
}

\details{
 In order to compute the null distributions for a test statistic  (with a specific aggregation and score type, and all partition sizes), the only necessary information is the sample size (the test statistic is "distribution free"). The accuracy of the quantiles of the null distribution depend on the number of replicates used for constructing the null tables. The necessary accuracy depends on the threshold used for rejection of the null hypotheses.
 
  This function creates an object for efficiently storing the null distribution of the test statistics (by partition size \code{m}). Use the returned object, together with \code{\link{hhg.univariate.ind.pvalue}} to compute the P-value for the statistics computed by \code{\link{hhg.univariate.ind.stat}}
  
  Generated null tables also hold the distribution of statistics for combination types (\code{comb.type=='MinP'} and \code{comb.type=='Fisher'}), used by \code{\link{hhg.univariate.ind.combined.test}}.
  
  Variant types \code{"ADP-EQP"} and \code{"ADP-EQP-ML"}, are the computationally efficient versions of the \code{"ADP"} and \code{"ADP-ML"}. EQP type variants reduce calculation time by summing over a subset of partitions, where a split between cells may be performed only every \eqn{n/nr.atoms} observations. This allows for a complexity of O(nr.atoms^4). These variants are only available for \code{aggregation.type=='sum'} type aggregation. 
 
 Null tables may be compressed, using the \code{compress} argument. For each of the partition sizes (i.e. \code{m} or   \code{mXm}), the null distribution is held at a \code{compress.p0} resolution up to the \code{compress.p} percentile. Beyond that value, the distribution is held at a finer resolution defined by \code{compress.p1} (since higher values are attained when a relation exists in the data, this is required for computing the p-value accurately in the tail of the null disribution.)
  
  For large data (n>100), it is recommended to used \code{Fast.ADP.test}, which is an optimized version of the \code{hhg.univariate.ind.stat} and \code{hhg.univariate.ind.combined} tests. Null Tables for \code{Fast.ADP.test} can be constructed using \code{Fast.ADP.nulltable}.
 
}

\value{
   \item{m.stats}{If keep.simulation.data= TRUE, \code{m.stats} a matrix with \code{nr.replicates} rows and  \code{mmax-mmin+1} columns of null test statistics.}
  
  \item{univariate.object}{A useful format of the null tables for computing p-values efficiently.}
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
#Testing for independance, sample size = 35
N=35

#null table for aggregation by summation (on ADP): 
ADP.null = hhg.univariate.ind.nulltable(N)

#create a null table, using aggregation by summation over DDP partitions,
#with partitions sizes up to 5, over Pearson scores, with 1000 bootstrap repetitions.
DDP.null = hhg.univariate.ind.nulltable(N,mmax = 5,variant = 'DDP',
score.type = 'Pearson', nr.replicates = 1000)

#Create a null table for the ADP-EQP and ADP-EQP-ML variants, 
#which are tailored for independece testing, with larger n (n>100):

N_large =200
ADP.EQP.null = hhg.univariate.ind.nulltable(N_large, variant = 'ADP-EQP',nr.atoms = 40)

#Null table for the sum of log likelihood scores over all possible M X L tables:
ADP.EQP.ML.null = hhg.univariate.ind.nulltable(N_large, variant = 'ADP-EQP-ML',
  nr.atoms = 30,nr.replicates=200)

}

}
