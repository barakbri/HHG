\name{hhg.univariate.ks.nulltable}
\alias{hhg.univariate.ks.nulltable}


\title{The K-sample test null tables for all partition sizes}

\description{Functions for creating null table objects, for the omnibus distribution-free test of equality of distributions among K groups, as described in Heller et al. (2016). To be used for the p-value computation, see examples in \code{\link{hhg.univariate.ks.pvalue}}. }

\usage{
hhg.univariate.ks.nulltable(group.sizes,mmin=2,
mmax=max(4,round(min(group.sizes)/3)),variant = 'KSample-Variant',
aggregation.type='sum',score.type='LikelihoodRatio',
nr.replicates=1000,keep.simulation.data=F,
nr.atoms = nr_bins_equipartition(sum(group.sizes)),
compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001)
}

\arguments{
  \item{group.sizes}{the number of observations in each group.}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is 1/3 the number of observations in the smallest group.}
  \item{variant}{Default value is \code{'KSample-Variant'}. Setting the variant to \code{'KSample-Equipartition'} performs the K-sample tests over partitions of the data where splits between cells are at least \eqn{n/nr.atoms} apart.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default) or \code{"max"}.}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default) or \code{"Pearson"}.}
  \item{nr.replicates}{The number of permutations for the null distribution.}
  \item{keep.simulation.data}{TRUE/FALSE. If TRUE, then in addition to the sorted statistics per column, the original matrix of size nr.replicates by mmax-mmin+1 is also stored.}
  \item{nr.atoms}{For \code{variant=='KSample-Equipartition'} type tests, sets the number of possible split points in the data. The default value is the minimum between \eqn{n} and \eqn{60+0.5*\sqrt{n}}.}
  \item{compress}{TRUE or FALSE. If enabled, null tables are compressed: The lower \code{compress.p} part of the null statistics is kept at a \code{compress.p0} resolution, while the upper part is kept at a \code{compress.p1} resolution (which is finer).}.
\item{compress.p0}{Parameter for compression. This is the resolution for the lower \code{compress.p} part of the null distribution.}
\item{compress.p}{Parameter for compression. Part of the null distribution to compress.}
\item{compress.p1}{Parameter for compression. This is the resolution for the upper value of the null distribution.}
}

\details{
  In order to compute the null distributions for a test statistic  (with a specific aggregation and score type, and all partition sizes), the only necessary information is the group sizes (the test statistic is "distribution free"). The accuracy of the quantiles of the null distribution depend on the number of replicates used for constructing the null tables. The necessary accuracy depends on the threshold used for rejection of the null hypotheses.
  
  This function creates an object for efficiently storing the null distribution of the test statistics (by partition size \code{m}). Use the returned object, together with \code{\link{hhg.univariate.ks.pvalue}} to compute the P-value for the statistics computed by \code{\link{hhg.univariate.ks.stat}}
  
  Generated null tables also hold the distribution of statistics for combination types (\code{comb.type=='MinP'} and \code{comb.type=='Fisher'}), used by \code{\link{hhg.univariate.ks.combined.test}}.
  
  Variant type \code{"KSample-Equipartition"} is the computationally efficient version of the K-sample test. calculation time is reducing by aggregating over a subset of partitions, where a split between cells may be performed only every \eqn{n/nr.atoms} observations. This allows for a complexity of O(nr.atoms^2) (instead of O(n^2)). Computationly efficient versions are available for \code{aggregation.type=='sum'} and \code{aggregation.type=='max'} variants.
  
  Null tables may be compressed, using the \code{compress} argument. For each of the partition sizes (i.e. \code{m} or   \code{mXm}), the null distribution is held at a \code{compress.p0} resolution up to the \code{compress.p} percentile. Beyond that value, the distribution is held at a finer resolution defined by \code{compress.p1} (since higher values are attained when a relation exists in the data, this is required for computing the p-value accurately.)
}

\value{
  \item{m.stats}{If keep.simulation.data= TRUE, a matrix with \code{nr.replicates} rows and  \code{mmax-mmin+1} columns of null test statistics.}
  
  \item{univariate.object}{A useful format of the null tables for computing p-values efficiently.}

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
\dontrun{
#Testing for the difference between two groups, each from a normal mixture:
N0=30
N1=30

#null table for aggregation by summation: 
sum.nulltable = hhg.univariate.ks.nulltable(c(N0,N1), nr.replicates=100)
#default nr. of replicates is 1000,
#but may take several seconds. For illustration only, we use 100 replicates,
#but it is highly recommended to use at least 1000 in practice. 

#null table for aggregation by maximization: 
max.nulltable = hhg.univariate.ks.nulltable(c(N0,N1), aggregation.type = 'max',
  score.type='LikelihoodRatio', mmin = 3, mmax = 5, nr.replicates = 100)
#default nr. of replicates is 1000, but may take several seconds. For illustration only,
#we use 100 replicates, but it is highly recommended to use at least 1000 in practice.

#null tables for aggregation by summation and maximization, for large data variants:
#make sure to change mmax, such that mmax<= nr.atoms

N0_large = 5000
N1_large = 5000

Sm.EQP.null.table = hhg.univariate.ks.nulltable(c(N0_large,N1_large),
nr.replicates=200, variant = 'KSample-Equipartition', mmax = 30)
Mm.EQP.null.table = hhg.univariate.ks.nulltable(c(N0_large,N1_large),
nr.replicates=200, aggregation.type='max', variant = 'KSample-Equipartition', mmax = 30)
}

}
