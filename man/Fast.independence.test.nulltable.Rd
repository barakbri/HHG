\name{Fast.independence.test.nulltable}
\alias{Fast.independence.test.nulltable}

\title{Null tables for the atoms based distribution-free test of independence.}

\description{Functions for creating null table objects for the atoms based omnibus distribution-free test of independence between two univariate random variables.}

\usage{
Fast.independence.test.nulltable(n,mmin=2,mmax=min(10,n),
variant = 'ADP-EQP-ML',nr.atoms = min(40,n),
score.type='LikelihoodRatio',nr.perm=200,compress=T,
compress.p0=0.001, compress.p=0.99, compress.p1=0.000001)
}

\arguments{
  \item{n}{The sample size}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is the minimum between 10 and the data size.}
  \item{variant}{a character string specifying the partition type, must be one of  \code{"ADP-EQP"}, \code{"ADP-EQP-ML"} (default).}
  \item{nr.atoms}{the number of atoms (i.e., possible split points in the data). Ignored if \code{NullTable} is non-null. The default value is the minimum between \eqn{n} and \eqn{40}.}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}.}
  \item{nr.perm}{The number of permutations for the null distribution.}
  \item{compress}{a logical variable indicating whether you want to compress the null tables. If TRUE, null tables are compressed: The lower \code{compress.p} part of the null statistics is kept at a \code{compress.p0} resolution, while the upper part is kept at a \code{compress.p1} resolution (which is finer).}
  \item{compress.p0}{Parameter for compression. This is the resolution for the lower \code{compress.p} part of the null distribution.}
  \item{compress.p}{Parameter for compression. Part of the null distribution to compress.}
  \item{compress.p1}{Parameter for compression. This is the resolution for the upper value of the null distribution.}
}

\details{
 In order to compute the null distributions for a test statistic  (with a specific aggregation and score type, and all partition sizes), the only necessary information is the sample size, since the test statistic is distribution-free. The accuracy of the quantiles of the null distribution depend on the number of replicates used for constructing the null tables. The necessary accuracy depends on the threshold used for rejection of the null hypotheses.
 
  This function creates an object for efficiently storing the null distribution of the test statistics.
  Generated null tables  hold the null distribution of statistics for the two combination types, i.e. for \code{comb.type} value  (\code{'MinP'} and \code{'Fisher'}), as well as for fixed partition sizes.
  
  Variant types \code{"ADP-EQP"} and \code{"ADP-EQP-ML"}, are the atom-based generalizations of the \code{"ADP"} and \code{"ADP-ML"}. EQP type variants reduce calculation time by summing over a subset of partitions, where a split between cells may be performed only every \eqn{n/nr.atoms} observations. This allows for a complexity of O(nr.atoms^4). These variants are only available for \code{aggregation.type=='sum'} type aggregation. 
 
 Null tables may be compressed, using the \code{compress} argument. For each of the partition sizes, the null distribution is held at a \code{compress.p0} resolution up to the \code{compress.p} percentile. Beyond that value, the distribution is held at a finer resolution defined by \code{compress.p1} (since higher values are attained when a relation exists in the data, this is required for computing the p-value accurately in the tail of the null distribution).
 
}

\value{
   \item{m.stats}{If keep.simulation.data= TRUE, \code{m.stats} a matrix with \code{nr.replicates} rows and  \code{mmax-mmin+1} columns of null test statistics.}
  
  \item{univariate.object}{A useful format of the null tables for computing p-values efficiently.}
}

\references{

Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, JMLR 17(29):1-54
\url{http://www.jmlr.org/papers/volume17/14-441/14-441.pdf}


Brill B. (2016) Scalable Non-Parametric Tests of Independence (master's thesis)
\url{http://primage.tau.ac.il/libraries/theses/exeng/free/2899741.pdf}

}

\author{
  Barak Brill.
}


\examples{
\dontrun{


N_Large = 1000
data_Large = hhg.example.datagen(N_Large, 'W')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)

NullTable_for_N_Large_MXM_tables = Fast.independence.test.nulltable(N_Large,
variant = 'ADP-EQP',  nr.atoms = 30,nr.perm=200)
NullTable_for_N_Large_MXL_tables = Fast.independence.test.nulltable(N_Large,
variant = 'ADP-EQP-ML', nr.atoms = 30,nr.perm=200)

ADP_EQP_Result = Fast.independence.test(X_Large,Y_Large,
NullTable_for_N_Large_MXM_tables)
ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large,
NullTable_for_N_Large_MXL_tables)

ADP_EQP_Result
ADP_EQP_ML_Result

#null distribution depends only on data size (length(X)),
#so same null table can be used many times.
#For example, another data set:
data_Large = hhg.example.datagen(N_Large, 'Circle')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)

#you may use Fisher type scores:
ADP_EQP_Result = Fast.independence.test(X_Large,Y_Large,
NullTable_for_N_Large_MXM_tables, combining.type='Fisher')
#or both MinP and Fisher:
ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large,
NullTable_for_N_Large_MXL_tables, combining.type='Both')


ADP_EQP_Result
ADP_EQP_ML_Result

}

}
