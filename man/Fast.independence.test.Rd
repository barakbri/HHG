\name{Fast.independence.test}
\alias{Fast.independence.test}


\title{Atoms based distribution-free test of independence}

\description{Performs the atoms based distribution-free test for independence of two univariate random variables, which is computationally efficient for large data sets (recommended for sample size greater than 100).}

\usage{
Fast.independence.test(X,Y,NullTable=NULL,mmin=2,
mmax=min(10,length(X)), variant='ADP-EQP-ML',nr.atoms = min(40,length(X)),
combining.type='MinP',score.type='LikelihoodRatio',nr.perm=200,
compress=T, compress.p0=0.001, compress.p=0.99, compress.p1=0.000001)
}

\arguments{
  \item{X}{a numeric vector with observed \code{X} values.}
  \item{Y}{a numeric vector with observed \code{Y} values.}
  \item{NullTable}{The null table of the statistic, which can be downloaded from the software website or computed by the function \code{\link{Fast.independence.test.nulltable}}.}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2. Ignored if \code{NullTable} is non-null.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is the minimum between 10 and the data size.}
  \item{variant}{a character string specifying the partition type, must be one of  \code{"ADP-EQP"} or \code{"ADP-EQP-ML"} (default). Ignored if \code{NullTable} is non-null.}
  \item{nr.atoms}{the number of atoms (i.e., possible split points in the data). Ignored if \code{NullTable} is non-null. The default value is the minimum between \eqn{n} and \eqn{40}.}
  \item{combining.type}{a character string specifying the combining type, must be one of \code{"MinP"} (default), \code{"Fisher"}, or \code{"both"}.}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}. Ignored if \code{NullTable} is non-null.}
  \item{nr.perm}{The number of permutations for the null distribution. Ignored if \code{NullTable} is non-null.}
  \item{compress}{a logical variable indicating whether you want to compress the null tables. If TRUE, null tables are compressed: The lower \code{compress.p} part of the null statistics is kept at a \code{compress.p0} resolution, while the upper part is kept at a \code{compress.p1} resolution (which is finer).}
  \item{compress.p0}{Parameter for compression. This is the resolution for the lower \code{compress.p} part of the null distribution.}
  \item{compress.p}{Parameter for compression. Part of the null distribution to compress.}
  \item{compress.p1}{Parameter for compression. This is the resolution for the upper value of the null distribution.}
}

\details{
  This function is a smart wrapper for the \code{\link{hhg.univariate.ind.combined.test}} function, with parameters optimized for a large number of observations.
  The function first calls  \code{\link{hhg.univariate.ind.stat}} to compute the vector of test statistics. Test statistics are the sum of log-likelihood
  scores, for All Derived Partitions (ADP) of the data (computed as explained in Heller et al. (2014)).
  
  For the \code{'ADP-EQP-ML'} variant, the base test statistics are:
  
    \eqn{S_{2X2}, S_{2X3} ,S_{3X2}, ... ,S_{mmax X mmax}}.
    
  For the \code{'ADP-EQP'} varint, only the sum of symmetric tables (same number of cell on both axis) is considered:
  
    \eqn{S_{2X2}, S_{3X3} ,S_{4X4}, ... ,S_{mmax X mmax}}
    
  
  Other variant are described in \code{\link{hhg.univariate.ind.combined.test}}. The above varaiants are the ones to be used for a large number of observations (n>100).
  
  Test functions are capable of handling large datasets by attempting a split only every \eqn{N/nr.atoms} observations. An atom is a sequence of observations which cannot be split when performing a partition of the data (i.e. setting \code{nr.atoms}, the number of sequences which cannot be split, sets the number of equidistant partition points). For the above variants, 'EQP' stands for equipartition over atoms. Brill (2016) suggests a minimum of 40 atoms, with an increase of up to 60 for alternatives which are more difficult to detect (on the expense of computational complexity. Algorithm complexity is O(nr.atoms^4)). Very few alternatives require over 80 atoms.
  
  The vector of \eqn{S_{mXl}} statistics is then combined according to the method suggested in Heller et al. (2014). The default combining type in the minimum p-value, so the test statistic is the minimum p-value over the range of partition sizes m from \code{mmin} to \code{mmax}, where the p-value for a fixed partition size m is defined by the aggregation type and score type.  The combination is done over the statistics computed by \code{\link{hhg.univariate.ind.stat}}. The second type of combination method for statistics, is via a Fisher type statistic, \eqn{-\Sigma log(p_m)} (with the sum going from \eqn{mmin} to \eqn{mmax}). The returned result may include the test statistic for the \code{MinP} combination, the \code{Fisher} combination, or both (see \code{comb.type}).
  
  If the argument \code{NullTable} is supplied with a proper null table (constructed using
  
  \code{\link{Fast.independence.test.nulltable}} or \code{\link{hhg.univariate.ind.nulltable}}, for the data sample size), test parameters are taken from \code{NullTable} (\code{ mmax, mmin, variant, score.type, nr.atoms} ,...).  If \code{NullTable} is left \code{NULL}, a null table is generated by a call to \code{\link{Fast.independence.test.nulltable}} using the arguments supplied to this function. Null table is generated with \code{nr.perm} repetitions. It is stored in the returned object, under \code{generated_null_table}. When testing for multiple hypotheses, one may generate only one null table (using this function or \code{\link{Fast.independence.test.nulltable}}), and use it many times (thus, substantially reducing computation time). Generated null tables hold the distribution of statistics for both combination types, (\code{comb.type=='MinP'} and \code{comb.type=='Fisher'}). 
  
 Null tables may be compressed, using the \code{compress} argument. For each of the partition sizes (i.e. \code{m} or   \code{mXm}), the null distribution is held at a \code{compress.p0} resolution up to the \code{compress.p} percentile. Beyond that value, the distribution is held at a finer resolution defined by \code{compress.p1} (since higher values are attained when a relation exists in the data, this is required for computing the p-value accurately.)
  
}

\value{
  Returns a \code{UnivariateStatistic} class object, with the following entries:
  \item{MinP}{The test statistic when the combining type is \code{"MinP"}.}
  
  \item{MinP.pvalue}{The p-value when the combining type is \code{"MinP"}.}
  
  \item{MinP.m.chosen}{The partition size m for which the p-value was the smallest.}
  
  \item{Fisher}{The test statistic when the combining type is \code{"Fisher"}.}
  
  \item{Fisher.pvalue}{The p-value when the combining type is \code{"Fisher"}.}
  
  \item{m.stats}{The statistic for each m in the range \code{mmin} to \code{mmax}.}
  
  \item{pvalues.of.single.m}{The p-values for each m in the range \code{mmin} to \code{mmax}.}
  
  \item{generated_null_table}{The null table object. Null if \code{NullTable} is non-null.}
  
  \item{stat.type}{"Independence-Combined"}
  
  \item{variant}{a character string specifying the partition type used in the test, one of \code{"ADP"} or \code{"DDP"}.}
   
  \item{aggregation.type}{ \code{"sum"} or the aggregation type used by \code{NullTable}}
  
  \item{score.type}{a character string specifying the score typeused in the test, one of \code{"LikelihoodRatio"} or \code{"Pearson"}. }
  
  \item{mmax}{The maximum partition size of the ranked observations used for MinP or Fisher test statistic.}
  
  \item{mmin}{The minimum partition size of the ranked observations used for MinP or Fisher test statistic.}
  
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}.}
  \item{nr.atoms}{The input \code{nr.atoms}.}

}

\references{

Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, JMLR 17(29):1-54
\url{https://www.jmlr.org/papers/volume17/14-441/14-441.pdf}

Brill B., Heller Y., and Heller R. (2018) Nonparametric Independence Tests and k-sample Tests for Large Sample Sizes Using Package HHG, R Journal 10.1
\url{https://journal.r-project.org/archive/2018/RJ-2018-008/RJ-2018-008.pdf}

Brill B. (2016) Scalable Non-Parametric Tests of Independence (master's thesis)
\url{http://primage.tau.ac.il/libraries/theses/exeng/free/2899741.pdf}

}

\author{
  Barak Brill
}


\examples{

\dontrun{

N_Large = 1000
data_Large = hhg.example.datagen(N_Large, 'W')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)

NullTable_for_N_Large_MXL_tables = Fast.independence.test.nulltable(N_Large,
variant = 'ADP-EQP-ML', nr.atoms = 30,nr.perm=200)

ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large,
NullTable_for_N_Large_MXL_tables)

ADP_EQP_ML_Result

#the null distribution depends only on the sample size, so the same
#null table can be used for testing different hypotheses with the same sample size.
#For example, for another data set with N_Large sample size:
data_Large = hhg.example.datagen(N_Large, 'Circle')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)

#The MinP combining method p-values may be reported:
ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large,
NullTable_for_N_Large_MXL_tables,
  combining.type='MinP')
ADP_EQP_ML_Result

#or both MinP and Fisher combining methods p-values:
ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large,
NullTable_for_N_Large_MXL_tables,
  combining.type='Both')
ADP_EQP_ML_Result


}


}
