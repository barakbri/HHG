\name{hhg.univariate.ind.combined.test}
\alias{hhg.univariate.ind.combined.test}


\title{Distribution-free test of independence}

\description{Performs distribution-free tests for independence of two univariate random variables.}

\usage{
hhg.univariate.ind.combined.test(X,Y=NULL,NullTable=NULL,mmin=2,
mmax=max(floor(sqrt(length(X))/2),2),variant='ADP',aggregation.type='sum',
score.type='LikelihoodRatio', w.sum = 0, w.max = 2 ,combining.type='MinP',
nr.perm=100,nr.atoms = nr_bins_equipartition(length(X)),
compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001,
keep.simulation.data=T)
}

\arguments{
  \item{X}{a numeric vector with observed \code{X} values,  or the test statistic as output from \code{\link{hhg.univariate.ind.stat}}.}
  \item{Y}{a numeric vector with observed \code{Y} values.  Leave as Null if the input to \code{X} is the test statistic. }
  \item{NullTable}{The null table of the statistic, which can be downloaded from the software website or computed by the function \code{\link{hhg.univariate.ind.nulltable}}.}
  \item{mmin}{The minimum partition size of the ranked observations, default value is 2. Ignored if \code{NullTable} is non-null.}
  \item{mmax}{The maximum partition size of the ranked observations, default value is half the square root of the number of observations. For a max aggregation.type, this parameter cannot be more than 2 for the ADP variant and 4 for DDP variant. Ignored if \code{NullTable} is non-null.}
  \item{variant}{a character string specifying the partition type, must be one of \code{"ADP"} (default) or \code{"DDP"}, \code{"ADP-ML"}, \code{"ADP-EQP"}, \code{"ADP-EQP-ML"}. Ignored if \code{NullTable} is non-null.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"} (default), \code{"max"}, or \code{"both".}. Ignored if \code{NullTable} is non-null.}
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"} (default), \code{"Pearson"}, or \code{"both"}. Ignored if \code{NullTable} is non-null.}
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}, default value 0. Ignored if \code{NullTable} is non-null.}
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}, default value 2. Ignored if \code{NullTable} is non-null.}
  \item{combining.type}{a character string specifying the combining type, must be one of \code{"MinP"} (default), \code{"Fisher"}, or \code{"both"}.}
  \item{nr.perm}{The number of permutations for the null distribution. Ignored if \code{NullTable} is non-null.}
  \item{nr.atoms}{For \code{"ADP-EQP"} and \code{"ADP-EQP-ML"} type tests, sets the number of possible split points in the data. Ignored if \code{NullTable} is non-null. The default value is the minimum between \eqn{n} and \eqn{60+0.5*\sqrt{n}}.}
  \item{compress}{TRUE or FALSE. If enabled, null tables are compressed: The lower \code{compress.p} part of the null statistics is kept at a \code{compress.p0} resolution, while the upper part is kept at a \code{compress.p1} resolution (which is finer).}
  \item{compress.p0}{Parameter for compression. This is the resolution for the lower \code{compress.p} part of the null distribution.}
  \item{compress.p}{Parameter for compression. Part of the null distribution to compress.}
  \item{compress.p1}{Parameter for compression. This is the resolution for the upper value of the null distribution.}
  \item{keep.simulation.data}{TRUE/FALSE. If TRUE, then in addition to the sorted statistics per column, the original matrix of size nr.replicates by mmax-mmin+1 is also stored.}
}

\details{
  The test statistic and p-value of the recommended independence test between two univariate random variables in Heller et al. (2014). The default combining type in the minimum p-value, so the test statistic is the minimum p-value over the range of partition sizes m from \code{mmin} to \code{mmax}, where the p-value for a fixed partition size m is defined by the aggregation type and score type.  The combination is done over the statistics computed by \code{\link{hhg.univariate.ind.stat}}. The second type of combination method for statistics, is via a Fisher type statistic, \eqn{-\Sigma log(p_m)} (with the sum going from \eqn{mmin} to \eqn{mmax}). The returned result may include the test statistic for the \code{MinP} combination, the \code{Fisher} combination, or both (see \code{comb.type}).
  
  If the argument \code{NullTable} is supplied with a proper null table (constructed using \code{\link{hhg.univariate.ind.nulltable}}, for the data sample size), test parameters are taken from \code{NullTable}:
  
  (\code{ mmax, mmin}\code{, variant,aggregation.type,}\code{ score.type, nr.atoms} ,...).
  
  If \code{NullTable} is left \code{NULL}, a null table is generated by a call to \code{\link{hhg.univariate.ind.nulltable}} using the arguments supplied to this function. Null table is generated with \code{nr.perm} repetitions. It is stored in the returned object, under \code{generated_null_table}. When testing for multiple hypotheses, one may generate only one null table (using this function or \code{\link{hhg.univariate.ind.nulltable}}), and use it many times (thus, substantially reducing computation time). Generated null tables hold the distribution of statistics for both combination types, (\code{comb.type=='MinP'} and \code{comb.type=='Fisher'}). 
  
  If \code{X} is supplied with a statistic (\code{UnivariateStatistic} object, returned by \code{\link{hhg.univariate.ind.stat}}), X must have the statistics (by \code{m}), required by either \code{NullTable} or the user supplied arguments \code{mmin} and \code{mmax}. If \code{X} has a larger \code{mmax} arguemnt than the supplied null table object, \code{m} statistics which exceed the null table's \code{mmax} are not taken into consideration when computing the combined statistic.
  
  Variant types \code{"ADP-EQP"} and \code{"ADP-EQP-ML"}, are the computationally efficient versions of the \code{"ADP"} and \code{"ADP-ML"}. EQP type variants reduce calculation time by summing over a subset of partitions, where a split between cells may be performed only every \eqn{n/nr.atoms} observations. This allows for a complexity of O(nr.atoms^4). These variants are only available for \code{aggregation.type=='sum'} type aggregation. 
 
 Null tables may be compressed, using the \code{compress} argument. For each of the partition sizes (i.e. \code{m} or   \code{mXm}), the null distribution is held at a \code{compress.p0} resolution up to the \code{compress.p} percentile. Beyond that value, the distribution is held at a finer resolution defined by \code{compress.p1} (since higher values are attained when a relation exists in the data, this is required for computing the p-value accurately.)
  
 For large data (n>100), it is recommended to used \code{\link{Fast.independence.test}}, which is an optimized version of the \code{\link{hhg.univariate.ind.stat}} and \code{\link{hhg.univariate.ind.combined.test}} tests.
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
   
  \item{aggregation.type}{a character string specifying the aggregation type used in the , one of \code{"sum"} or \code{"max"}.}
  
  \item{score.type}{a character string specifying the score typeused in the test, one of \code{"LikelihoodRatio"} or \code{"Pearson"}. }
  
  \item{mmax}{The maximum partition size of the ranked observations used for MinP or Fisher test statistic.}
  
  \item{mmin}{The minimum partition size of the ranked observations used for MinP or Fisher test statistic.}
  
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}.}
  
  \item{nr.atoms}{The input \code{nr.atoms}.}

}

\references{

Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, JMLR 17(29):1-54

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

#I) Perform MinP & Fisher Tests - without existing null tables.
#Null tables are generated by the test function.
#using partitions sizes up to 5
results = hhg.univariate.ind.combined.test(X,Y,nr.perm = 100,mmax=5)
results


#The null table can then be accessed.
generated.null.table = results$generated_null_table


#II) Perform MinP & Fisher Tests - with existing null tables. 

#create null table for aggregation by summation (on ADP), with partitions sizes up to 5: 
ADP.null = hhg.univariate.ind.nulltable(N,mmax=5)

#create a null table, using aggregation by summation over DDP partitions,
#with partitions sizes up to 5, over Pearson scores,
#with 1000 bootstrap repetitions.
DDP.null = hhg.univariate.ind.nulltable(N,mmax = 5,variant = 'DDP',
score.type = 'Pearson', nr.replicates = 1000)

MinP.ADP.existing.null.table = hhg.univariate.ind.combined.test(X,Y, NullTable = ADP.null)
#Results 
MinP.ADP.existing.null.table

#using the other null table (DDP variant, with pearson scores):
MinP.DDP.existing.null.table = hhg.univariate.ind.combined.test(X,Y, NullTable = DDP.null)

MinP.DDP.existing.null.table

# combined test can also be performed by using the test statistic.
ADP.statistic = hhg.univariate.ind.stat(X,Y,mmax=5)
MinP.using.statistic.result = hhg.univariate.ind.combined.test(ADP.statistic,
NullTable = ADP.null)
# same result as above (as MinP.ADP.result.using.existing.null.table$MinP.pvalue)
MinP.using.statistic.result$MinP.pvalue

#III) Perform MinP & Fisher Tests - using the efficient variants for large N. 

N_Large = 1000
data_Large = hhg.example.datagen(N_Large, 'W')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)

NullTable_for_N_Large_MXM_tables = hhg.univariate.ind.nulltable(N_Large,
variant = 'ADP-EQP',nr.atoms = 30,nr.replicates=200)
NullTable_for_N_Large_MXL_tables = hhg.univariate.ind.nulltable(N_Large,
variant = 'ADP-EQP-ML',nr.atoms = 30,nr.replicates=200)

ADP_EQP_Result = hhg.univariate.ind.combined.test(X_Large,Y_Large,
NullTable_for_N_Large_MXM_tables)
ADP_EQP_ML_Result = hhg.univariate.ind.combined.test(X_Large,Y_Large,
NullTable_for_N_Large_MXL_tables)

ADP_EQP_Result
ADP_EQP_ML_Result

}


}
