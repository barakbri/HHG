\name{hhg.univariate.ind.pvalue}
\alias{hhg.univariate.ind.pvalue}


\title{The p-value computation for the test of independence using a fixed partition size}

\description{The p-value computation for the distribution free test of independence between two univariate random variables of Heller et al. (2016) ,using a fixed partition size \code{m}.}

\usage{
hhg.univariate.ind.pvalue(statistic, NullTable, m=min(statistic$mmax,4),l=m)
}

\arguments{
  \item{statistic}{The value of the computed statistic  by the function \code{\link{hhg.univariate.ind.stat}}. The statistic object includes the score type (one of \code{"LikelihoodRatio"} or \code{"Pearson"}), and the aggregation type (one of \code{"sum"} or \code{"max"}).}
  \item{NullTable}{The null table of the statistic, which can be downloaded from the software website (http://www.math.tau.ac.il/~ruheller/Software.html) or computed by the function
  \code{\link{hhg.univariate.ind.nulltable}}. See \code{vignette('HHG')} for a method of computing null tables on multiple cores.}
  \item{m}{The partition size.}
  \item{l}{For \code{"ADP-ML"} and \code{"ADP-EQP-ML"} test variants, sets the second parameter for the partition size.}
}

\details{
  For the test statistic, the function extracts the fraction of observations in the null table that are at least as large as the test statistic, i.e. the p-value. 
  
  For \code{'DDP'} , \code{'ADP'} and \code{'ADP-EQP'} variants, the partition size is described by a single parameter \code{m} (since partition size is \eqn{m X m}). For \code{'ADP-ML'} and \code{'ADP-EQP-ML'} variants, partition sizes of data are of sizes \eqn{m X l}, allowing for assymetric tables.
}


\value{
  The p-value. 
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


#I) Computing test statistics , with default parameters:

#statistic:
hhg.univariate.ADP.Likelihood.result = hhg.univariate.ind.stat(X,Y)
hhg.univariate.ADP.Likelihood.result

#null table:
ADP.null = hhg.univariate.ind.nulltable(N)
#pvalue:
hhg.univariate.ind.pvalue(hhg.univariate.ADP.Likelihood.result, ADP.null)

#II) Computing test statistics , with summation over Data Derived Partitions (DDP),
#using Pearson scores, and partition sizes up to 5:

#statistic:
hhg.univariate.DDP.Pearson.result = hhg.univariate.ind.stat(X,Y,variant = 'DDP',
  score.type = 'Pearson', mmax = 5)
hhg.univariate.DDP.Pearson.result

#null table:
DDP.null = hhg.univariate.ind.nulltable(N,mmax = 5,variant = 'DDP',
  score.type = 'Pearson', nr.replicates = 1000)
  
#pvalue , for different partition size:
hhg.univariate.ind.pvalue(hhg.univariate.DDP.Pearson.result, DDP.null, m =2)
hhg.univariate.ind.pvalue(hhg.univariate.DDP.Pearson.result, DDP.null, m =5)


#III) computing P-value for the variants used for large N:

N_Large = 1000
data_Large = hhg.example.datagen(N_Large, 'W')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)

NullTable_ADP_EQP = hhg.univariate.ind.nulltable(N_Large, variant = 'ADP-EQP',
  nr.atoms = 30,nr.replicates=200)
NullTable_ADP_EQP_ML = hhg.univariate.ind.nulltable(N_Large,
variant = 'ADP-EQP-ML',nr.atoms = 30,nr.replicates=200)

ADP_EQP_result = hhg.univariate.ind.stat(X_Large,Y_Large,variant = 'ADP-EQP',
nr.atoms =30)
ADP_EQP_ML_result = hhg.univariate.ind.stat(X_Large,Y_Large,variant='ADP-EQP-ML',
nr.atoms = 30)

#P-value for the S_(5X5) statistic, the sum over all 5X5 partitions:
hhg.univariate.ind.pvalue(ADP_EQP_result,NullTable_ADP_EQP,m=5 )

#P-value for the S_(5X3) statistic, the sum over all 5X3 partitions:
hhg.univariate.ind.pvalue(ADP_EQP_ML_result,NullTable_ADP_EQP_ML,m=5,l=3)

}

}


