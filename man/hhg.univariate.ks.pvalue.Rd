\name{hhg.univariate.ks.pvalue}
\alias{hhg.univariate.ks.pvalue}


\title{The p-value computation for the K-sample problem using a fixed partition size}

\description{The p-value computation for the K-sample test of Heller et al. (2014) using a fixed partition size \code{m}. }

\usage{
hhg.univariate.ks.pvalue(statistic, NullTable,m)
}

\arguments{
  \item{statistic}{The value of the computed statistic  by the function \code{\link{hhg.univariate.ks.stat}}. The statistic object includes the score type (one of \code{"LikelihoodRatio"} or \code{"Pearson"}), and the aggregation type (one of \code{"sum"} or \code{"max"}).}
  \item{NullTable}{The null table of the statistic, which can be downloaded from the software website (http://www.math.tau.ac.il/~ruheller/Software.html) or computed by the function
  
  \code{\link{hhg.univariate.ind.nulltable}}. See \code{vignette('HHG')} for a method of computing null tables on multiple cores.}
\item{m}{The partition size.}
}

\details{
 For the test statistic, the function extracts the fraction of observations in the null table that are at least as large as the test statistic, i.e. the p-value. 
}

\value{
  The p-value. 
}

\references{

Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, JMLR 17(29):1-54

Brill B. (2016) Scalable Non-Parametric Tests of Independence (master's thesis)

http://primage.tau.ac.il/libraries/theses/exeng/free/2899741.pdf

}

\author{
  Barak Brill Shachar Kaufman.
}



\examples{
\dontrun{

#Two groups, each from a different normal mixture:
N0=30
N1=30
X = c(c(rnorm(N0/2,-2,0.7),rnorm(N0/2,2,0.7)),c(rnorm(N1/2,-1.5,0.5),rnorm(N1/2,1.5,0.5)))
Y = (c(rep(0,N0),rep(1,N1)))
plot(Y,X)

#I)p-value for fixed partition size using the sum aggregation type
hhg.univariate.Sm.Likelihood.result = hhg.univariate.ks.stat(X,Y)
hhg.univariate.Sm.Likelihood.result


sum.nulltable = hhg.univariate.ks.nulltable(c(N0,N1), nr.replicates=100)
#default nr. of replicates is 1000, but may take several seconds.
#For illustration only, we use 100 replicates, but it is highly recommended
#to use at least 1000 in practice. 

#p-value for m=4 (the default):
hhg.univariate.ks.pvalue(hhg.univariate.Sm.Likelihood.result, sum.nulltable, m=4)


#p-value for m=2:
hhg.univariate.ks.pvalue(hhg.univariate.Sm.Likelihood.result, sum.nulltable, m=2)

#II) p-value for fixed partition size using the max aggregation type

hhg.univariate.Mm.likelihood.result = hhg.univariate.ks.stat(X,Y,aggregation.type = 'max')

hhg.univariate.Mm.likelihood.result

max.nulltable = hhg.univariate.ks.nulltable(c(N0,N1), aggregation.type = 'max',
  score.type='LikelihoodRatio', mmin = 3, mmax = 5, nr.replicates = 100)
#default nr. of replicates is 1000, but may take several seconds.
#For illustration only, we use 100 replicates,
#but it is highly recommended to use at least 1000 in practice.

#p-value for m=3:
hhg.univariate.ks.pvalue(hhg.univariate.Mm.likelihood.result, max.nulltable ,m = 3) 

#p-value for m=5:
hhg.univariate.ks.pvalue(hhg.univariate.Mm.likelihood.result, max.nulltable,m = 5)


#III) p-value for sum and max aggregation type, using large data variants:

#Two groups, each from a different normal mixture, total sample size is 10^4:
X_Large = c(c(rnorm(2500,-2,0.7),rnorm(2500,2,0.7)),
c(rnorm(2500,-1.5,0.5),rnorm(2500,1.5,0.5)))
Y_Large = (c(rep(0,5000),rep(1,5000)))
plot(Y_Large,X_Large)

# for these variants, make sure to change mmax so that mmax<= nr.atoms

hhg.univariate.Sm.EQP.Likelihood.result = hhg.univariate.ks.stat(X_Large,Y_Large,
variant = 'KSample-Equipartition',mmax=30)
hhg.univariate.Mm.EQP.likelihood.result = hhg.univariate.ks.stat(X_Large,Y_Large,
aggregation.type = 'max',variant = 'KSample-Equipartition',mmax=30)

N0_large = 5000
N1_large = 5000

Sm.EQP.null.table = hhg.univariate.ks.nulltable(c(N0_large,N1_large), nr.replicates=200,
variant = 'KSample-Equipartition', mmax = 30)
Mm.EQP.null.table = hhg.univariate.ks.nulltable(c(N0_large,N1_large), nr.replicates=200,
aggregation.type='max', variant = 'KSample-Equipartition', mmax = 30)

hhg.univariate.ks.pvalue(hhg.univariate.Sm.EQP.Likelihood.result, Sm.EQP.null.table, m=5)
hhg.univariate.ks.pvalue(hhg.univariate.Mm.EQP.likelihood.result, Mm.EQP.null.table, m=5)

}

}
