\name{hhg.univariate.nulltable.from.mstats}
\alias{hhg.univariate.nulltable.from.mstats}


\title{Constructor of Distribution Free Null Table Using Existing Statistics}

\description{This function converts null test statistics for different partition sizes into the null table object necessary for the computation of p-values efficiently. }

\usage{
hhg.univariate.nulltable.from.mstats(m.stats,minm,maxm,type,variant,
size,score.type,aggregation.type, w.sum = 0, w.max = 2,
keep.simulation.data=F,nr.atoms = nr_bins_equipartition(sum(size)),
compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001)
}

\arguments{
  \item{m.stats}{ A matrix with B rows and \code{maxm} - \code{minm}+1 columns, where each row contains the test statistics for partition sizes m from \code{minm} to \code{maxm} for the sample permutation of the input sample. }
  \item{minm}{The minimum partition size of the ranked observations, default value is 2.}
  \item{maxm}{The maximum partition size of the ranked observations.}
  \item{type}{A character string specifying the test type, must be one of \code{"KSample"}, \code{"Independence"}}
  \item{variant}{A character string specifying the partition type for the test of independence,  must be one of \code{"ADP"}, \code{"DDP"}, \code{"ADP-ML"}, \code{"ADP-EQP"},\code{"ADP-EQP-ML"} if \code{type="Independence"}. If \code{type="KSample"}, must be \code{"KSample-Variant"} or \code{"KSample-Equipartition"}.}
  \item{size}{The sample size if \code{type="Independence"}, and a vector of group sizes if \code{type="KSample"}. }
  \item{score.type}{a character string specifying the score type, must be one of \code{"LikelihoodRatio"}, or \code{"Pearson"}.}
  \item{aggregation.type}{a character string specifying the aggregation type, must be one of \code{"sum"}, or \code{"max"}.}
  \item{w.sum}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Sum"} and \code{score.type="Pearson"}, default value 0. }
  \item{w.max}{The minimum number of observations in a partition, only relevant for  \code{type="Independence"}, \code{aggregation.type="Max"} and \code{score.type="Pearson"}, default value 2.}
  \item{keep.simulation.data}{TRUE/FALSE.}
\item{nr.atoms}{For \code{"ADP-EQP"}, \code{"ADP-EQP-ML"} and \code{"KSample-Equipartition"} type tests, sets the number of possible split points in the data}
\item{compress}{TRUE or FALSE. If enabled, null tables are compressed: The lower \code{compress.p} part of the null statistics is kept at a \code{compress.p0} resolution, while the upper part is kept at a \code{compress.p1} resolution (which is finer).}.
\item{compress.p0}{Parameter for compression. This is the resolution for the lower \code{compress.p} part of the null distribution.}
\item{compress.p}{Parameter for compression. Part of the null distribution to compress.}
\item{compress.p1}{Parameter for compression. This is the resolution for the upper value of the null distribution.}
}


\details{
  For finding multiple quantiles, the null table object is more efficient than a matrix of a matrix with B rows and \code{maxm} - \code{minm}+1 columns, where each row contains the test statistics for partition sizes m from \code{minm} to \code{maxm} for the sample permutation of the input sample.
  
  
  Null tables may be compressed, using the \code{compress} argument. For each of the partition sizes (i.e. \code{m} or \code{mXm}), the null distribution is held at a \code{compress.p0} resolution up to the \code{compress.p} quantile. Beyond that value, the distribution is held at a finer resolution defined by \code{compress.p1} (since higher values are attained when a relation exists in the data, this is required for computing the p-value accurately.)
  
  See \code{vignette('HHG')} for a section on how to use this function, for computing a null tables using multiple cores.
}

\value{
  \code{m.stats}{The input \code{m.stats} if \code{keep.simulation.data=TRUE}} 
  
  \code{univariate.object}{A useful format of the null tables for computing p-values efficiently.}. 

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

# 1. Downloading a lookup table from site
# download from site http://www.math.tau.ac.il/~ruheller/Software.html
####################################################################
#using an already ready null table as object (for use in test functions)
#for example, ADP likelihood ratio statistics, for the independence problem,
#for sample size n=300
load('Object-ADP-n_300.Rdata') #=>null.table

#or using a matrix of statistics generated for the null distribution,
#to create your own table.
load('ADP-nullsim-n_300.Rdata') #=>mat
null.table = hhg.univariate.nulltable.from.mstats(m.stats = mat,minm = 2,
             maxm = 5,type = 'Independence', variant = 'ADP',size = 300,
             score.type = 'LikelihoodRatio',aggregation.type = 'sum')
             
# 2. generating an independence null table using multiple cores,
#and then compiling to object.
####################################################################
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)

#generate an independence null table
nr.cores = 4 #this is computer dependent
n = 30 #size of independence problem
nr.reps.per.core = 25
mmax =5
score.type = 'LikelihoodRatio'
aggregation.type = 'sum'
variant = 'ADP'

#generating null table of size 4*25

#single core worker function
generate.null.distribution.statistic =function(){
  library(HHG)
  null.table = matrix(NA,nrow=nr.reps.per.core,ncol = mmax-1)
  for(i in 1:nr.reps.per.core){
    #note that the statistic is distribution free (based on ranks),
    #so creating a null table (for the null distribution)
    #is essentially permuting over the ranks
    statistic = hhg.univariate.ind.stat(1:n,sample(1:n),
                                        variant = variant,
                                        aggregation.type = aggregation.type,
                                        score.type = score.type,
                                        mmax = mmax)$statistic
    null.table[i,]=statistic
  }
  rownames(null.table)=NULL
  return(null.table)
}

#parallelize over cores
cl = makeCluster(nr.cores)
registerDoParallel(cl)
res = foreach(core = 1:nr.cores, .combine = rbind, .packages = 'HHG',
              .export=c('variant','aggregation.type','score.type',
              'mmax','nr.reps.per.core','n'), .options.RNG=1234) \%dorng\% 
              { generate.null.distribution.statistic() }
stopCluster(cl)

#the null table:
head(res)

#as object to be used:
null.table = hhg.univariate.nulltable.from.mstats(res,minm=2,
  maxm = mmax,type = 'Independence',
  variant = variant,size = n,score.type = score.type,
  aggregation.type = aggregation.type)

#using the null table, checking for dependence in a linear relation
x=rnorm(n)
y=x+rnorm(n)
ADP.test = hhg.univariate.ind.combined.test(x,y,null.table)
ADP.test$MinP.pvalue #pvalue


# 3. generating a k-sample null table using multiple cores
# and then compiling to object.
####################################################################

library(parallel)
library(doParallel)
library(foreach)
library(doRNG)

#generate a k sample null table
nr.cores = 4 #this is computer dependent
n1 = 25 #size of first group
n2 = 25 #size of first group
nr.reps.per.core = 25
mmax =5
score.type = 'LikelihoodRatio'
aggregation.type = 'sum'

#generating null table of size 4*25

#single core worker function
generate.null.distribution.statistic =function(){
  library(HHG)
  null.table = matrix(NA,nrow=nr.reps.per.core,ncol = mmax-1)
  for(i in 1:nr.reps.per.core){
    #note that the statistic is distribution free (based on ranks),
    #so creating a null table (for the null distribution)
    #is essentially permuting over the ranks
    statistic = hhg.univariate.ks.stat(1:(n1+n2),sample(c(rep(0,n1),rep(1,n2))),
                                        aggregation.type = aggregation.type,
                                        score.type = score.type,
                                        mmax = mmax)$statistic
    null.table[i,]=statistic
  }
  rownames(null.table)=NULL
  return(null.table)
}

#parallelize over cores
cl = makeCluster(nr.cores)
registerDoParallel(cl)
res = foreach(core = 1:nr.cores, .combine = rbind, .packages = 'HHG',
              .export=c('n1','n2','aggregation.type','score.type','mmax',
              'nr.reps.per.core'), .options.RNG=1234) \%dorng\% 
              {generate.null.distribution.statistic()}
stopCluster(cl)

#the null table:
head(res)

#as object to be used:
null.table = hhg.univariate.nulltable.from.mstats(res,minm=2,
  maxm = mmax,type = 'KSample',
  variant = 'KSample-Variant',size = c(n1,n2),score.type = score.type,
  aggregation.type = aggregation.type)

#using the null table, checking for dependence in a case of two distinct samples
x=1:(n1+n2)
y=c(rep(0,n1),rep(1,n2))
Sm.test = hhg.univariate.ks.combined.test(x,y,null.table)
Sm.test$MinP.pvalue #pvalue
}

}
