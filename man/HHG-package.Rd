\name{HHG-package}
\alias{HHG-package}
\alias{HHG}
\docType{package}
\title{Heller-Heller-Gorfine (HHG) Tests of Independence and Equality of Distributions}
\description{This R package implements the  permutation test of independnece between two random vectors of arbitrary dimensions, and equality of two or more multivariate distributions, introduced in Heller et al. (2013), as well as the distribution-free tests of independence and equality of distribution between two univariate random variables introduced in Heller et al. (2016).}
\details{
\tabular{ll}{
Package: \tab HHG\cr
Type: \tab Package\cr
Version: \tab 2.3.6\cr
Date: \tab 2023-12-08\cr
License: \tab GPL-2\cr
}

The package contains six major functions:

\code{\link{hhg.test}} - the permutation test for independence of two multivariate (or univariate) vectors.

\code{\link{hhg.test.k.sample}} - the permutation test for equality of a multivariate (or univariate) distribution across K groups.

\code{\link{hhg.test.2.sample}} - the permutation test for equality of a multivariate (or univariate) distribution across 2 groups.

\code{\link{hhg.univariate.ind.combined.test}} - the distribution-free test for independence of two univariate random variables (due to the computational complexity of this function, for large  sample sizes we recommend the atom based test \code{\link{Fast.independence.test}} instead).

\code{\link{hhg.univariate.ks.combined.test}} - the distribution-free test for equality of a univariate distribution across K groups. 

\code{\link{Fast.independence.test}} -  the atom based distribution-free test for independence of two univariate random variables, which is computationally efficient for large data sets (recommended for sample sizes greater than 100).

See \code{vignette('HHG')} for additional information. 
}
\author{
Barak Brill & Shachar Kaufman, based in part on an earlier implementation of the original HHG test by Ruth Heller <ruheller@post.tau.ac.il> and Yair Heller <heller.yair@gmail.com>.
Maintainer: Barak Brill <barakbri@mail.tau.ac.il>
}
\references{
Heller, R., Heller, Y., and Gorfine, M. (2013). A consistent multivariate test of association based on ranks of distances. \emph{Biometrika}, 100(2), 503-510.

Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, \emph{JMLR} 17(29):1-54
\url{https://www.jmlr.org/papers/volume17/14-441/14-441.pdf}

Brill B., Heller Y., and Heller R. (2018) Nonparametric Independence Tests and k-sample Tests for Large Sample Sizes Using Package HHG, R Journal 10.1
\url{https://journal.r-project.org/archive/2018/RJ-2018-008/RJ-2018-008.pdf}

Brill B. (2016) Scalable Non-Parametric Tests of Independence (master's thesis).
\url{https://tau.userservices.exlibrisgroup.com/discovery/delivery/972TAU_INST:TAU/12397000130004146?lang=he&viewerServiceCode=AlmaViewer}

}

\examples{

\dontrun{

# Some examples, for more see the vignette('HHG') and specific help pages

#######################################
#1. Univariate Independence Example
#######################################
#For (N<100):

N = 30
data = hhg.example.datagen(N, 'Parabola')
X = data[1,]
Y = data[2,]
plot(X,Y)

#For (N<100) , Option 1: Perform the ADP combined test
#using partitions sizes up to 4. see documentation for other parameters of the combined test 
#(it is recommended to use mmax >= 4, or the default parameter for large data sets)
combined = hhg.univariate.ind.combined.test(X,Y,nr.perm = 200,mmax=4)
combined


#For (N<100) , Option 2: Perform the hhg test:

## Compute distance matrices, on which the HHG test will be based
Dx = as.matrix(dist((X), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((Y), diag = TRUE, upper = TRUE))

hhg = hhg.test(Dx, Dy, nr.perm = 1000)

hhg

#For N>100, Fast.independence.test is the reccomended option:

N_Large = 1000
data_Large = hhg.example.datagen(N_Large, 'W')
X_Large = data_Large[1,]
Y_Large = data_Large[2,]
plot(X_Large,Y_Large)


NullTable_for_N_Large_MXL_tables = Fast.independence.test.nulltable(N_Large, variant = 'ADP-EQP-ML',
nr.atoms = 30,nr.perm=200)


ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large, NullTable_for_N_Large_MXL_tables)

ADP_EQP_ML_Result

#######################################
#2. Univariate K-Sample Example
#######################################

N0=50
N1=50
X = c(c(rnorm(N0/2,-2,0.7),rnorm(N0/2,2,0.7)),c(rnorm(N1/2,-1.5,0.5),rnorm(N1/2,1.5,0.5)))
Y = (c(rep(0,N0),rep(1,N1)))
#plot the two distributions by group index (0 or 1)
plot(Y,X)


#Option 1: Perform the distribution-free test for equality of a univariate distribution
combined.test = hhg.univariate.ks.combined.test(X,Y)
combined.test



#Option 2: Perform the permutation test for equality of distributions.


Dx = as.matrix(dist(X, diag = TRUE, upper = TRUE))

hhg = hhg.test.k.sample(Dx, Y, nr.perm = 1000)

hhg


#######################################
#3. Multivariate Independence Example:
#######################################

n=30 #number of samples
dimensions_x=5 #dimension of X matrix
dimensions_y=5 #dimension of Y matrix
X=matrix(rnorm(n*dimensions_x,mean = 0, sd = 1),nrow = n,ncol = dimensions_x) #generate noise
Y=matrix(rnorm(n*dimensions_y,mean =0, sd = 3),nrow = n,ncol = dimensions_y)

Y[,1] = Y[,1] + X[,1] + 4*(X[,1])^2 #add in the relations
Y[,2] = Y[,2] + X[,2] + 4*(X[,2])^2

#compute the distance matrix between observations.
#User may use other distance metrics.
Dx = as.matrix(dist((X)), diag = TRUE, upper = TRUE) 
Dy = as.matrix(dist((Y)), diag = TRUE, upper = TRUE)

#run test
hhg = hhg.test(Dx, Dy, nr.perm = 1000)

hhg


#######################################
#4. Multivariate K-Sample Example
#######################################

#multivariate k-sample, with k=3 groups
n=100 #number of samples in each group
x1 = matrix(rnorm(2*n),ncol = 2) #group 1
x2 = matrix(rnorm(2*n),ncol = 2) #group 2
x2[,2] = 1*x2[,1] + x2[,2]
x3 = matrix(rnorm(2*n),ncol = 2) #group 3
x3[,2] = -1*x3[,1] + x3[,2]
x= rbind(x1,x2,x3)
y=c(rep(0,n),rep(1,n),rep(2,n)) #group numbers, starting from 0 to k-1

plot(x[,1],x[,2],col = y+1,xlab = 'first component of X',ylab = 'second component of X',
     main = 'Multivariate K-Sample Example with K=3 \n Groups Marked by Different Colors')

Dx = as.matrix(dist(x, diag = TRUE, upper = TRUE)) #distance matrix

hhg = hhg.test.k.sample(Dx, y, nr.perm = 1000) 

hhg

}

}
