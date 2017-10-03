\name{nr_bins_equipartition}
\alias{nr_bins_equipartition}

\title{Function for selecting the default number of atoms for equipartition}

\description{Function selects the default number of atoms by sample size, being  min(n,60+0.5*floor(sqrt(n)))}

\usage{
  nr_bins_equipartition(n)
}
\arguments{
  \item{n}{Sample size}
}

\details{
Function returns the default number of atoms (bins of sample space partition) for tests performed over a large number of observations. The default number of atoms is min(n,60+0.5*floor(sqrt(n))).
}

\value{
Default number of atoms by sample size.
}


\author{
  Barak Brill
}


\examples{
  nr_bins_equipartition(100)
}
