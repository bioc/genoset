\name{chrIndices-methods}
\alias{chrIndices-methods}
\alias{chrIndices}
\title{Get a matrix of first and last index of features in each chromosome...}
\description{Get a matrix of first and last index of features in each chromosome}
\details{\code{chrIndices-methods}: Sometimes it is handy to know the first and last index for each chr.
This is like chrInfo but for feature indices rather than chromosome
locations. If chr is specified, the function will return a sequence
of integers representing the row indices of features on that chromosome.

}
\value{\code{chrIndices-methods}: data.frame with "first" and "last" columns

}
\author{Peter M. Haverty}
\arguments{\item{object}{GenoSet or RangedData}
\item{chr}{character, specific chromosome name}
}
\examples{data(genoset)
chrIndices(genoset.ds)
chrIndices(locData(genoset.ds))  # The same}
\alias{chrIndices,-method}
\alias{chrIndices}

