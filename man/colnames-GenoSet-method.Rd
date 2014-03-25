\docType{methods}
\name{colnames,GenoSet-method}
\alias{colnames,GenoSet-method}
\title{Get colnames from a GenoSet}
\usage{
\S4method{colnames}{GenoSet}(x)
}
\arguments{
  \item{object}{GenoSet}

  \item{x}{ A matrix-like object.  }
}
\value{
character vector with names of samples
}
\description{
Get colnames from a GenoSet
}
\examples{
data(genoset)
  head(colnames(genoset.ds))
}

