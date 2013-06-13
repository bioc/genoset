\name{colnames}
\alias{colnames}
\alias{colnames,GenoSet-method}
\alias{sampleNames}
\alias{sampleNames,GenoSet-method}
\title{Get colnames from a GenoSet}
\arguments{
  \item{object}{GenoSet}
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

