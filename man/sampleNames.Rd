\name{sampleNames}
\alias{sampleNames}
\alias{sampleNames,GenoSet-method}
\title{Get sampleNames from a GenoSet}
\arguments{
  \item{object}{GenoSet}
}
\value{
  character vector with names of samples
}
\description{
  Get sampleNames from a GenoSet
}
\examples{
data(genoset)
  head(sampleNames(genoset.ds))
}

