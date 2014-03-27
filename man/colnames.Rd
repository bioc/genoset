\name{colnames}
\alias{colnames}
\alias{sampleNames}
\alias{sampleNames<-}
\title{Get colnames from a GenoSet}
\arguments{
  \item{x}{GenoSet}
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

