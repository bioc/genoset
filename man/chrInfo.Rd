\name{chrInfo}
\alias{chrInfo}
\title{Chromosome Information}
\arguments{
  \item{object}{A GenoSet object or similar}
}
\value{
  list with start and stop position, by ordered chr
}
\description{
  Get chromosome start and stop positions
}
\details{
  Provides a matrix of start, stop and offset, in base
  numbers for each chromosome.
}
\examples{
data(genoset)
  chrInfo(genoset.ds)
  chrInfo(locData(genoset.ds))  # The same
}

