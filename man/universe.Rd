\name{universe}
\alias{universe}
\alias{universe<-}
\title{Get and set the genome universe annotation.}
\arguments{
  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{value}{character, new universe string, e.g. hg19}
}
\value{
  character, e.g. hg19

  A GenoSet object
}
\description{
  Genome universe for locData

  Set genome universe
}
\details{
  The genome positions of the features in locData. The UCSC
  notation (e.g. hg18, hg19, etc.) should be used.
}
\examples{
data(genoset)
  universe(genoset.ds)
  universe(genoset.ds) = "hg19"
}
\author{
  Peter M. Haverty

  Peter Haverty
}

