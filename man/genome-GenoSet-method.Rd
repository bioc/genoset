\docType{methods}
\name{genome,GenoSet-method}
\alias{genome,GenoSet-method}
\alias{universe,GenoSet-method}
\title{Get and set the genome universe annotation.}
\usage{
\S4method{genome}{GenoSet}(x)
}
\arguments{
  \item{x}{GenoSet}
}
\value{
character, e.g. hg19
}
\description{
Genome version
}
\details{
The genome positions of the features in locData. The UCSC
notation (e.g. hg18, hg19, etc.) should be used.
}
\examples{
data(genoset)
  genome(genoset.ds)
  genome(genoset.ds) = "hg19"
}
\author{
Peter M. Haverty
}

