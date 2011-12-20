\name{readGenoSet}
\alias{readGenoSet}
\title{Load a GenoSet from a RData file}
\usage{
  readGenoSet(path)
}
\arguments{
  \item{path}{character, path to RData file}
}
\value{
  GenoSet or related object (only object in RData file)
}
\description{
  Given a RData file with one object (a GenoSet or related
  object), load it, and return.
}
\examples{
\dontrun{ ds = readGenoSet("/path/to/genoset.RData") }
}
\author{
  Peter M. Haverty \email{phaverty@gene.com}
}

