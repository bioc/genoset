\name{readGenoSet}
\alias{readGenoSet}
\title{Load a GenoSet from a RData file}
\usage{
readGenoSet(path)
}
\arguments{
  \item{path}{character, path to rds or rda file}
}
\value{
GenoSet or related object (only object in RData file)
}
\description{
Given a rds file or a rda file with one object (a GenoSet
or related object), load it, and return.
}
\examples{
\dontrun{ ds = readGenoSet("/path/to/genoset.RData") }
\dontrun{ ds = readGenoSet("/path/to/genoset.rda") }
\dontrun{ ds = readGenoSet("/path/to/genoset.rds") }
}
\author{
Peter M. Haverty \email{phaverty@gene.com}
}

