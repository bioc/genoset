\name{loadGenoSet}
\alias{loadGenoSet}
\title{Load a GenoSet from a RData file...}
\usage{loadGenoSet(path)
}
\description{Load a GenoSet from a RData file}
\details{Given a RData file with one object (a GenoSet or related object), load it, attach bigmatrix
objects as necessary, and return.}
\value{GenoSet or related object (only object in RData file)}
\author{Peter M. Haverty \email{phaverty@gene.com}}
\arguments{\item{path}{character, path to RData file}
}
\examples{\dontrun{ ds = loadGenoSet("/path/to/genoset.RData") }}

