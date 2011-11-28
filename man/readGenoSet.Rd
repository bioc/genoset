\name{readGenoSet}
\alias{readGenoSet}
\title{Load a GenoSet from a RData file...}
\usage{readGenoSet(path)}
\description{Load a GenoSet from a RData file}
\details{Given a RData file with one object (a GenoSet or related object), load it, attach bigmatrix
objects as necessary, and return.}
\value{GenoSet or related object (only object in RData file)}
\author{Peter M. Haverty \email{phaverty@gene.com}}
\arguments{\item{path}{character, path to RData file}}
\examples{\dontrun{ ds = readGenoSet("/path/to/genoset.RData") }}
