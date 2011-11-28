\name{relocateAssayData}
\alias{relocateAssayData}
\title{Update "desc" attributes for big...}
\usage{relocateAssayData(ds, new.bigmat.dir)}
\description{Update "desc" attributes for big.matrix assayDataElement to new location}
\details{Update "desc" attributes for big.matrix assayDataElement to new location. Assumes files have already
been moved on the filesystem. Assumes names of description and data files are the same.}
\value{eSet}
\author{Peter M. Haverty \email{phaverty@gene.com}}
\arguments{\item{ds}{eSet}
\item{new.bigmat.dir}{character, path to directory holding desc and data files}}
