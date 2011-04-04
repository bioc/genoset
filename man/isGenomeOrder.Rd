\name{isGenomeOrder}
\alias{isGenomeOrder}
\title{Check if a RangedData or GenoSet is in genome order...}
\usage{isGenomeOrder(ds, strict=FALSE)
}
\description{Check if a RangedData or GenoSet is in genome order}
\details{Checks that rows in each chr are ordered by start.  If strict=TRUE, then chromsomes
must be in order specified by chrOrder.}
\value{logical}
\author{Peter M. Haverty}
\arguments{\item{ds}{RangedData}
\item{strict}{logical, should space/chromosome order be identical to that from chrOrder?}
}
\examples{data(genoset)
isGenomeOrder( locData(genoset.ds) )}

