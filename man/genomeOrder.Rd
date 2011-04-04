\name{genomeOrder}
\alias{genomeOrder}
\title{Get indices to set a RangedData or GenoSet to genome order...}
\usage{genomeOrder(ds, strict=FALSE)
}
\description{Get indices to set a RangedData or GenoSet to genome order}
\details{Returns a vector of idices to use in re-ordering a RangedData or
GenoSet to genome order. If strict=TRUE, then chromsomes must be in order specified by chrOrder.}
\value{numeric vector of indices for re-ordering}
\author{Peter M. Haverty}
\arguments{\item{ds}{RangedData or GenoSet}
\item{strict}{logical, should chromosomes be in order specified by chrOrder?}
}
\examples{data(genoset)
genomeOrder( baf.ds )
genomeOrder( baf.ds, strict=TRUE )}

