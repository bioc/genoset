\name{genomeorder}
\alias{toGenomeOrder,RangedData-method}
\alias{toGenomeOrder}
\title{Get indices to set a RangedData or GenoSet to genome order...}
\usage{\S4method{toGenomeOrder}{RangedData}(ds, strict=FALSE)
\S4method{toGenomeOrder}{GenoSet}(ds, strict=TRUE)
}
\description{Get indices to set a RangedData or GenoSet to genome order}
\details{\code{toGenomeOrder,RangedData-method}: Returns a vector of idices to use in re-ordering a RangedData or
GenoSet to genome order. If strict=TRUE, then chromsomes must be in order specified by chrOrder.

}
\value{\code{toGenomeOrder,RangedData-method}: numeric vector of indices for re-ordering

}
\author{Peter M. Haverty}
\arguments{\item{ds}{RangedData or GenoSet}
\item{strict}{logical, should chromosomes be in order specified by chrOrder?}
}
\examples{data(genoset)
toGenomeOrder( baf.ds, strict=TRUE )
toGenomeOrder( baf.ds )
toGenomeOrder( locData(baf.ds) )}
\alias{toGenomeOrder,GenoSet-method}
\alias{toGenomeOrder}

