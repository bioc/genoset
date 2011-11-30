\name{toGenomeOrder}
\alias{toGenomeOrder}
\title{Get indices to set a RangedData or GenoSet to genome order}
\arguments{
  \item{ds}{RangedData or GenoSet}

  \item{strict}{logical, should chromosomes be in order
  specified by chrOrder?}
}
\value{
  numeric vector of indices for re-ordering
}
\description{
  Returns a vector of idices to use in re-ordering a
  RangedData or GenoSet to genome order. If strict=TRUE,
  then chromsomes must be in order specified by chrOrder.
}
\examples{
data(genoset)
  toGenomeOrder( baf.ds, strict=TRUE )
  toGenomeOrder( baf.ds )
  toGenomeOrder( locData(baf.ds) )
}
\author{
  Peter M. Haverty
}

