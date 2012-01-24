\docType{methods}
\name{toGenomeOrder}
\alias{toGenomeOrder}
\alias{toGenomeOrder,GenoSet-method}
\alias{toGenomeOrder,RangedData-method}
\title{Get indices to set a RangedData or GenoSet to genome order}
\arguments{
  \item{ds}{RangedData or GenoSet}

  \item{strict}{logical, should chromosomes be in order
  specified by chrOrder?}
}
\value{
  re-ordered RangedData or GenoSet
}
\description{
  Returns a vector of indices to use in re-ordering a
  RangedData or GenoSet to genome order. If strict=TRUE,
  then chromosomes must be in order specified by chrOrder.
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
\seealso{
  Other "genome ordering": \code{\link{chrOrder}},
  \code{\link{isGenomeOrder}}
}

