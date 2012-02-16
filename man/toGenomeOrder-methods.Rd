\docType{methods}
\name{toGenomeOrder}
\alias{toGenomeOrder}
\alias{toGenomeOrder,GenoSet-method}
\alias{toGenomeOrder,GRanges-method}
\alias{toGenomeOrder,RangedData-method}
\title{Get indices to set a RangedData or GenoSet to genome order}
\arguments{
  \item{ds}{GenoSet, GRanges, or RangedData}

  \item{strict}{logical, should chromosomes be in order
  specified by chrOrder?}
}
\value{
  re-ordered ds
}
\description{
  Returns a vector of indices to use in re-ordering a
  GRanges, RangedData or GenoSet to genome order. If
  strict=TRUE, then chromosomes must be in order specified
  by chrOrder. If ds is already ordered, no re-ordering is
  done. Therefore, checking order with isGenomeOrder, is
  unnecessary if order will be corrected if isGenomeOrder
  is FALSE.
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

