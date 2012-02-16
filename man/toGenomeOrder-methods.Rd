\docType{methods}
\name{toGenomeOrder}
\alias{toGenomeOrder}
\alias{toGenomeOrder,GenoSet-method}
\alias{toGenomeOrder,GRanges-method}
\alias{toGenomeOrder,RangedData-method}
\title{Set a GRanges, GenoSet, or RangedData to genome order}
\arguments{
  \item{ds}{GenoSet, GRanges, or RangedData}

  \item{strict}{logical, should chromosomes be in order
  specified by chrOrder?}
}
\value{
  re-ordered ds
}
\description{
  Returns a re-ordered object sorted by chromosome and
  start position. If strict=TRUE, then chromosomes must be
  in order specified by chrOrder. If ds is already ordered,
  no re-ordering is done. Therefore, checking order with
  isGenomeOrder, is unnecessary if order will be corrected
  if isGenomeOrder is FALSE.
}
\details{
  toGenomeOrder for GRanges differs from sort in that it
  orders by chromsome and start position only, rather than
  chromsome, strand, start, and width.
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
  \code{\link{isGenomeOrder}}, \code{\link{isGenomeOrder}},
  \code{\link{isGenomeOrder}}
}

