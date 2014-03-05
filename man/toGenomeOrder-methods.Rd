\docType{methods}
\name{toGenomeOrder}
\alias{toGenomeOrder}
\alias{toGenomeOrder,GRanges-method}
\alias{toGenomeOrder,GenoSet-method}
\alias{toGenomeOrder,RangedData-method}
\title{Set a GRanges, GenoSet, or RangedData to genome order}
\usage{
toGenomeOrder(ds, ...)

\S4method{toGenomeOrder}{RangedData}(ds, strict = TRUE)

\S4method{toGenomeOrder}{GRanges}(ds, strict = TRUE)

\S4method{toGenomeOrder}{GenoSet}(ds, strict = TRUE)
}
\arguments{
  \item{ds}{GenoSet, GRanges, or RangedData}

  \item{strict}{logical, should chromosomes be in order
  specified by chrOrder?}
}
\value{
re-ordered ds
}
\description{
Returns a re-ordered object sorted by chromosome and start
position. If strict=TRUE, then chromosomes must be in order
specified by chrOrder. If ds is already ordered, no
re-ordering is done. Therefore, checking order with
isGenomeOrder, is unnecessary if order will be corrected if
isGenomeOrder is FALSE.
}
\details{
toGenomeOrder for GRanges differs from sort in that it
orders by chromsome and start position only, rather than
chromsome, strand, start, and width.
}
\examples{
data(genoset)
  toGenomeOrder( genoset.ds, strict=TRUE )
  toGenomeOrder( genoset.ds, strict=FALSE )
  toGenomeOrder( locData(genoset.ds) )
}
\author{
Peter M. Haverty
}
\seealso{
Other "genome ordering": \code{\link{chrOrder}};
\code{\link{isGenomeOrder}},
\code{\link{isGenomeOrder,GRanges-method}},
\code{\link{isGenomeOrder,GRanges-method}},
\code{\link{isGenomeOrder,RangedDataOrGenoSet-method}},
\code{\link{isGenomeOrder,RangedDataOrGenoSet-method}}
}

