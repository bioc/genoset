\docType{methods}
\name{isGenomeOrder}
\alias{isGenomeOrder}
\alias{isGenomeOrder,GRanges-method}
\alias{isGenomeOrder,RangedDataOrGenoSet-method}
\title{Check if a GRanges, GenoSet or RangedData is in genome order}
\usage{
isGenomeOrder(ds, ...)

\S4method{isGenomeOrder}{RangedDataOrGenoSet}(ds, strict = TRUE)

\S4method{isGenomeOrder}{GRanges}(ds, strict = TRUE)
}
\arguments{
  \item{ds}{GenoSet, GRanges, or RangedData}

  \item{strict}{logical, should space/chromosome order be
  identical to that from chrOrder?}
}
\value{
logical
}
\description{
Checks that rows in each chr are ordered by start.  If
strict=TRUE, then chromosomes must be in order specified by
chrOrder. isGenomeOrder for GRanges differs from order in
that it orders by chromsome and start position only, rather
than chromsome, strand, start, and width.
}
\examples{
data(genoset)
  isGenomeOrder( locData(genoset.ds) )
}
\author{
Peter M. Haverty
}
\seealso{
Other "genome ordering": \code{\link{chrOrder}};
\code{\link{toGenomeOrder}},
\code{\link{toGenomeOrder,GenoSet-method}},
\code{\link{toGenomeOrder,GenoSet-method}},
\code{\link{toGenomeOrder,GRanges-method}},
\code{\link{toGenomeOrder,GRanges-method}},
\code{\link{toGenomeOrder,RangedData-method}},
\code{\link{toGenomeOrder,RangedData-method}}
}

