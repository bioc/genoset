\name{isGenomeOrder}
\alias{isGenomeOrder}
\title{Check if a RangedData or GenoSet is in genome order}
\usage{
  isGenomeOrder(ds, strict = FALSE)
}
\arguments{
  \item{ds}{RangedData or GenoSet}

  \item{strict}{logical, should space/chromosome order be
  identical to that from chrOrder?}
}
\value{
  logical
}
\description{
  Checks that rows in each chr are ordered by start.  If
  strict=TRUE, then chromsomes must be in order specified
  by chrOrder.
}
\examples{
data(genoset)
  isGenomeOrder( locData(genoset.ds) )
}
\author{
  Peter M. Haverty
}
\seealso{
  Other "genome ordering": \code{\link{chrOrder}},
  \code{\link{toGenomeOrder}}, \code{\link{toGenomeOrder}},
  \code{\link{toGenomeOrder}}
}

