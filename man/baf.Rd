\name{baf}
\alias{baf}
\alias{baf,BAFSet-method}
\alias{baf<-}
\alias{baf<-,BAFSet,matrix-method}
\title{Get baf data}
\arguments{
  \item{object}{A BAFset object}
}
\value{
  matrix
}
\description{
  Get or Set the baf assayData slot
}
\details{
  Get or Set the baf assayData slot
}
\examples{
data(genoset)
  baf(baf.ds)  # Returns assayDataElement called "baf"
  baf(baf.ds) <- baf2mbaf( baf(baf.ds) )
}
\author{
  Peter M. Haverty
}

