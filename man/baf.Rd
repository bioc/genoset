\name{baf}
\alias{baf}
\alias{baf<-}
\title{Get or Set the baf assayData slot}
\arguments{
  \item{object}{A BAFset object}
}
\value{
  matrix
}
\description{
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

