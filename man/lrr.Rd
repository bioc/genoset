\name{lrr}
\alias{lrr}
\alias{lrr,BAFSet-method}
\alias{lrr<-}
\alias{lrr<-,BAFSet,matrix-method}
\title{Get lrr data}
\arguments{
  \item{object}{A BAFset object}
}
\value{
  matrix
}
\description{
  Get or Set the lrr assayData slot
}
\details{
  Get or Set the lrr assayData slot
}
\examples{
data(genoset)
  lrr(baf.ds)  # Returns assayDataElement called "lrr"
  lrr(baf.ds) <- lrr(baf.ds) + 0.1
}
\author{
  Peter M. Haverty
}

