\name{cn}
\alias{cn}
\alias{cn<-}
\title{Get or Set the cn assayData slot}
\arguments{
  \item{object}{A BAFset object}
}
\value{
  matrix
}
\description{
  Get or Set the cn assayData slot
}
\examples{
data(genoset)
  cn(cn.ds)  # Returns assayDataElement called "cn"
  cn(cn.ds) <- cn(cn.ds) + 5
}
\author{
  Peter M. Haverty
}

